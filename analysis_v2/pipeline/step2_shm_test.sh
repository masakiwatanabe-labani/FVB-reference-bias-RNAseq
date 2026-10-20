#!/bin/bash
# Speed-up Step 2 equivalence test (one animal, run once). Compares the production alignment (sequential STAR,
# built-in coordinate sort, --sjdbGTFfile at mapping time) with the shared-memory variant
# (--genomeLoad LoadAndKeep, no --sjdbGTFfile, --outSAMtype BAM Unsorted + samtools sort), two mappings concurrently
# with 4 threads each. Holds the STAR lock while running so that the production pipeline never runs a STAR beside it.
# Output: revision_v2/speedup/02_shm_equivalence.md (+ files under the scratch dir).
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; SID=${1:-hypothalamus_B6_CWHyp4}
eval "$(python3 -I $PROJ/analysis_v2/pipeline/manifest_row.py "$SID")"
D=$PROJ/analysis_v2/data/$M_tissue/$SID
W=${2:-/tmp/claude-0/-usr-local-jupyter-FVB-bias/e900eed5-e9ea-4c20-a6d8-90c422bc83a4/scratchpad/step2}; mkdir -p $W
OUT=$PROJ/revision_v2/speedup/02_shm_equivalence.md
STAR=/usr/local/bin/STAR; SAMTOOLS=/usr/local/bin/samtools; FC=/usr/local/bin/featureCounts
GTF=$PROJ/ref/b6_gtf/GRCm39.110.gtf; VCF=$PROJ/ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf
STAR_B6=$PROJ/ref/star_b6; STAR_FVB=$PROJ/ref/star_fvb; LOCK=$PROJ/analysis_v2/locks/star.lock
log(){ echo "[$(date '+%F %T')] step2test $*"; }
PE=0; [ "$M_library_layout" = "PAIRED" ] && PE=1; PF=""; [ $PE -eq 1 ] && PF="-p"
RG="ID:$SID SM:$SID PL:ILLUMINA"

# 1) keep copies of the trimmed reads (cp: the scratch dir is on another filesystem, hard links failed 2026-10-16 18:52) and of the FVB BAM before the production pipeline deletes them
until [ -e $D/trim/done.fastp ]; do sleep 30; done
if [ -s $W/r1.fq.gz ] && [ ! -s $D/trim/${SID}_1.fq.gz ] && [ ! -s $D/trim/${SID}.fq.gz ]; then [ $PE -eq 1 ] && READS="$W/r1.fq.gz $W/r2.fq.gz" || READS="$W/r1.fq.gz"; log "trimmed reads already in scratch (production copies deleted)"; elif [ $PE -eq 1 ]; then cp -f $D/trim/${SID}_1.fq.gz $W/r1.fq.gz; cp -f $D/trim/${SID}_2.fq.gz $W/r2.fq.gz; READS="$W/r1.fq.gz $W/r2.fq.gz"; else cp -f $D/trim/${SID}.fq.gz $W/r1.fq.gz; READS="$W/r1.fq.gz"; fi
log "trimmed reads copied"
until [ -e $D/bam/done.fvb ]; do sleep 20; done
{ [ -s $W/prod.fvb.bam ] && [ ! -s $D/bam/${SID}.fvb.bam ]; } || { cp -f $D/bam/${SID}.fvb.bam $W/prod.fvb.bam; log "production fvb.bam copied"; }
until [ -e $D/done.all ]; do sleep 60; done          # production counts complete -> reference values exist
cp -f $D/bam/${SID}.naive.bam $W/prod.naive.bam
log "production run complete; waiting for STAR lock"

sig(){ # content signature: header removed, records sorted by name, flag, chr, pos, CIGAR, vW tag (vW absent -> "-")
  $SAMTOOLS view $1 | awk 'BEGIN{OFS="\t"}{t="-"; for(i=12;i<=NF;i++) if($i ~ /^vW:i:/) t=$i; print $1,$2,$3,$4,$6,t}' | LC_ALL=C sort -S 2G | md5sum | cut -d' ' -f1; }
sigfull(){ $SAMTOOLS view $1 | LC_ALL=C sort -S 2G | md5sum | cut -d' ' -f1; }   # whole records (all tags), stricter
fcmd5(){ grep -v '^#' $1 | cut -f1,7 | md5sum | cut -d' ' -f1; }

run_arm(){ # $1 arm (naive|fvb)  $2 genomeDir  $3 extra args
  local arm=$1 g=$2 extra=$3
  log "$arm: loading genome into shared memory"
  local t0=$(date +%s)
  $STAR --genomeLoad LoadAndExit --genomeDir $g --outFileNamePrefix $W/$arm.load. > $W/$arm.load.log 2>&1 || { log "$arm: LoadAndExit FAILED: $(grep -h EXITING $W/$arm.load.log $W/$arm.load.Log.out 2>/dev/null | head -1)"; return 1; }   # LoadAndExit pre-loads without mapping (LoadAndKeep alone demanded readFilesIn: bug found 2026-10-16 22:11)
  local t1=$(date +%s); log "$arm: genome loaded in $((t1-t0)) s; two concurrent mappings, 4 threads each"
  for k in 1 2; do
    $STAR --runThreadN 4 --genomeLoad LoadAndKeep --genomeDir $g --readFilesIn $READS --readFilesCommand zcat \
      --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical $extra \
      --outSAMattrRGline $RG --outFileNamePrefix $W/$arm.$k. > $W/$arm.$k.log 2>&1 &
  done
  wait; local t2=$(date +%s); log "$arm: both mappings done in $((t2-t1)) s"
  for k in 1 2; do $SAMTOOLS sort -@ 4 -m 1G -o $W/$arm.$k.bam $W/$arm.$k.Aligned.out.bam 2> $W/$arm.$k.sort.log || log "$arm.$k: sort failed"; done
  local t3=$(date +%s); log "$arm: sorted (sequential, 4 threads) in $((t3-t2)) s"
  $STAR --genomeLoad Remove --genomeDir $g --outFileNamePrefix $W/$arm.remove. > $W/$arm.remove.log 2>&1
  echo "$arm load=$((t1-t0)) map2x=$((t2-t1)) sort2x=$((t3-t2))" >> $W/timings.txt
  rm -f $W/$arm.[12].Aligned.out.bam
}
compare_arm(){ # $1 arm ; production reference $W/prod.$arm.bam and counts in $D/counts
  local arm=$1
  {
  echo; echo "### $arm arm"; echo
  echo "| item | production | test run 1 | test run 2 |"; echo "|---|---|---|---|"
  echo "| STAR Log.final.out: uniquely mapped % | $(grep 'Uniquely mapped reads %' $D/bam/$SID.$arm.Log.final.out | cut -d'|' -f2 | xargs) | $(grep 'Uniquely mapped reads %' $W/$arm.1.Log.final.out | cut -d'|' -f2 | xargs) | $(grep 'Uniquely mapped reads %' $W/$arm.2.Log.final.out | cut -d'|' -f2 | xargs) |"
  echo "| md5 of records sorted by (name, flag, chr, pos, CIGAR, vW) | $(sig $W/prod.$arm.bam) | $(sig $W/$arm.1.bam) | $(sig $W/$arm.2.bam) |"
  echo "| md5 of whole records (all tags), sorted | $(sigfull $W/prod.$arm.bam) | $(sigfull $W/$arm.1.bam) | $(sigfull $W/$arm.2.bam) |"
  for k in 1 2; do $SAMTOOLS flagstat $W/$arm.$k.bam > $W/$arm.$k.flagstat; done; $SAMTOOLS flagstat $W/prod.$arm.bam > $W/prod.$arm.flagstat
  echo "| flagstat identical to production | — | $(diff -q $W/prod.$arm.flagstat $W/$arm.1.flagstat >/dev/null && echo yes || echo NO) | $(diff -q $W/prod.$arm.flagstat $W/$arm.2.flagstat >/dev/null && echo yes || echo NO) |"
  for k in 1 2; do
    $FC -a $GTF -o $W/$arm.$k.unique.txt -s 0 -T 4 $PF $W/$arm.$k.bam > $W/$arm.$k.fcu.log 2>&1
    $FC -a $GTF -o $W/$arm.$k.multi.txt -M --fraction -s 0 -T 4 $PF $W/$arm.$k.bam > $W/$arm.$k.fcm.log 2>&1
  done
  echo "| featureCounts unique (gene, count) md5 | $(fcmd5 $D/counts/$SID.$arm.unique.txt) | $(fcmd5 $W/$arm.1.unique.txt) | $(fcmd5 $W/$arm.2.unique.txt) |"
  echo "| featureCounts -M --fraction md5 | $(fcmd5 $D/counts/$SID.$arm.multi.txt) | $(fcmd5 $W/$arm.1.multi.txt) | $(fcmd5 $W/$arm.2.multi.txt) |"
  } >> $OUT
}

( /usr/bin/flock 9
  log "STAR lock acquired"
  {
  echo "# Speed-up Step 2 — shared-memory STAR equivalence test"; echo
  echo "Animal: $SID ($M_tissue, $M_strain, $M_library_layout). Generated by analysis_v2/pipeline/step2_shm_test.sh on $(date '+%F %T')."; echo
  echo "Production: \`run_sample.sh\` (STAR 8 threads, --outSAMtype BAM SortedByCoordinate, --sjdbGTFfile at mapping time, --genomeLoad NoSharedMemory)."
  echo "Test: --genomeLoad LoadAndKeep, no --sjdbGTFfile (both indices were built with the GTF, sjdbOverhang 100), --outSAMtype BAM Unsorted + samtools sort -@4 -m 1G, two mappings of the same reads concurrently with 4 threads each; all other options identical (naive: --waspOutputMode SAMtag --varVCFfile, --outSAMattributes NH HI AS NM MD vW; RG line)."
  } > $OUT
  run_arm naive $STAR_B6 "--waspOutputMode SAMtag --varVCFfile $VCF --outSAMattributes NH HI AS NM MD vW" && compare_arm naive
  run_arm fvb $STAR_FVB "" && compare_arm fvb
  {
  echo; echo "### Timings (seconds)"; echo; echo '```'; cat $W/timings.txt; echo '```'; echo
  echo "Production STAR times for this animal: naive $(awk '/STAR naive$/{a=$2} /STAR naive ok/{b=$2} END{gsub(/\]/,"",a); gsub(/\]/,"",b); print a" -> "b}' $D/logs/run_sample.log), fvb $(awk '/STAR fvb$/{a=$2} /STAR fvb ok/{b=$2} END{gsub(/\]/,"",a); gsub(/\]/,"",b); print a" -> "b}' $D/logs/run_sample.log) (8 threads, sequential, incl. built-in sort)."
  echo; echo "STAR messages (first lines of test run logs):"; echo '```'; head -n 3 $W/naive.1.log $W/fvb.1.log 2>/dev/null; grep -h "EXITING\|FATAL\|WARNING" $W/*.log | head; echo '```'
  } >> $OUT
  for g in $STAR_B6 $STAR_FVB; do $STAR --genomeLoad Remove --genomeDir $g --outFileNamePrefix $W/cleanup.$(basename $g). > /dev/null 2>&1; done   # safety: never leave a 25 GB genome in shared memory
  log "test complete; lock released"
) 9>$LOCK
