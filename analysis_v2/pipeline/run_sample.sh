#!/bin/bash
# analysis_v2: process ONE biological sample (manifest row) end to end.
#   download (ENA, md5) -> concat technical runs -> fastp -> STAR naive (WASP tag, RG)
#   -> strain-identity QC (stop on mismatch) -> WASP filter -> STAR FVB arm
#   -> featureCounts (naive/wasp/fvb unique; naive/fvb -M --fraction) -> cleanup
# Every step is skipped when its done-marker exists, so the script is resumable.
# Parameters: analysis_v2/params.yaml. Input: ONLY revision_v2/manifest/manifest.tsv.
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias
SID=$1
STAGE=${2:-all}            # run up to this stage and exit: download|fastp|naive|qc|wasp|fvb|counts|all (used by scheduler.sh)
THREADS=${THREADS:-8}      # STAR threads
FASTP_T=${FASTP_T:-4}; SAM_T=${SAM_T:-4}; FC_T=${FC_T:-$THREADS}   # non-STAR thread counts (do not change outputs; scheduler.sh sets 2/2/2)
FC_PAR=${FC_PAR:-1}        # featureCounts calls run concurrently (scheduler.sh sets 2)
STARLOCK=$PROJ/analysis_v2/locks/star.lock; STOPFLAG=$PROJ/analysis_v2/locks/STOP; mkdir -p $PROJ/analysis_v2/locks
STAR_MODE=${STAR_MODE:-lock}   # lock = one 8-thread STAR under flock with mapping-time GTF insertion (animals 1-17, 2026-10-14..17);
                               # shm  = request to star_server.sh (shared-memory genome, rebuilt FVB index, 2 x 4 threads; speed-up Step 2,
                               #        byte-identical outputs shown in revision_v2/speedup/06_shm_retest_b6gtf_index.md)
star_request(){ # $1 arm: queue a mapping for star_server.sh and wait; leaves <prefix>Aligned.sortedByCoord.out.bam + STAR logs like the lock mode
  local arm=$1 q=$PROJ/analysis_v2/locks/starq why; mkdir -p $q
  rm -f $q/$SID.$arm.done $q/$SID.$arm.fail
  printf 'sid=%s\narm=%s\nreads="%s"\nprefix=%s\nslog=%s\n' "$SID" "$arm" "$READS" "$D/bam/${SID}.$arm." "$D/logs/star_$arm.log" > $q/$SID.$arm.req.tmp && mv $q/$SID.$arm.req.tmp $q/$SID.$arm.req
  log "STAR $arm queued (shared-memory server)"
  while [ ! -e $q/$SID.$arm.done ] && [ ! -e $q/$SID.$arm.fail ]; do sleep 20; done
  if [ -e $q/$SID.$arm.fail ]; then why=$(cat $q/$SID.$arm.fail); rm -f $q/$SID.$arm.fail
    [ "$why" = STOP ] && { log "STAR $arm not started (STOP flag); paused"; exit 3; }; die "STAR $arm ($why)"; fi
  rm -f $q/$SID.$arm.done; return 0
}
eval "$(python3 -I $PROJ/analysis_v2/pipeline/manifest_row.py "$SID")"
D=$PROJ/analysis_v2/data/$M_tissue/$SID
mkdir -p $D/fastq $D/trim $D/bam $D/counts $D/qc $D/logs
LOG=$D/logs/run_sample.log
FASTP=$PROJ/env_v2/bin/fastp; STAR=/usr/local/bin/STAR; SAMTOOLS=/usr/local/bin/samtools; FC=/usr/local/bin/featureCounts; GATK=$PROJ/bin/gatk
GTF=$PROJ/ref/b6_gtf/GRCm39.110.gtf; REF=$PROJ/ref/b6_fasta/GRCm39.primary.fa
STAR_B6=$PROJ/ref/star_b6; STAR_FVB=$PROJ/ref/star_fvb; VCF=$PROJ/ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf
SITES=$PROJ/analysis_v2/qc/fvbnj_homalt_exonic.hetGT.sub10pct.seed20261016.vcf.gz   # fixed 10 % random subset (seed 20261016) of the 203,510 exonic hom-alt sites; speed-up Step 3
SITES_BED=$PROJ/analysis_v2/qc/fvbnj_homalt_exonic.hetGT.sub10pct.seed20261016.intervals.bed   # the same sites merged within 1 kb, restricts the GATK traversal only
log(){ echo "[$(date '+%F %T')] $SID $*" | tee -a $LOG; }
die(){ log "ERROR: $*"; exit 1; }
stop_if_flag(){ [ -e $STOPFLAG ] && { log "paused at step boundary (STOP flag present)"; exit 3; }; }   # graceful stop between steps
wait_ram(){ while [ "$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)" -lt 6 ]; do sleep 60; done; }   # non-STAR steps wait when <6 GB free
stage_done(){ [ "$STAGE" = "$1" ] && { log "stage $1 reached"; exit 0; }; }
PE=0; [ "$M_library_layout" = "PAIRED" ] && PE=1
log "start tissue=$M_tissue strain=$M_strain layout=$M_library_layout runs=$M_SRR"

# 1) download + md5 (ENA over HTTPS since 2026-10-20 11:50: the single-connection FTP retry after an md5 mismatch stalled > 30 min at 0 bytes on the F1 re-download; same host/paths, md5-checked, outputs unaffected) ---------------------------------------------------------------------------
if [ ! -e $D/fastq/done.download ]; then
  stop_if_flag
  IFS=';' read -ra URLS <<< "$M_fastq_ftp"; IFS=';' read -ra MD5S <<< "$M_fastq_md5"
  [ ${#URLS[@]} -eq ${#MD5S[@]} ] || die "ftp/md5 count mismatch"
  : > $D/fastq/md5.txt
  for i in "${!URLS[@]}"; do
    u=${URLS[$i]}; f=$(basename $u)
    # files already obtained via the NCBI SRA fallback (re-gzipped, different md5) are neither re-downloaded nor md5-checked against ENA
    if [ -s $D/fastq/ncbi_fallback.txt ] && cut -f1 $D/fastq/ncbi_fallback.txt | grep -qx "$(echo $f | sed -E 's/(_[12])?\.fastq\.gz$//')"; then continue; fi
    echo "${MD5S[$i]}  $f" >> $D/fastq/md5.txt
    [ -s $D/fastq/$f ] && continue
    log "download $f"
    # transient FTP failures: up to 3 attempts, 2 min apart (partial file resumed via .aria2 control file)
    ok=0; for try in 1 2 3; do
      /usr/bin/aria2c -q -x 4 -s 4 --retry-wait=10 --max-tries=8 -d $D/fastq -o $f "https://$u" && { ok=1; break; }
      log "download attempt $try failed $f"; sleep 120
    done
    [ $ok -eq 1 ] || die "download failed $f"
  done
  for attempt in 1 2 3; do
    [ -s $D/fastq/md5.txt ] || { echo 'all files from NCBI fallback' > $D/fastq/md5.check; break; }
    (cd $D/fastq && md5sum -c md5.txt > md5.check 2>&1) && break
    bad=$(grep FAILED $D/fastq/md5.check | cut -d: -f1)
    [ $attempt -eq 3 ] && die "md5 mismatch after ENA re-download and NCBI fallback: $bad"
    if [ $attempt -eq 1 ]; then
      log "md5 mismatch (attempt 1): $bad -> re-download from ENA with a single connection"
      for f in $bad; do rm -f $D/fastq/$f; u=$(echo "$M_fastq_ftp" | tr ';' '\n' | grep "/$f$"); /usr/bin/aria2c -q -x 1 -s 1 --retry-wait=15 --max-tries=8 -d $D/fastq -o $f "https://$u" || die "re-download failed $f"; done
    else
      # ENA copy still corrupt -> fetch the whole run from NCBI SRA (prefetch + fasterq-dump), gzip, and verify by read count
      runs=$(for f in $bad; do echo ${f%%_*} | sed 's/\.fastq\.gz//'; done | sort -u)
      log "md5 mismatch (attempt 2): $bad -> NCBI SRA fallback for run(s): $runs"
      for r in $runs; do
        rm -f $D/fastq/${r}_1.fastq.gz $D/fastq/${r}_2.fastq.gz $D/fastq/${r}.fastq.gz; mkdir -p $D/fastq/sra
        $PROJ/env_v2/bin/prefetch --max-size 200G -O $D/fastq/sra $r > $D/logs/prefetch_$r.log 2>&1 || die "prefetch $r"
        $PROJ/env_v2/bin/fasterq-dump --split-3 -e 4 -O $D/fastq/sra -t $D/fastq/sra $D/fastq/sra/$r > $D/logs/fasterq_$r.log 2>&1 || die "fasterq-dump $r"
        if [ $PE -eq 1 ]; then gzip -1 -c $D/fastq/sra/${r}_1.fastq > $D/fastq/${r}_1.fastq.gz && gzip -1 -c $D/fastq/sra/${r}_2.fastq > $D/fastq/${r}_2.fastq.gz; else gzip -1 -c $D/fastq/sra/${r}.fastq > $D/fastq/${r}.fastq.gz; fi
        n=$(( $(wc -l < $D/fastq/sra/${r}$([ $PE -eq 1 ] && echo _1).fastq) / 4 )); echo -e "$r\tNCBI_SRA_prefetch_fasterq-dump\treads=$n\t$(md5sum $D/fastq/${r}_1.fastq.gz 2>/dev/null | cut -d' ' -f1)\t$(md5sum $D/fastq/${r}_2.fastq.gz 2>/dev/null | cut -d' ' -f1)" >> $D/fastq/ncbi_fallback.txt
        rm -rf $D/fastq/sra; log "NCBI fallback $r ok (reads=$n)"
        grep -v "^[0-9a-f]*  ${r}_" $D/fastq/md5.txt > $D/fastq/md5.tmp; mv $D/fastq/md5.tmp $D/fastq/md5.txt   # ENA md5 no longer applies to these files
      done
    fi
  done
  touch $D/fastq/done.download; log "download ok ($(du -sh $D/fastq | cut -f1))"
fi
stage_done download

# 2) concatenate technical runs ------------------------------------------------------------------
if [ ! -e $D/fastq/done.concat ]; then
  IFS=',' read -ra RUNS <<< "$M_SRR"
  if [ $PE -eq 1 ]; then
    for r in 1 2; do : > $D/fastq/${SID}_$r.fastq.gz; for s in "${RUNS[@]}"; do cat $D/fastq/${s}_$r.fastq.gz >> $D/fastq/${SID}_$r.fastq.gz || die "concat $s"; done; done
  else
    : > $D/fastq/${SID}.fastq.gz; for s in "${RUNS[@]}"; do cat $D/fastq/${s}.fastq.gz >> $D/fastq/${SID}.fastq.gz || die "concat $s"; done
  fi
  for s in "${RUNS[@]}"; do rm -f $D/fastq/${s}_1.fastq.gz $D/fastq/${s}_2.fastq.gz $D/fastq/${s}.fastq.gz; done
  echo "concatenated runs: ${RUNS[*]} (cat, order as listed)" > $D/fastq/concat.txt
  touch $D/fastq/done.concat; log "concat ok (${#RUNS[@]} runs)"
fi

# 3) fastp -----------------------------------------------------------------------------------------
if [ ! -e $D/trim/done.fastp ]; then
  stop_if_flag; wait_ram; log "fastp"
  if [ $PE -eq 1 ]; then
    $FASTP -i $D/fastq/${SID}_1.fastq.gz -I $D/fastq/${SID}_2.fastq.gz -o $D/trim/${SID}_1.fq.gz -O $D/trim/${SID}_2.fq.gz --detect_adapter_for_pe -w $FASTP_T -j $D/qc/fastp.json -h $D/qc/fastp.html > $D/logs/fastp.log 2>&1 || die "fastp"
  else
    $FASTP -i $D/fastq/${SID}.fastq.gz -o $D/trim/${SID}.fq.gz -w $FASTP_T -j $D/qc/fastp.json -h $D/qc/fastp.html > $D/logs/fastp.log 2>&1 || die "fastp"
  fi
  rm -f $D/fastq/${SID}_1.fastq.gz $D/fastq/${SID}_2.fastq.gz $D/fastq/${SID}.fastq.gz
  touch $D/trim/done.fastp; log "fastp ok"
fi
stage_done fastp
if [ $PE -eq 1 ]; then READS="$D/trim/${SID}_1.fq.gz $D/trim/${SID}_2.fq.gz"; else READS="$D/trim/${SID}.fq.gz"; fi
RG="ID:$SID SM:$SID PL:ILLUMINA"

# 4) STAR naive (B6 reference, WASP tagged) ----------------------------------------------------------
if [ ! -e $D/bam/done.naive ]; then
  stop_if_flag
  if [ "$STAR_MODE" = shm ]; then star_request naive; log "STAR naive (shared-memory server) ok"
  else log "STAR naive (waiting for STAR lock)"
  ( /usr/bin/flock 9; log "STAR naive"
  $STAR --runThreadN $THREADS --genomeDir $STAR_B6 --sjdbGTFfile $GTF --readFilesIn $READS --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical \
    --waspOutputMode SAMtag --varVCFfile $VCF --outSAMattributes NH HI AS NM MD vW --outSAMattrRGline $RG \
    --outFileNamePrefix $D/bam/${SID}.naive. > $D/logs/star_naive.log 2>&1 ) 9>$STARLOCK || die "STAR naive"
  fi
  mv $D/bam/${SID}.naive.Aligned.sortedByCoord.out.bam $D/bam/${SID}.naive.bam && $SAMTOOLS index $D/bam/${SID}.naive.bam || die "index naive"
  rm -rf $D/bam/${SID}.naive._STARgenome $D/bam/${SID}.naive._STARtmp
  touch $D/bam/done.naive; log "STAR naive ok: $(grep 'Uniquely mapped reads %' $D/bam/${SID}.naive.Log.final.out | cut -d'|' -f2 | xargs)"
fi
stage_done naive

# 5) strain-identity QC (must match manifest strain) --------------------------------------------------
if [ ! -e $D/qc/done.identity ]; then
  stop_if_flag; wait_ram; log "identity QC"
  $GATK ASEReadCounter -R $REF -I $D/bam/${SID}.naive.bam -V $SITES -L $SITES_BED --min-depth-of-non-filtered-base 20 --min-mapping-quality 20 --min-base-quality 20 \
     -O $D/qc/${SID}.ase.tsv > $D/logs/ase_qc.log 2>&1 || die "ASEReadCounter"
  python3 -I $PROJ/analysis_v2/pipeline/identity_verdict.py $D/qc/${SID}.ase.tsv "$M_strain" > $D/qc/identity.tsv || { log "STRAIN IDENTITY MISMATCH: $(tail -n1 $D/qc/identity.tsv)"; touch $D/qc/MISMATCH; exit 2; }
  touch $D/qc/done.identity; log "identity ok: $(tail -n1 $D/qc/identity.tsv)"
fi
stage_done qc

# 6) WASP arm ------------------------------------------------------------------------------------------
if [ ! -e $D/bam/done.wasp ]; then
  stop_if_flag; wait_ram; log "WASP filter"
  $SAMTOOLS view -@ $SAM_T -b -e '![vW] || [vW]==1' $D/bam/${SID}.naive.bam > $D/bam/${SID}.wasp.bam && $SAMTOOLS index $D/bam/${SID}.wasp.bam || die "wasp"
  $SAMTOOLS view -c $D/bam/${SID}.naive.bam > $D/qc/naive_reads.txt; $SAMTOOLS view -c $D/bam/${SID}.wasp.bam > $D/qc/wasp_reads.txt
  $SAMTOOLS view $D/bam/${SID}.naive.bam | awk '{t="none"; for(i=12;i<=NF;i++) if($i ~ /^vW:i:/) t=substr($i,6); c[t]++} END{for(k in c) print k"\t"c[k]}' > $D/qc/vW_tag_counts.tsv
  touch $D/bam/done.wasp; log "WASP ok"
fi
stage_done wasp

# 7) STAR FVB arm -----------------------------------------------------------------------------------------
if [ ! -e $D/bam/done.fvb ]; then
  stop_if_flag
  if [ "$STAR_MODE" = shm ]; then star_request fvb; log "STAR fvb (shared-memory server, ref/star_fvb_b6gtf) ok"
  else log "STAR fvb (waiting for STAR lock)"
  ( /usr/bin/flock 9; log "STAR fvb"
  $STAR --runThreadN $THREADS --genomeDir $STAR_FVB --sjdbGTFfile $GTF --readFilesIn $READS --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMattrRGline $RG \
    --outFileNamePrefix $D/bam/${SID}.fvb. > $D/logs/star_fvb.log 2>&1 ) 9>$STARLOCK || die "STAR fvb"
  fi
  mv $D/bam/${SID}.fvb.Aligned.sortedByCoord.out.bam $D/bam/${SID}.fvb.bam && $SAMTOOLS index $D/bam/${SID}.fvb.bam || die "index fvb"
  rm -rf $D/bam/${SID}.fvb._STARgenome $D/bam/${SID}.fvb._STARtmp
  rm -f $D/trim/*.fq.gz    # trimmed FASTQ no longer needed once both arms are aligned
  touch $D/bam/done.fvb; log "STAR fvb ok: $(grep 'Uniquely mapped reads %' $D/bam/${SID}.fvb.Log.final.out | cut -d'|' -f2 | xargs)"
fi
stage_done fvb

# 8) featureCounts ------------------------------------------------------------------------------------------
if [ ! -e $D/counts/done.counts ]; then
  stop_if_flag; wait_ram; log "featureCounts"
  PF=""; [ $PE -eq 1 ] && PF="-p"
  fc_unique(){ $FC -a $GTF -o $D/counts/${SID}.$1.unique.txt -s 0 -T $FC_T $PF $D/bam/${SID}.$1.bam > $D/logs/fc_$1_unique.log 2>&1; }
  fc_multi(){  $FC -a $GTF -o $D/counts/${SID}.$1.multi.txt -M --fraction -s 0 -T $FC_T $PF $D/bam/${SID}.$1.bam > $D/logs/fc_$1_multi.log 2>&1; }
  FCRC=0; n=0
  for job in "fc_unique naive" "fc_unique wasp" "fc_unique fvb" "fc_multi naive" "fc_multi fvb"; do
    $job & n=$((n+1))
    if [ $n -ge $FC_PAR ]; then wait -n || FCRC=1; n=$((n-1)); fi
  done
  while [ $n -gt 0 ]; do wait -n || FCRC=1; n=$((n-1)); done
  [ $FCRC -eq 0 ] || die "featureCounts (see $D/logs/fc_*.log)"
  touch $D/counts/done.counts; log "featureCounts ok"
fi
stage_done counts

# 9) cleanup ---------------------------------------------------------------------------------------------------
if [ ! -e $D/done.all ]; then
  rm -f $D/trim/*.fq.gz $D/bam/${SID}.fvb.bam $D/bam/${SID}.fvb.bam.bai $D/bam/${SID}.wasp.bam $D/bam/${SID}.wasp.bam.bai
  touch $D/done.all; log "DONE (kept naive.bam, logs, counts, qc)"
fi
