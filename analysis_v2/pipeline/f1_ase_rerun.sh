#!/bin/bash
# Phase H-1b (2026-10-20): F1 allele-specific expression recomputed at the identity-QC site set.
# The old F1 naive/WASP BAMs were deleted with the legacy data, so the 12 F1 animals (manifest role f1_ase) are re-downloaded
# and re-aligned with the production pipeline (run_sample.sh up to stage "naive": ENA download + md5 -> fastp -> STAR naive arm
# with WASP tags, 8 threads under the STAR lock). The identity-QC stage of run_sample.sh is NOT used (its verdict table expects
# an inbred strain and would stop on an F1); instead the WASP BAM is made here with the production command
# (samtools view -e '![vW] || [vW]==1', same as run_sample.sh step 6) and GATK ASEReadCounter v4.6.2.0 is run on the naive and
# the WASP BAM at the FULL 203,510-site set (FVB/NJ hom-alt exonic SNPs, GT set to 0/1; -L = sites merged within 1 kb) with
# --min-depth-of-non-filtered-base 20 --min-mapping-quality 20 --min-base-quality 20. Per-site filtering (total >= 20, ref >= 3,
# alt >= 3) and the per-sample statistics are done afterwards by f1_ase_summary_v2.py.
# Outputs: analysis_v2/results/f1_ase/v2/<sid>.{naive,wasp}.ase.tsv (per-site), <sid>.identity.tsv (verdict at the subset for the
# record), analysis_v2/data/<tissue>/<sid>/bam/<sid>.naive.bam kept; raw/trimmed FASTQ and the WASP BAM deleted after use.
# Downloads of all 12 run in the background (2 at a time) while alignment proceeds sample by sample; GATK of sample i runs in the
# background while STAR of sample i+1 is aligning. STOP flag (analysis_v2/locks/STOP) honoured at step boundaries.
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; MAN=$PROJ/revision_v2/manifest/manifest.tsv; OUT=$PROJ/analysis_v2/results/f1_ase/v2; mkdir -p $OUT
GATK=$PROJ/bin/gatk; REF=$PROJ/ref/b6_fasta/GRCm39.primary.fa; SAMTOOLS=/usr/local/bin/samtools
SITES=$PROJ/analysis_v2/qc/fvbnj_homalt_exonic.hetGT.vcf.gz; BED=$PROJ/analysis_v2/qc/qc_sites_merged_gap1000.bed; STOP=$PROJ/analysis_v2/locks/STOP
export FASTP_T=8 SAM_T=8 THREADS=8 STAR_MODE=lock
log(){ echo "[$(date '+%F %T')] f1_rerun $*"; }
LIST=$(awk -F'\t' 'NR>1 && $16=="f1_ase" && $13=="Y"{print $1}' $MAN); log "animals: $(echo $LIST | wc -w)"
ddir(){ echo $PROJ/analysis_v2/data/$(awk -F'\t' -v s=$1 'NR>1 && $1==s{print $7}' $MAN)/$1; }
# 1) downloads in the background, 2 at a time
( for s in $LIST; do echo $s; done | xargs -P 2 -I{} bash -c 'd=$(bash -c "echo $0" {}); mkdir -p '"$PROJ"'/analysis_v2/data/$(awk -F"\t" -v s={} "NR>1 && \$1==s{print \$7}" '"$MAN"')/{}/logs; bash '"$PROJ"'/analysis_v2/pipeline/run_sample.sh {} download >> '"$PROJ"'/analysis_v2/data/$(awk -F"\t" -v s={} "NR>1 && \$1==s{print \$7}" '"$MAN"')/{}/logs/worker.log 2>&1; echo "[$(date "+%F %T")] f1_rerun download {} rc=$?"' ) &
DL=$!
gatk_one(){ local sid=$1 d=$2 arm=$3
  $GATK ASEReadCounter -R $REF -I $d/bam/$sid.$arm.bam -V $SITES -L $BED --min-depth-of-non-filtered-base 20 --min-mapping-quality 20 --min-base-quality 20 \
     -O $OUT/$sid.$arm.ase.tsv.tmp > $d/logs/ase_f1_$arm.log 2>&1 && mv $OUT/$sid.$arm.ase.tsv.tmp $OUT/$sid.$arm.ase.tsv && log "$sid GATK $arm ok ($(($(wc -l < $OUT/$sid.$arm.ase.tsv)-1)) sites)" || log "$sid GATK $arm FAILED"
}
post(){ local sid=$1 d=$2   # WASP BAM + GATK both arms + verdict for the record; then delete the WASP BAM
  [ -s $OUT/$sid.naive.ase.tsv ] || gatk_one $sid $d naive
  if [ ! -s $OUT/$sid.wasp.ase.tsv ]; then
    $SAMTOOLS view -@ $SAM_T -b -e '![vW] || [vW]==1' $d/bam/$sid.naive.bam > $d/bam/$sid.wasp.bam && $SAMTOOLS index $d/bam/$sid.wasp.bam || { log "$sid WASP filter FAILED"; return 1; }
    $SAMTOOLS view -c $d/bam/$sid.naive.bam > $d/qc/naive_reads.txt; $SAMTOOLS view -c $d/bam/$sid.wasp.bam > $d/qc/wasp_reads.txt
    $SAMTOOLS view $d/bam/$sid.naive.bam | awk '{t="none"; for(i=12;i<=NF;i++) if($i ~ /^vW:i:/) t=substr($i,6); c[t]++} END{for(k in c) print k"\t"c[k]}' > $d/qc/vW_tag_counts.tsv
    gatk_one $sid $d wasp; rm -f $d/bam/$sid.wasp.bam $d/bam/$sid.wasp.bam.bai
  fi
  python3 -I $PROJ/analysis_v2/pipeline/identity_verdict.py $OUT/$sid.naive.ase.tsv "F1" > $OUT/$sid.identity.tsv 2>/dev/null; log "$sid verdict (full set): $(tail -n1 $OUT/$sid.identity.tsv | cut -f1-6)"
  touch $d/done.f1_ase
}
PIDS=""
for s in $LIST; do
  d=$(ddir $s); mkdir -p $d/logs
  [ -e $d/done.f1_ase ] && { log "$s already done"; continue; }
  [ -e $STOP ] && { log "STOP flag; halting before $s"; break; }
  until [ -e $d/fastq/done.download ]; do sleep 30; [ -e $STOP ] && break; done
  [ -e $STOP ] && { log "STOP flag; halting before $s"; break; }
  log "$s fastp + STAR naive"; t0=$(date +%s)
  bash $PROJ/analysis_v2/pipeline/run_sample.sh $s naive >> $d/logs/worker.log 2>&1; rc=$?
  [ $rc -eq 0 ] || { log "$s run_sample.sh naive rc=$rc; halting"; touch $STOP; break; }
  rm -f $d/trim/${s}_1.fq.gz $d/trim/${s}_2.fq.gz $d/trim/${s}.fq.gz
  log "$s naive BAM ready in $(( ($(date +%s)-t0)/60 )) min; GATK in background"
  ( post $s $d ) & PIDS="$PIDS $!"
done
for p in $PIDS; do wait $p; done; wait $DL 2>/dev/null
n=$(ls $OUT/*.wasp.ase.tsv 2>/dev/null | wc -l); log "DONE: $n of 12 animals have naive+WASP per-site tables in $OUT"
