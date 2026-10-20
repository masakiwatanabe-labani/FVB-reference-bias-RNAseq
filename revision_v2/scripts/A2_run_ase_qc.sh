#!/bin/bash
# Phase A-2: run GATK ASEReadCounter at FVB/NJ hom-alt exonic SNP sites on every old naive BAM.
# Old STAR BAMs carry no @RG, which makes GATK discard all reads, so a temporary read-grouped
# copy is made per sample (samtools addreplacerg) and deleted after counting.
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias
QC=$PROJ/analysis_v2/qc
SITES=$QC/fvbnj_homalt_exonic.hetGT.vcf.gz
REF=$PROJ/ref/b6_fasta/GRCm39.primary.fa
GATK=$PROJ/bin/gatk
NPAR=${NPAR:-4}
mkdir -p $QC/ase $QC/tmp $QC/ase_logs
until [ -s $SITES.tbi ] && grep -q "A2 sites v2 DONE" $PROJ/analysis_v2/logs/A2_build_qc_sites.log; do sleep 60; done
echo "[$(date)] sites ready: $(bcftools index -n $SITES) sites"

run_one() {
  srr=$1; gse=$2
  out=$QC/ase/${srr}.tsv
  [ -s $out ] && { echo "[$(date)] $srr exists, skip"; return 0; }
  case $gse in
    GSE200632) d=$PROJ/data/GSE200632_F1/processed/$srr ;;
    GSE293440) d=$PROJ/data/GSE293440_F1/processed/$srr ;;
    *)         d=$PROJ/data/$gse/processed/$srr ;;
  esac
  bam=$d/${srr}_b6Aligned.sortedByCoord.out.bam
  [ -s $bam ] || bam=$d/${srr}_naive.bam
  [ -s $bam ] || { echo "[$(date)] $srr MISSING BAM $bam"; return 1; }
  tmp=$QC/tmp/${srr}.rg.bam
  echo "[$(date)] $srr start ($gse)"
  samtools addreplacerg -r "@RG\tID:${srr}\tSM:${srr}\tPL:ILLUMINA" -@ 2 -o $tmp $bam && samtools index $tmp || { echo "[$(date)] $srr addreplacerg FAILED"; return 1; }
  $GATK ASEReadCounter -R $REF -I $tmp -V $SITES \
     --min-depth-of-non-filtered-base 20 --min-mapping-quality 20 --min-base-quality 20 \
     -O $out > $QC/ase_logs/${srr}.log 2>&1
  rc=$?
  rm -f $tmp $tmp.bai
  echo "[$(date)] $srr done rc=$rc sites=$(( $(wc -l < $out 2>/dev/null || echo 0) - 1 ))"
}
export -f run_one; export PROJ QC SITES REF GATK
tail -n +2 $PROJ/revision_v2/audit/old_group_assignment.tsv | cut -f1,2 | xargs -P $NPAR -n 2 bash -c 'run_one "$0" "$1"'
echo "[$(date)] A2 ASE DONE"
