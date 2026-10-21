#!/bin/bash
# Phase H-1b, reproduction check: ASEReadCounter on the NEW F1 BAMs at the OLD site set (94,351 positions where MGP calls FVB_NJ
# heterozygous; revision_v2/audit/f1_old/old_sites.FVB_NJ_het.vcf.gz, built from the v1 VCF; H-1a). Same GATK parameters. The WASP BAM
# is regenerated transiently (f1_ase_rerun.sh deletes it after its own GATK run). Outputs <sid>.{naive,wasp}.oldsites.ase.tsv in
# analysis_v2/results/f1_ase/v2/. Run after f1_ase_rerun.sh has finished: bash f1_oldsites_recount.sh (4 animals in parallel)
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; MAN=$PROJ/revision_v2/manifest/manifest.tsv; OUT=$PROJ/analysis_v2/results/f1_ase/v2
GATK=$PROJ/bin/gatk; REF=$PROJ/ref/b6_fasta/GRCm39.primary.fa; SAMTOOLS=/usr/local/bin/samtools
VCF=$PROJ/revision_v2/audit/f1_old/old_sites.FVB_NJ_het.vcf.gz; BED=$PROJ/revision_v2/audit/f1_old/old_sites.merged_gap1000.bed
log(){ echo "[$(date '+%F %T')] oldsites $*"; }
one(){ s=$1
  d=$PROJ/analysis_v2/data/$(awk -F'\t' -v x=$s 'NR>1 && $1==x{print $7}' $MAN)/$s
  [ -s $d/bam/$s.naive.bam ] || { log "$s naive BAM missing"; return; }
  if [ ! -s $OUT/$s.naive.oldsites.ase.tsv ]; then
    $GATK ASEReadCounter -R $REF -I $d/bam/$s.naive.bam -V $VCF -L $BED --min-depth-of-non-filtered-base 20 --min-mapping-quality 20 --min-base-quality 20 -O $OUT/$s.naive.oldsites.ase.tsv > $d/logs/ase_oldsites_naive.log 2>&1 && log "$s naive ok" || log "$s naive FAILED"
  fi
  if [ ! -s $OUT/$s.wasp.oldsites.ase.tsv ]; then
    $SAMTOOLS view -@ 2 -b -e '![vW] || [vW]==1' $d/bam/$s.naive.bam > $d/bam/$s.wasp.tmp.bam && $SAMTOOLS index $d/bam/$s.wasp.tmp.bam || { log "$s WASP filter FAILED"; return; }
    $GATK ASEReadCounter -R $REF -I $d/bam/$s.wasp.tmp.bam -V $VCF -L $BED --min-depth-of-non-filtered-base 20 --min-mapping-quality 20 --min-base-quality 20 -O $OUT/$s.wasp.oldsites.ase.tsv > $d/logs/ase_oldsites_wasp.log 2>&1 && log "$s wasp ok" || log "$s wasp FAILED"
    rm -f $d/bam/$s.wasp.tmp.bam $d/bam/$s.wasp.tmp.bam.bai
  fi
}
export -f one log; export PROJ MAN OUT GATK REF SAMTOOLS VCF BED
# 4 animals in parallel (GATK is single-threaded; samtools limited to 2 threads each); estimated ~3 h for 12 animals
awk -F'\t' 'NR>1 && $16=="f1_ase" && $13=="Y"{print $1}' $MAN | xargs -P 4 -I{} bash -c 'one {}'
log "DONE"
