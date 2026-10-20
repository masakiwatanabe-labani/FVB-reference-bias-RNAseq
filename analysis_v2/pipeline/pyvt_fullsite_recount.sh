#!/bin/bash
# Follow-up approved 2026-10-20 (Phase C report, item 2): re-count the FULL 203,510-site set restricted to chr 3, 10 and 16 for the
# 13 PyVT animals (rule 0.4-5: identity re-check before any wording about residual FVB-derived segments in the congenic background).
# Comparators (same sites, same command) so that a site-level artefact can be told from a real FVB-allele signal: C57BL/6J animals
# of the single-end skin set (same layout as PyVT) and of the paired-end hypothalamus set. Identical GATK parameters to the
# production identity QC (run_sample.sh step 5): min depth 20, MAPQ 20, BQ 20, naive-arm BAM. Outputs next to the production files:
#   <animal>/qc/<sid>.ase.full_chr3_10_16.tsv   (never overwrites <sid>.ase.tsv; the production verdicts are untouched)
# Summary: analysis_v2/results/negative_control/pyvt_fullsite_chr3_10_16.tsv is written by followup_analysis.py.
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; MAN=$PROJ/revision_v2/manifest/manifest.tsv
GATK=$PROJ/bin/gatk; REF=$PROJ/ref/b6_fasta/GRCm39.primary.fa
SITES=$PROJ/analysis_v2/qc/fvbnj_homalt_exonic.hetGT.vcf.gz                 # full 203,510 exonic FVB/NJ hom-alt sites (GT rewritten 0/1 for ASEReadCounter)
BED=$PROJ/analysis_v2/qc/qc_sites_merged_gap1000.bed                        # all sites merged within 1 kb (restricts the traversal only)
W=$PROJ/analysis_v2/qc/tmp; mkdir -p $W
awk '$1=="3"||$1=="10"||$1=="16"' $BED > $W/chr3_10_16.intervals.bed
log(){ echo "[$(date '+%F %T')] recount $*"; }
log "intervals on chr 3/10/16: $(wc -l < $W/chr3_10_16.intervals.bed); sites: $(bcftools view -H -r 3,10,16 $SITES | wc -l)"
sel(){ awk -F'\t' -v t=$1 -v s="$2" 'NR>1 && $13=="Y" && $7==t && $6==s {print $1}' $MAN; }   # $6 strain, $7 tissue, $13 include
LIST="$(sel mammary_tumour_PyVT 'B6.FVB-Tg(MMTV-PyVT)634Mul/LellJ') $(sel skin C57BL/6J) $(sel hypothalamus C57BL/6J)"
log "animals: $(echo $LIST | wc -w) (PyVT first, then B6 skin SE, then B6 hypothalamus PE)"
run_one(){ local sid=$1; local t=$(awk -F'\t' -v s=$1 'NR>1 && $1==s{print $7}' $MAN); local d=$PROJ/analysis_v2/data/$t/$sid; local out=$d/qc/$sid.ase.full_chr3_10_16.tsv   # one local per line: under set -u a variable assigned earlier in the same `local` statement is still unbound
  [ -s $out ] && { log "$sid already done"; return 0; }
  [ -s $d/bam/$sid.naive.bam ] || { log "$sid: naive BAM missing, skipped"; return 0; }
  local t0=$(date +%s)
  $GATK ASEReadCounter -R $REF -I $d/bam/$sid.naive.bam -V $SITES -L $W/chr3_10_16.intervals.bed --min-depth-of-non-filtered-base 20 --min-mapping-quality 20 --min-base-quality 20 \
     -O $out.tmp > $d/logs/ase_full_chr3_10_16.log 2>&1 && mv $out.tmp $out && log "$sid done in $(( $(date +%s)-t0 )) s ($(($(wc -l < $out)-1)) sites with depth >= 20)" || log "$sid FAILED (see $d/logs/ase_full_chr3_10_16.log)"
}
N=0
for sid in $LIST; do
  run_one $sid & N=$((N+1))
  [ $((N % 3)) -eq 0 ] && wait
done
wait
log "ALL DONE"
