#!/bin/bash
# Phase A-2: build FVB/NJ homozygous-alt, PASS, autosomal, biallelic SNP sites from MGP REL-2021,
# restricted to Ensembl 110 exons. Output: analysis_v2/qc/fvbnj_homalt_{auto,exonic}.vcf.gz
#
# NOTE (2026-10-13): in bcftools 1.21 `view -s FVB_NJ -i 'GT="AA"'` evaluates the -i expression on the
# full multi-sample record BEFORE the sample subset, so GT="AA" matched any strain and 72,197,703
# sites came through (FVB_NJ GT mostly 0/0). The sample subset and the genotype filter are therefore
# run as two separate stages. Stage 1 output is kept as fvbnj_subset_allsites.vcf.gz.
set -euo pipefail
PROJ=/usr/local/jupyter/FVB_bias
QC=$PROJ/analysis_v2/qc
MGP=$PROJ/ref/fvb_vcf/mgp_REL2021_snps.vcf.gz
GTF=$PROJ/ref/b6_gtf/GRCm39.110.gtf
mkdir -p $QC
if [ ! -s $QC/exons.GRCm39.110.bed ]; then
  echo "[$(date)] exon bed"
  awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5}' $GTF | sort -k1,1 -k2,2n -S 2G | bedtools merge > $QC/exons.GRCm39.110.bed
fi
if [ ! -s $QC/fvbnj_subset_allsites.vcf.gz.tbi ]; then
  echo "[$(date)] stage 1: subset sample FVB_NJ, PASS, autosomal SNPs (all genotypes)"
  bcftools view -s FVB_NJ -f PASS -v snps -t ^X,Y,MT $MGP -Oz -o $QC/fvbnj_subset_allsites.vcf.gz
  bcftools index -t $QC/fvbnj_subset_allsites.vcf.gz
fi
echo "[$(date)] stage 2: FVB_NJ GT hom-alt, biallelic"
bcftools view -i 'GT="AA"' -m2 -M2 -v snps $QC/fvbnj_subset_allsites.vcf.gz -Oz -o $QC/fvbnj_homalt_auto.vcf.gz
bcftools index -t $QC/fvbnj_homalt_auto.vcf.gz
echo "[$(date)] n_auto: $(bcftools index -n $QC/fvbnj_homalt_auto.vcf.gz)"
echo "[$(date)] exon intersect"
bedtools intersect -header -u -a $QC/fvbnj_homalt_auto.vcf.gz -b $QC/exons.GRCm39.110.bed | bgzip > $QC/fvbnj_homalt_exonic.vcf.gz
bcftools index -t $QC/fvbnj_homalt_exonic.vcf.gz
echo "[$(date)] n_exonic: $(bcftools index -n $QC/fvbnj_homalt_exonic.vcf.gz)"
echo "[$(date)] GT check (must be only 1/1): $(bcftools query -f '[%GT]\n' $QC/fvbnj_homalt_exonic.vcf.gz | sort | uniq -c | tr '\n' ' ')"
md5sum $MGP $GTF $QC/fvbnj_homalt_auto.vcf.gz $QC/fvbnj_homalt_exonic.vcf.gz > $QC/qc_sites.md5
echo "[$(date)] A2 sites v2 DONE"
