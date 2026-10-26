#!/bin/bash
# Phase W-3 (2026-10-26): sex-marker read counts on the naive BAMs of the F1 (GSE200632, GSE293440) and PyVT (GSE175625) samples, to
# support the animal counts of these sets: uniquely mapped reads (MAPQ 255) in Xist (X:102503972-102526860) and in the chrY genes
# Kdm5d..Ddx3y (Y:897788-1286629), GRCm39 / Ensembl 110 coordinates from analysis_v2/results/gene_meta.tsv.
# Output: revision_v2/audit/sex_marker_counts.tsv. Usage: sex_marker_check.sh
set -euo pipefail
D=/usr/local/jupyter/FVB_bias/analysis_v2/data; OUT=/usr/local/jupyter/FVB_bias/revision_v2/audit/sex_marker_counts.tsv
{ echo -e "sample_id\ttotal_reads_idxstats\txist_reads\tchrY_gene_reads\tinferred_sex"
for b in $D/keratinocytes/*/bam/*.naive.bam $D/cerebral_cortex/*/bam/*.naive.bam $D/mammary_tumour_PyVT/*/bam/*.naive.bam; do
  tot=$(samtools idxstats $b | awk '{s+=$3}END{print s}'); x=$(samtools view -c -q 255 $b X:102503972-102526860); y=$(samtools view -c -q 255 $b Y:897788-1286629)
  sex=$(awk -v x=$x -v y=$y 'BEGIN{print (x>10*y && x>500)?"F":((y>10*x && y>500)?"M":"unclear")}')
  echo -e "$(basename $b .naive.bam)\t$tot\t$x\t$y\t$sex"; done; } > $OUT
cut -f5 $OUT | tail -n+2 | sort | uniq -c
