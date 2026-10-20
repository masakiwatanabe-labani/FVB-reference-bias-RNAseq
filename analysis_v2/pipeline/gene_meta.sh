#!/bin/bash
# Gene-level metadata for the predictor analysis and the catalog filter (Methods: SNP density = FVB/NJ SNPs in the gene body
# [Ensembl start-end inclusive] per kb from the FVB/NJ SNP VCF; gene length; biotype; presence in the g2gtools-converted FVB GTF).
set -euo pipefail
PROJ=/usr/local/jupyter/FVB_bias; OUT=$PROJ/analysis_v2/results; mkdir -p $OUT; T=$OUT/tmp_gene_meta; mkdir -p $T
GTF=$PROJ/ref/b6_gtf/GRCm39.110.gtf; FGTF=$PROJ/ref/g2g/GRCm39.110.FVB_NJ.gtf; VCF=$PROJ/ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf.gz
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10); match($9,/gene_name "[^"]+"/); nm=(RSTART>0)?substr($9,RSTART+11,RLENGTH-12):id; match($9,/gene_biotype "[^"]+"/); bt=substr($9,RSTART+14,RLENGTH-15); print $1,$4-1,$5,id,nm,bt,$7}' $GTF | sort -k1,1 -k2,2n > $T/genes.bed
bcftools query -f '%CHROM\t%POS0\t%POS\n' $VCF > $T/snps.bed
bedtools intersect -a $T/genes.bed -b $T/snps.bed -c -sorted 2>/dev/null > $T/genes_snpcount.bed || bedtools intersect -a $T/genes.bed -b $T/snps.bed -c > $T/genes_snpcount.bed
awk -F'\t' '$3=="gene"{match($9,/gene_id "[^"]+"/); print substr($9,RSTART+9,RLENGTH-10)}' $FGTF | sort -u > $T/fvb_gtf_genes.txt
python3 -I - <<PY
T="$T"; fv=set(open(f"{T}/fvb_gtf_genes.txt").read().split())
with open("$OUT/gene_meta.tsv","w") as o:
    o.write("gene_id\tgene_name\tbiotype\tchrom\tstart\tend\tstrand\tlength_bp\tn_snps\tsnp_density_per_kb\tin_fvb_gtf\n")
    for l in open(f"{T}/genes_snpcount.bed"):
        c,s,e,id,nm,bt,st,n=l.rstrip("\n").split("\t"); L=int(e)-int(s); o.write(f"{id}\t{nm}\t{bt}\t{c}\t{int(s)+1}\t{e}\t{st}\t{L}\t{n}\t{int(n)/L*1000:.4f}\t{'Y' if id in fv else 'N'}\n")
PY
wc -l $OUT/gene_meta.tsv; awk -F'\t' 'NR>1{c[$11]++} END{for(k in c) print "in_fvb_gtf="k, c[k]}' $OUT/gene_meta.tsv
rm -rf $T; echo "[$(date)] gene_meta DONE"
