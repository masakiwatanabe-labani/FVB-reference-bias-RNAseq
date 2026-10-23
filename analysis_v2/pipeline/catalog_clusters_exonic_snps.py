#!/usr/bin/env python3
"""Phase P-1 (2026-10-23): (b) exonic FVB/NJ SNP count per gene, and (c) genomic clusters of catalog genes.
(b) Exons = union of the exons of all Ensembl 110 transcripts of the gene (the featureCounts meta-feature; as for Rps2 in
    rps2_exonic_variants_union.tsv); SNPs = every record of ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf.gz (the SNP VCF from which the gene-body
    SNP density of the Methods is computed). Output: analysis_v2/results/predictors/gene_exonic_snp_counts.tsv for all Ensembl 110 genes.
(c) Rule fixed before looking at the result: single-linkage of catalog genes on the same chromosome, consecutive Ensembl gene starts
    less than 1 Mb apart; clusters of >= 3 catalog genes are kept. Output: analysis_v2/results/catalog/catalog_clusters.tsv (one row per
    catalog gene; cluster = "chr:start-end Mb" of the cluster, blank if not in a cluster) and printed summary.
Deterministic recount from the VCF, the GTF and the catalog table; nothing of the pipeline is re-run. Usage: catalog_clusters_exonic_snps.py"""
import subprocess, os, tempfile
import pandas as pd
P = "/usr/local/jupyter/FVB_bias"; GTF = f"{P}/ref/b6_gtf/GRCm39.110.gtf"; VCF = f"{P}/ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf.gz"
tmp = tempfile.mkdtemp()
# (b) union of exons per gene, merged; SNP count per gene
subprocess.run(f"""awk -F'\\t' '$3=="exon"{{id=$9; sub(/.*gene_id "/,"",id); sub(/".*/,"",id); print $1"\\t"$4-1"\\t"$5"\\t"id}}' {GTF} | sort -k4,4 -k1,1 -k2,2n -S 2G > {tmp}/ex.bed""", shell=True, check=True)
rows = []; cur = None
for l in open(f"{tmp}/ex.bed"):
    c, s, e, g = l.rstrip("\n").split("\t"); s, e = int(s), int(e)
    if cur and cur[3] == g and cur[0] == c and s <= cur[2]: cur[2] = max(cur[2], e)
    else:
        if cur: rows.append(cur)
        cur = [c, s, e, g]
rows.append(cur)
with open(f"{tmp}/exu.bed", "w") as o:
    for c, s, e, g in sorted(rows, key=lambda r: (r[0], r[1])): o.write(f"{c}\t{s}\t{e}\t{g}\n")
subprocess.run(f"bcftools query -f '%CHROM\\t%POS0\\t%POS\\n' {VCF} | sort -k1,1 -k2,2n -S 2G > {tmp}/snp.bed", shell=True, check=True)
out = subprocess.run(f"bedtools intersect -c -sorted -a {tmp}/exu.bed -b {tmp}/snp.bed", shell=True, check=True, capture_output=True, text=True).stdout
d = pd.DataFrame([l.split("\t") for l in out.strip().split("\n")], columns=["chrom", "start0", "end", "gene_id", "n"])
d["n"] = d.n.astype(int); d["exon_bp"] = d.end.astype(int) - d.start0.astype(int)
g = d.groupby("gene_id").agg(n_exonic_snps=("n", "sum"), exonic_bp=("exon_bp", "sum")).reset_index()
g.to_csv(f"{P}/analysis_v2/results/predictors/gene_exonic_snp_counts.tsv", sep="\t", index=False)
print("genes:", len(g), "with zero exonic SNPs:", int((g.n_exonic_snps == 0).sum()))
# (c) clusters
cat = pd.read_csv(f"{P}/analysis_v2/results/catalog/catalog.primary.tsv", sep="\t")
meta = pd.read_csv(f"{P}/analysis_v2/results/gene_meta.tsv", sep="\t", index_col=0, dtype={"chrom": str})
cat = cat.join(meta[["chrom", "start", "end"]], on="gene_id").sort_values(["chrom", "start"])
cat["cluster"] = ""
for ch, s in cat.groupby("chrom"):
    s = s.sort_values("start"); grp = []; cur = []
    for gid_, st in zip(s.gene_id, s.start):
        if cur and st - cur[-1][1] >= 1_000_000: grp.append(cur); cur = []
        cur.append((gid_, st))
    grp.append(cur)
    for gr in grp:
        if len(gr) >= 3:
            ids = [x for x, _ in gr]; sub = cat[cat.gene_id.isin(ids)]
            lab = f"chr{ch}:{sub.start.min()/1e6:.1f}-{sub.end.max()/1e6:.1f} Mb"
            cat.loc[cat.gene_id.isin(ids), "cluster"] = lab
cat[["gene_id", "gene_name", "chrom", "start", "end", "mechanism_class", "cluster"]].to_csv(f"{P}/analysis_v2/results/catalog/catalog_clusters.tsv", sep="\t", index=False)
for lab, s in cat[cat.cluster != ""].groupby("cluster", sort=False):
    print(f"{lab}\t{len(s)}\t{', '.join(s.gene_name)}\tmulti-mapping {int((s.mechanism_class=='SNP-induced multi-mapping').sum())}, direct {int((s.mechanism_class!='SNP-induced multi-mapping').sum())}")
print("in clusters:", int((cat.cluster != "").sum()), "clusters:", cat[cat.cluster != ""].cluster.nunique())
