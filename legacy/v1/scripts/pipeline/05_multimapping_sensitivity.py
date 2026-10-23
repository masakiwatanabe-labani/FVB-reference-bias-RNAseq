#!/usr/bin/env python3
# Pipeline step 5: Multi-mapping sensitivity analysis (Supplementary Note 2)
#
# Runs featureCounts -M --fraction on naive and FVB arm BAMs.
# Computes Δ = bias_multi - bias_unique and classifies catalog genes by mechanism.
#
# Input:  $FVB_PROJ/results/strain_comp/*.counts.txt  (primary counts from step 2)
#         BAM paths extracted from featureCounts command-comment headers
# Output: $FVB_PROJ/results/paralog_analysis/
#
# Usage:
#   export FVB_PROJ=/path/to/working/directory
#   python 05_multimapping_sensitivity.py
#
"""
Paralog / multi-copy bias analysis
===================================
Part B: unique vs multi-mapping featureCounts comparison

Strategy:
  - Existing counts (no -M flag) = unique-only baseline → bias_unique (from all_datasets_bias.tsv)
  - New counts with -M --fraction  = multi-mapping included → bias_multi
  - Delta = bias_multi − bias_unique
  - Large positive Δ: multi-mappers amplify FVB undercount → paralog-driven component
  - Δ ≈ 0: bias unchanged by multi-mapper inclusion → SNP-mediated only

Output: $PROJ/results/paralog_analysis/
  bias_unique_vs_multi.tsv   — per-gene per-dataset bias_unique, bias_multi, delta
  catalog_mechanism.tsv      — catalog genes classified by mean |Δ| and SNP density
"""

import os
import pathlib, re, subprocess, sys
import numpy as np
import pandas as pd

PROJ = os.environ.get("FVB_PROJ", str(pathlib.Path(__file__).parent.parent.parent.resolve()))
COUNTS = f"{PROJ}/results/strain_comp"
STAGE6 = f"{PROJ}/results/stage6"
OUT    = f"{PROJ}/results/paralog_analysis"
GTF    = f"{PROJ}/ref/b6_gtf/GRCm39.110.gtf"
os.makedirs(OUT, exist_ok=True)

PSEUDO  = 0.5
MIN_CPM = 1.0

TISSUES = {
    "GSE123875": "adipose",
    "GSE123893": "liver",
    "GSE123894": "kidney",
    "GSE135230": "heart",
    "GSE175625": "mammary_PyVT",
}
IS_PE = {"GSE123875": True, "GSE123893": True, "GSE123894": True,
         "GSE135230": False, "GSE175625": False}

FVB_SRRS = {
    "GSE123875": ["SRR8321484","SRR8321485","SRR8321486",
                  "SRR8321487","SRR8321488","SRR8321489"],
    "GSE123893": ["SRR8321791","SRR8321792","SRR8321793",
                  "SRR8321794","SRR8321795","SRR8321796"],
    "GSE123894": ["SRR8322104","SRR8322105","SRR8322106",
                  "SRR8322107","SRR8322108","SRR8322109"],
    "GSE135230": ["SRR9879444","SRR9879445","SRR9879446","SRR9879447",
                  "SRR9879448","SRR9879449","SRR9879450","SRR9879451",
                  "SRR9879452","SRR9879453"],
    "GSE175625": ["SRR14664845","SRR14664846","SRR14664847","SRR14664848",
                  "SRR14664849","SRR14664850","SRR14664851"],
}
B6_SRRS = {
    "GSE123875": ["SRR8321439","SRR8321440","SRR8321441",
                  "SRR8321442","SRR8321443","SRR8321444"],
    "GSE123893": ["SRR8321732","SRR8321733","SRR8321734",
                  "SRR8321735","SRR8321736","SRR8321737"],
    "GSE123894": ["SRR8322030","SRR8322031","SRR8322032",
                  "SRR8322033","SRR8322034","SRR8322035"],
    "GSE135230": ["SRR9879454","SRR9879455","SRR9879456","SRR9879457",
                  "SRR9879458","SRR9879459","SRR9879460","SRR9879461",
                  "SRR9879462","SRR9879463","SRR9879464","SRR9879465"],
    "GSE175625": ["SRR14664793","SRR14664794","SRR14664795","SRR14664796",
                  "SRR14664797","SRR14664798"],
}


# ── helpers ────────────────────────────────────────────────────────────────────

def extract_bams_from_header(counts_path):
    """Read the featureCounts command header line to extract BAM paths."""
    bams = []
    with open(counts_path) as fh:
        for line in fh:
            if line.startswith("#"):
                # BAMs are quoted paths ending in .bam
                bams += re.findall(r'"([^"]+\.bam)"', line)
            else:
                break
    return bams


def load_counts(path):
    df = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    df = df.drop(columns=["Chr","Start","End","Strand","Length"], errors="ignore")
    rename = {}
    for col in df.columns:
        m = re.search(r'(SRR\d+)', os.path.basename(col))
        if m:
            rename[col] = m.group(1)
    df.rename(columns=rename, inplace=True)
    return df.astype(float)


def cpm(mat):
    lib = mat.sum(axis=0)
    return mat.divide(lib, axis=1) * 1e6


def mean_log2fc(fvb_cpm, b6_cpm):
    mf = fvb_cpm.mean(axis=1) + PSEUDO
    mb = b6_cpm.mean(axis=1)  + PSEUDO
    return np.log2(mf / mb)


def run_featurecounts(out_path, bams, is_pe, threads=8):
    if os.path.exists(out_path):
        print(f"  [skip] {os.path.basename(out_path)} already exists")
        return
    cmd = ["featureCounts", "-T", str(threads),
           "-a", GTF,
           "-o", out_path,
           "-s", "0",
           "-M", "--fraction"]
    if is_pe:
        cmd.append("-p")
    cmd += bams
    print(f"  Running: {' '.join(cmd[:8])} ... [{len(bams)} BAMs]")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR:\n{result.stderr[-2000:]}", file=sys.stderr)
        sys.exit(1)
    print(f"  Done: {os.path.basename(out_path)}")


# ── Step 1: run featureCounts -M --fraction for all datasets ──────────────────

print("\n=== Step 1: featureCounts -M --fraction runs ===")

for gse in ["GSE123875","GSE123893","GSE123894","GSE135230","GSE175625"]:
    print(f"\n{gse} ({TISSUES[gse]})")
    pe = IS_PE[gse]

    # naive_multi: same BAMs as naive.counts.txt
    naive_bams = extract_bams_from_header(f"{COUNTS}/{gse}.naive.counts.txt")
    run_featurecounts(f"{OUT}/{gse}.naive_multi.counts.txt", naive_bams, pe)

    # fvb_multi: same BAMs as fvb_b6gtf.counts.txt
    fvb_bams = extract_bams_from_header(f"{COUNTS}/{gse}.fvb_b6gtf.counts.txt")
    run_featurecounts(f"{OUT}/{gse}.fvb_multi.counts.txt", fvb_bams, pe)


# ── Step 2: compute bias_multi per dataset ────────────────────────────────────

print("\n=== Step 2: compute bias_multi ===")

all_multi = []

for gse in ["GSE123875","GSE123893","GSE123894","GSE135230","GSE175625"]:
    tissue   = TISSUES[gse]
    fvb_ids  = FVB_SRRS[gse]
    b6_ids   = B6_SRRS[gse]

    naive_multi = load_counts(f"{OUT}/{gse}.naive_multi.counts.txt")
    fvb_multi   = load_counts(f"{OUT}/{gse}.fvb_multi.counts.txt")

    fvb_naive_m = naive_multi[[c for c in fvb_ids if c in naive_multi.columns]]
    b6_naive_m  = naive_multi[[c for c in b6_ids  if c in naive_multi.columns]]
    fvb_fvb_m   = fvb_multi  [[c for c in fvb_ids if c in fvb_multi.columns]]

    # expression filter: CPM > 1 in ≥ half of samples in either group
    fvb_cpm_m = cpm(fvb_naive_m)
    b6_cpm_m  = cpm(b6_naive_m)
    n_fvb, n_b6 = fvb_naive_m.shape[1], b6_naive_m.shape[1]
    ok = ((fvb_cpm_m > MIN_CPM).sum(axis=1) >= max(1, n_fvb // 2)) | \
         ((b6_cpm_m  > MIN_CPM).sum(axis=1) >= max(1, n_b6  // 2))
    genes = ok[ok].index

    # borrow B6 normalisation from naive_multi for fvb arm
    b6_naive_cpm_m = b6_cpm_m.loc[genes]
    fvb_naive_cpm_m = fvb_cpm_m.loc[genes]
    fvb_fvb_cpm_m   = cpm(fvb_fvb_m).loc[genes.intersection(fvb_fvb_m.index)]
    common = fvb_fvb_cpm_m.index

    fc_naive_m = mean_log2fc(fvb_naive_cpm_m.loc[common], b6_naive_cpm_m.loc[common])
    fc_fvb_m   = mean_log2fc(fvb_fvb_cpm_m,               b6_naive_cpm_m.loc[common])
    bias_m     = fc_fvb_m - fc_naive_m

    df = pd.DataFrame({
        "gene_id":       common,
        "gse":           gse,
        "tissue":        tissue,
        "log2fc_naive_multi": fc_naive_m.values,
        "log2fc_fvb_multi":   fc_fvb_m.values,
        "bias_multi":         bias_m.values,
    })
    all_multi.append(df)
    print(f"  {gse}: {len(common)} genes, mean bias_multi={bias_m.mean():+.4f}")

multi_df = pd.concat(all_multi, ignore_index=True)
multi_df.to_csv(f"{OUT}/bias_multi_all.tsv", sep="\t", index=False)
print(f"\nSaved: {OUT}/bias_multi_all.tsv")


# ── Step 3: join bias_unique (existing) with bias_multi ──────────────────────

print("\n=== Step 3: join unique vs multi ===")

unique_df = pd.read_csv(f"{STAGE6}/all_datasets_bias.tsv", sep="\t")
unique_df = unique_df[["gene_id","gse","tissue","log2fc_naive","log2fc_fvb","bias"]].copy()
unique_df.rename(columns={"bias":"bias_unique"}, inplace=True)

merged = unique_df.merge(
    multi_df[["gene_id","gse","bias_multi"]],
    on=["gene_id","gse"], how="inner"
)
merged["delta_bias"] = merged["bias_multi"] - merged["bias_unique"]

merged.to_csv(f"{OUT}/bias_unique_vs_multi.tsv", sep="\t", index=False)
print(f"Saved: {OUT}/bias_unique_vs_multi.tsv  ({len(merged)} rows)")

# summary stats
print(f"\nOverall Δ stats across {len(merged)} gene-dataset pairs:")
print(f"  mean Δ  = {merged['delta_bias'].mean():+.4f}")
print(f"  median Δ = {merged['delta_bias'].median():+.4f}")
print(f"  std Δ   = {merged['delta_bias'].std():.4f}")
print(f"  |Δ| > 0.5 in {(merged['delta_bias'].abs() > 0.5).sum()} pairs ({(merged['delta_bias'].abs() > 0.5).mean()*100:.1f}%)")


# ── Step 4: catalog gene classification ──────────────────────────────────────

print("\n=== Step 4: catalog gene classification ===")

catalog = pd.read_csv(f"{STAGE6}/catalog_final.tsv", sep="\t")
catalog_ids = set(catalog["gene_id"])

cat_df = merged[merged["gene_id"].isin(catalog_ids)].copy()

# per-gene summary across tissues (excluding PyVT)
cat_4t = cat_df[cat_df["tissue"] != "mammary_PyVT"]
per_gene = cat_4t.groupby("gene_id").agg(
    bias_unique_mean   = ("bias_unique","mean"),
    bias_multi_mean    = ("bias_multi","mean"),
    delta_mean         = ("delta_bias","mean"),
    delta_abs_mean     = ("delta_bias", lambda x: x.abs().mean()),
    n_tissues          = ("tissue","nunique"),
).reset_index()

# merge with catalog gene name
per_gene = per_gene.merge(catalog[["gene_id","gene_name","snp_density_per_kb",
                                    "mean_bias","n_datasets_4t"]],
                           on="gene_id", how="left")

# classify: SNP-driven vs paralog-driven
DELTA_THRESH = 0.3
per_gene["mechanism"] = "SNP-driven"
per_gene.loc[per_gene["delta_abs_mean"] > DELTA_THRESH, "mechanism"] = "paralog-influenced"

per_gene = per_gene.sort_values("delta_abs_mean", ascending=False)
per_gene.to_csv(f"{OUT}/catalog_mechanism.tsv", sep="\t", index=False)
print(f"Saved: {OUT}/catalog_mechanism.tsv")

print(f"\nCatalog gene Δ distribution (4 non-PyVT tissues):")
print(per_gene[["gene_name","snp_density_per_kb","bias_unique_mean",
                "bias_multi_mean","delta_mean","delta_abs_mean","mechanism"]].to_string(index=False))

n_snp = (per_gene["mechanism"] == "SNP-driven").sum()
n_par = (per_gene["mechanism"] == "paralog-influenced").sum()
print(f"\nClassification: {n_snp} SNP-driven, {n_par} paralog-influenced (|Δ| threshold = {DELTA_THRESH})")

# highlight ribosomal genes
ribo = per_gene[per_gene["gene_name"].str.startswith(("Rpl","Rps","Fau","H4"))]
if len(ribo):
    print(f"\nRibosomal/H4 genes in catalog:")
    print(ribo[["gene_name","bias_unique_mean","delta_mean","mechanism"]].to_string(index=False))

print("\n=== Done ===")
