#!/usr/bin/env python3
# Pipeline step 4: Bias quantification and catalog construction
#
# Input:  $FVB_PROJ/results/strain_comp/*.counts.txt  (from step 2)
# Output: $FVB_PROJ/results/stage6/
#
# Usage:
#   export FVB_PROJ=/path/to/working/directory
#   python 04_bias_analysis.py
#
"""
Stage 6: 3-arm strain comparison analysis  (v4 — B6-GTF-unified FVB arm)
Input:  results/strain_comp/*.counts.txt
Output: results/stage6/  (per-dataset TSVs + combined bias table)

Normalization: standard per-group CPM, B6 shared across arms.
  - naive arm: FVB_naive_CPM vs B6_naive_CPM
  - wasp  arm: FVB_wasp_CPM  vs B6_wasp_CPM
  - fvb   arm: FVB_b6gtf_CPM vs B6_naive_CPM

FVB arm uses fvb_b6gtf.counts.txt:
  FVB-aligned BAMs counted with B6 GTF (not FVB GTF).
  The VCI is SNP-only (no indels) → FVB genome coordinates are identical to
  B6 coordinates → B6 GTF applies directly to FVB-aligned BAMs.
  This eliminates the FVB-GTF / B6-GTF gene-set asymmetry (1,727 missing genes
  in FVB GTF) and makes all three arms use the same annotation.
  PE datasets: counted with -p flag (same as naive/wasp).
  SE datasets: counted without -p flag (same as naive/wasp).
  B6 samples in the fvb arm are borrowed from naive (identical B6 alignment).
  bias = log2(FVB_fvb_CPM+0.5) - log2(FVB_naive_CPM+0.5) [B6 cancels].
"""
import os, re
import pathlib
import numpy as np
import pandas as pd

PROJ = os.environ.get("FVB_PROJ", str(pathlib.Path(__file__).parent.parent.parent.resolve()))
COUNTS  = f"{PROJ}/results/strain_comp"
OUT     = f"{PROJ}/results/stage6"
os.makedirs(OUT, exist_ok=True)

MIN_CPM   = 1.0
PSEUDO    = 0.5

# ── Sample metadata ───────────────────────────────────────────────────────────
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
TISSUES = {
    "GSE123875": "adipose",
    "GSE123893": "liver",
    "GSE123894": "kidney",
    "GSE135230": "heart",
    "GSE175625": "mammary_PyVT",
}

# ── Helpers ───────────────────────────────────────────────────────────────────
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
    """Per-group CPM: each column divided by its own column sum."""
    lib = mat.sum(axis=0)
    return mat.divide(lib, axis=1) * 1e6

def expr_filter(fvb_cpm, b6_cpm):
    n_fvb = fvb_cpm.shape[1]
    n_b6  = b6_cpm.shape[1]
    ok_fvb = (fvb_cpm > MIN_CPM).sum(axis=1) >= max(1, n_fvb // 2)
    ok_b6  = (b6_cpm  > MIN_CPM).sum(axis=1) >= max(1, n_b6  // 2)
    return ok_fvb | ok_b6

def mean_log2fc(fvb_cpm, b6_cpm):
    mf = fvb_cpm.mean(axis=1) + PSEUDO
    mb = b6_cpm.mean(axis=1)  + PSEUDO
    return np.log2(mf / mb)

# ── Main loop ─────────────────────────────────────────────────────────────────
all_bias = []

for gse in ["GSE123875","GSE123893","GSE123894","GSE135230","GSE175625"]:
    tissue  = TISSUES[gse]
    fvb_ids = FVB_SRRS[gse]
    b6_ids  = B6_SRRS[gse]
    print(f"\n{'='*60}")
    print(f"{gse} ({tissue})")

    # ── load raw counts ────────────────────────────────────────────────────────
    naive_raw = load_counts(f"{COUNTS}/{gse}.naive.counts.txt")
    wasp_raw  = load_counts(f"{COUNTS}/{gse}.wasp.counts.txt")
    fvb_raw   = load_counts(f"{COUNTS}/{gse}.fvb_b6gtf.counts.txt")

    fvb_naive = naive_raw[[c for c in fvb_ids if c in naive_raw.columns]]
    b6_naive  = naive_raw[[c for c in b6_ids  if c in naive_raw.columns]]
    fvb_wasp  = wasp_raw [[c for c in fvb_ids if c in wasp_raw.columns]]
    b6_wasp   = wasp_raw [[c for c in b6_ids  if c in wasp_raw.columns]]
    fvb_fvb   = fvb_raw  [[c for c in fvb_ids if c in fvb_raw.columns]]

    # ── library size diagnostics ──────────────────────────────────────────────
    lib_naive_fvb = fvb_naive.sum(axis=0).mean()
    lib_fvb_b6    = fvb_fvb.sum(axis=0).mean()
    print(f"  FVB lib: naive={lib_naive_fvb:,.0f}  fvb_b6gtf={lib_fvb_b6:,.0f}  "
          f"ratio={lib_fvb_b6/lib_naive_fvb:.4f}")

    # ── per-group CPM ─────────────────────────────────────────────────────────
    fvb_naive_cpm = cpm(fvb_naive)
    b6_naive_cpm  = cpm(b6_naive)
    fvb_wasp_cpm  = cpm(fvb_wasp)
    b6_wasp_cpm   = cpm(b6_wasp)
    # fvb arm FVB: divide by NAIVE FVB library size (not fvb_fvbgtf total).
    # FVB GTF covers fewer genes than B6 GTF → fvb_fvbgtf total < naive total
    # → using fvb_fvbgtf total as denominator inflates ALL gene CPMs uniformly.
    # Using naive library size as denominator makes unbiased genes land on the
    # diagonal; only genuinely biased genes (more reads in fvb_fvbgtf) shift up.
    naive_fvb_lib = fvb_naive.sum(axis=0)          # per-sample, from naive arm
    fvb_fvb_cpm   = fvb_fvb.divide(naive_fvb_lib, axis=1) * 1e6
    # B6 in fvb arm: same alignment as naive → same CPM

    # ── expression filter (based on naive arm) ────────────────────────────────
    keep = expr_filter(fvb_naive_cpm, b6_naive_cpm)
    print(f"  Genes after filter: {keep.sum()} / {len(keep)}")

    # ── log2FC per arm ────────────────────────────────────────────────────────
    fc_naive = mean_log2fc(fvb_naive_cpm[keep], b6_naive_cpm[keep])
    fc_wasp  = mean_log2fc(fvb_wasp_cpm[keep],  b6_wasp_cpm[keep])

    common = keep[keep].index.intersection(fvb_fvb_cpm.index)
    # fvb arm: FVB from FVB-ref alignment, B6 from naive (same B6 GTF)
    fc_fvb = mean_log2fc(fvb_fvb_cpm.loc[common], b6_naive_cpm.loc[common])

    # ── compile per-gene result ───────────────────────────────────────────────
    genes = keep[keep].index.intersection(common)

    result = pd.DataFrame({
        "gene_id":      genes,
        "gse":          gse,
        "tissue":       tissue,
        "log2fc_naive": fc_naive.loc[genes].values,
        "log2fc_wasp":  fc_wasp.loc[genes].values,
        "log2fc_fvb":   fc_fvb.loc[genes].values,
    })
    # bias > 0: FVB undercounted in naive → reference bias
    result["bias"] = result["log2fc_fvb"] - result["log2fc_naive"]

    n_above = (result["log2fc_fvb"] > result["log2fc_naive"]).sum()
    print(f"  Naive  mean log2FC: {result['log2fc_naive'].mean():+.4f}")
    print(f"  WASP   mean log2FC: {result['log2fc_wasp'].mean():+.4f}")
    print(f"  FVB    mean log2FC: {result['log2fc_fvb'].mean():+.4f}")
    print(f"  Bias   mean       : {result['bias'].mean():+.4f}")
    print(f"  Above diagonal    : {100*n_above/len(result):.1f}%")

    out_path = f"{OUT}/{gse}_bias.tsv"
    result.to_csv(out_path, sep="\t", index=False)
    print(f"  → {out_path}")
    all_bias.append(result)

# ── Combined table ────────────────────────────────────────────────────────────
combined = pd.concat(all_bias, ignore_index=True)
combined.to_csv(f"{OUT}/all_datasets_bias.tsv", sep="\t", index=False)
print(f"\nCombined → {OUT}/all_datasets_bias.tsv  ({len(combined)} rows)")

# ── Bias catalog ──────────────────────────────────────────────────────────────
BIAS_THRESH = 0.5
biased      = combined[combined["bias"] >= BIAS_THRESH].groupby("gene_id")["gse"].nunique()
catalog_ids = biased[biased >= 2].index
catalog     = combined[combined["gene_id"].isin(catalog_ids)].copy()
catalog_summ = catalog.groupby("gene_id").agg(
    n_datasets=("gse",    "nunique"),
    mean_bias= ("bias",   "mean"),
    tissues=   ("tissue", lambda x: ",".join(sorted(set(x))))
).sort_values("mean_bias", ascending=False)
catalog_summ.to_csv(f"{OUT}/bias_catalog.tsv", sep="\t")
print(f"Bias catalog (≥2 datasets, bias≥{BIAS_THRESH}): {len(catalog_summ)} genes")
print(catalog_summ.head(10).to_string())
