#!/usr/bin/env python3
"""
Fig 2: F1 ASE bias — FVB allele fraction at heterozygous sites
Panel A: Per-sample dot + mean±SE (naive vs WASP, both datasets)
Panel B: Site-level scatter — naive vs WASP alt_frac (GSE293440 representative)
Panel C: Per-site histogram — naive vs WASP distribution (GSE293440, all 8 samples pooled)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats
import os
import pathlib, glob

PROJ = os.environ.get("FVB_PROJ", str(pathlib.Path(__file__).parent.parent.parent.resolve()))
COMB    = f"{PROJ}/results/f1_ase_combined"
ASE293  = f"{PROJ}/results/gse293440_f1_test"
ASE200  = f"{PROJ}/results/f1_ase_v2"
OUT     = f"{PROJ}/results/figures"
os.makedirs(OUT, exist_ok=True)

MIN_COV   = 20
MIN_MINOR = 3

BLUE     = "#2166AC"
LBLUE    = "#92C5DE"
ORANGE   = "#E08214"
LORANGE  = "#FDB863"
GRAY50   = "#888888"
RED_LINE = "#CC0000"

# ── helpers ────────────────────────────────────────────────────────────────────
def load_filt(path):
    if not os.path.exists(path):
        return None
    raw = pd.read_csv(path, sep="\t", comment="#")
    raw.columns = raw.columns.str.lower().str.replace(" ", "_")
    for old, new in [("refcount","ref_count"), ("altcount","alt_count"),
                     ("totalcount","total_count")]:
        if old in raw.columns:
            raw.rename(columns={old: new}, inplace=True)
    if "total_count" not in raw.columns:
        raw["total_count"] = raw["ref_count"] + raw["alt_count"]
    raw["alt_frac"] = raw["alt_count"] / raw["total_count"].replace(0, np.nan)
    filt = raw[(raw["total_count"] >= MIN_COV) &
               (raw["alt_count"]   >= MIN_MINOR) &
               (raw["ref_count"]   >= MIN_MINOR)].copy()
    return filt

# ── Load per-sample summary ────────────────────────────────────────────────────
df = pd.read_csv(f"{COMB}/f1_ase_all_samples.tsv", sep="\t")
df = df[df["mean_filt"].notna()].copy()

# ── Load Panel B data (representative sample: SRR32903247, 872 sites) ──────────
srr_b = "SRR32903247"
naive_b = load_filt(f"{ASE293}/{srr_b}.naive.ase.tsv")
wasp_b  = load_filt(f"{ASE293}/{srr_b}.wasp.ase.tsv")

# ── Load Panel C data (pool all GSE293440 samples) ─────────────────────────────
srrs_293 = [os.path.basename(p).replace(".naive.ase.tsv", "")
            for p in glob.glob(f"{ASE293}/*.naive.ase.tsv")]
naive_pool, wasp_pool = [], []
for srr in sorted(srrs_293):
    n = load_filt(f"{ASE293}/{srr}.naive.ase.tsv")
    w = load_filt(f"{ASE293}/{srr}.wasp.ase.tsv")
    if n is not None and w is not None and len(n) > 0 and len(w) > 0:
        merged = n.merge(w, on=["contig", "position"], suffixes=("_naive", "_wasp"))
        naive_pool.extend(merged["alt_frac_naive"].tolist())
        wasp_pool.extend(merged["alt_frac_wasp"].tolist())
naive_pool = np.array(naive_pool)
wasp_pool  = np.array(wasp_pool)

# ── significance: one-sample t-test vs 0.5 for WASP means ─────────────────────
def sig_star(vals):
    t, p = stats.ttest_1samp(vals, 0.5)
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "n.s."

wasp293_vals = df[(df["dataset"]=="GSE293440") & (df["mode"]=="wasp")]["mean_filt"].values
wasp200_vals = df[(df["dataset"]=="GSE200632") & (df["mode"]=="wasp")]["mean_filt"].values
star293 = sig_star(wasp293_vals)
star200 = sig_star(wasp200_vals)

# ── Figure layout ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 5.5))
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1.1, 1.2],
                      wspace=0.38, left=0.07, right=0.97,
                      top=0.86, bottom=0.13)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])
axC = fig.add_subplot(gs[2])

# ── Panel A: per-sample dot plot ───────────────────────────────────────────────
positions = {
    ("GSE293440","naive"): 0.0,
    ("GSE293440","wasp"):  0.7,
    ("GSE200632","naive"): 1.7,
    ("GSE200632","wasp"):  2.4,
}
colors = {
    ("GSE293440","naive"): BLUE,
    ("GSE293440","wasp"):  LBLUE,
    ("GSE200632","naive"): ORANGE,
    ("GSE200632","wasp"):  LORANGE,
}

for (gse, mode), xpos in positions.items():
    sub  = df[(df["dataset"]==gse) & (df["mode"]==mode)]
    vals = sub["mean_filt"].values
    np.random.seed(42)
    xs   = xpos + np.random.uniform(-0.08, 0.08, len(vals))
    axA.scatter(xs, vals, s=55, color=colors[(gse,mode)], zorder=5,
                edgecolors="white", linewidths=0.6, alpha=0.95)
    m, se = vals.mean(), vals.std(ddof=1)/np.sqrt(len(vals))
    axA.errorbar(xpos, m, yerr=se, fmt="_", color="black",
                 capsize=5, capthick=1.5, elinewidth=1.5, ms=12, zorder=6)

    # residual bias star for WASP columns
    if mode == "wasp":
        v   = wasp293_vals if gse == "GSE293440" else wasp200_vals
        st  = sig_star(v)
        axA.text(xpos, vals.mean() + vals.std(ddof=1)/np.sqrt(len(vals)) + 0.012,
                 st, ha="center", va="bottom", fontsize=10, color="black", zorder=7)

axA.axhline(0.5, color=RED_LINE, lw=1.2, ls="--", zorder=2, alpha=0.8)
axA.text(2.72, 0.502, "Expected\n(0.5)", va="bottom", ha="left",
         fontsize=8, color=RED_LINE)

for xc, label, color in [(0.35, "GSE293440\nfrontal cortex\n(B6J×FVB/NJ, N=8)", BLUE),
                          (2.05, "GSE200632\nkeratinocyte\n(B6×FVB/N, N=4)", ORANGE)]:
    axA.text(xc, 0.195, label, ha="center", va="top", fontsize=8.5,
             color=color, fontweight="bold")

axA.set_xlim(-0.45, 2.9)
axA.set_ylim(0.22, 0.58)
axA.set_xticks([0.0, 0.7, 1.7, 2.4])
axA.set_xticklabels(["Naive", "WASP", "Naive", "WASP"], fontsize=10)
axA.set_ylabel("FVB allele fraction\n(sample mean, filtered sites)", fontsize=10.5)
axA.spines[["top","right"]].set_visible(False)
axA.set_title("Per-sample FVB allele fraction", fontsize=11, pad=6)
axA.text(-0.13, 1.04, "A", transform=axA.transAxes,
         fontsize=15, fontweight="bold")

# 凡例: 下左に配置 ("Expected" 注記と重ならない)
leg_elements = [
    mpatches.Patch(facecolor=BLUE,    label="GSE293440 naive"),
    mpatches.Patch(facecolor=LBLUE,   label="GSE293440 WASP"),
    mpatches.Patch(facecolor=ORANGE,  label="GSE200632 naive"),
    mpatches.Patch(facecolor=LORANGE, label="GSE200632 WASP"),
]
axA.legend(handles=leg_elements, fontsize=8, loc="lower left",
           frameon=True, framealpha=0.9, edgecolor="#CCCCCC")
axA.text(1.45, 0.572, "* p<0.05 vs 0.5 (one-sample t-test)",
         fontsize=8, color="#555555", ha="right", va="top")

# ── Panel B: site-level scatter (naive vs WASP, representative sample) ─────────
if naive_b is not None and wasp_b is not None:
    merged = naive_b.merge(wasp_b, on=["contig","position"],
                           suffixes=("_naive","_wasp"))
    if len(merged) > 0:
        axB.scatter(merged["alt_frac_naive"], merged["alt_frac_wasp"],
                    s=12, alpha=0.35, color=BLUE, edgecolors="none", zorder=3)
        axB.plot([0,1],[0,1], "k--", lw=0.8, alpha=0.5, zorder=2)
        axB.axhline(0.5, color=RED_LINE, lw=0.9, ls=":", alpha=0.7, zorder=2)
        axB.axvline(0.5, color=RED_LINE, lw=0.9, ls=":", alpha=0.7, zorder=2)
        n  = len(merged)
        mn = merged["alt_frac_naive"].mean()
        mw = merged["alt_frac_wasp"].mean()
        axB.text(0.05, 0.96,
                 f"n = {n} sites\nNaive: {mn:.3f}\nWASP:  {mw:.3f}",
                 transform=axB.transAxes, fontsize=8.5, va="top",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                           edgecolor="#CCCCCC", alpha=0.9))

axB.set_xlim(0,1); axB.set_ylim(0,1)
axB.set_xlabel("FVB allele fraction (naive)",  fontsize=9.5)
axB.set_ylabel("FVB allele fraction (WASP)",   fontsize=9.5)
axB.spines[["top","right"]].set_visible(False)
axB.set_title(f"Site-level correction\n(GSE293440, {srr_b[:11]}…, representative)",
              fontsize=9.5, pad=6)
axB.text(-0.20, 1.04, "B", transform=axB.transAxes,
         fontsize=15, fontweight="bold")

# ── Panel C: per-site histogram (GSE293440, all 8 samples pooled) ──────────────
if len(naive_pool) > 0:
    bins = np.linspace(0, 1, 41)
    axC.hist(naive_pool, bins=bins, color=BLUE,  alpha=0.6,
             label=f"Naive  (mean={naive_pool.mean():.3f})", density=True, zorder=3)
    axC.hist(wasp_pool,  bins=bins, color=LBLUE, alpha=0.7,
             label=f"WASP   (mean={wasp_pool.mean():.3f})", density=True, zorder=4)

    axC.axvline(0.5, color=RED_LINE, lw=1.2, ls="--", zorder=5, alpha=0.9)
    axC.text(0.51, axC.get_ylim()[1]*0.01, "0.5", color=RED_LINE,
             fontsize=8, va="bottom")

    # mean lines
    axC.axvline(naive_pool.mean(), color=BLUE,  lw=1.2, ls=":", alpha=0.85, zorder=5)
    axC.axvline(wasp_pool.mean(),  color=LBLUE, lw=1.2, ls=":", alpha=0.85, zorder=5)

axC.set_xlabel("FVB allele fraction (per site)", fontsize=9.5)
axC.set_ylabel("Density", fontsize=9.5)
axC.spines[["top","right"]].set_visible(False)
axC.legend(fontsize=8, loc="upper left", frameon=True,
           framealpha=0.9, edgecolor="#CCCCCC")
n_sites = len(naive_pool)
axC.set_title(f"Per-site distribution\n(GSE293440, {len(srrs_293)} samples, {n_sites:,} site×sample)",
              fontsize=9.5, pad=6)
axC.text(-0.20, 1.04, "C", transform=axC.transAxes,
         fontsize=15, fontweight="bold")

# ── y-axis fix after hist fill ─────────────────────────────────────────────────
if len(naive_pool) > 0:
    axC.axvline(0.5, color=RED_LINE, lw=1.2, ls="--", zorder=5, alpha=0.9)

fig.suptitle(
    "Reference bias in F1 hybrid RNA-seq: FVB allele under-represented\n"
    "relative to expected 0.5 (naive B6 mapping), partially corrected by WASP\n"
    "(WASP partial correction reflects sparse FVB/NJ SNP coverage [~1–2% reads overlap variant sites];\n"
    " pure-strain attrition 0.2–2.3% [vW=2 only]; Table S4)",
    fontsize=9.5, y=0.995, color="#333333"
)

for ext in ["png","pdf"]:
    path = f"{OUT}/fig2_f1_ase.{ext}"
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {path}")

plt.close(fig)
print("Fig 2 done.")
