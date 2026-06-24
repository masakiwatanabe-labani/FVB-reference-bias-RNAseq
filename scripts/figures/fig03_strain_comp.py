#!/usr/bin/env python3
"""
Fig 3: 3-arm strain comparison — per-gene scatter as primary panel
Panel A: Scatter naive vs fvb log2FC (GSE123893 liver) + catalog gene overlay
Panel B: Bias (log2FC_fvb − log2FC_naive) distribution per dataset (violin)
Panel C: Dataset-level mean log2FC ± SE, 3 arms (secondary context; WASP note)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats
import os
import pathlib

PROJ = os.environ.get("FVB_PROJ", str(pathlib.Path(__file__).parent.parent.parent.resolve()))
OUT  = f"{PROJ}/results/figures"
os.makedirs(OUT, exist_ok=True)

df  = pd.read_csv(f"{PROJ}/results/stage6/all_datasets_bias.tsv", sep="\t")
df  = df.dropna(subset=["log2fc_naive","log2fc_wasp","log2fc_fvb"])
cat = pd.read_csv(f"{PROJ}/results/stage6/catalog_final.tsv", sep="\t")
# catalog_final.tsv = 32 clean genes (4-tissue criterion; PyVT-dependent excluded)

TISSUES = {
    "GSE123875": "Adipose", "GSE123893": "Liver", "GSE123894": "Kidney",
    "GSE135230": "Heart",   "GSE175625": "Mammary\n(PyVT*)",
}
GSES   = ["GSE123875","GSE123893","GSE123894","GSE135230","GSE175625"]
BLUE   = "#2166AC"
ORANGE = "#E08214"
GREEN  = "#1A9850"
GRAY   = "#888888"
RED    = "#CC0000"
LGRAY  = "#CCCCCC"

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 5.5))
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(1, 3, width_ratios=[2.0, 1.0, 1.0],
                      wspace=0.38, left=0.06, right=0.97,
                      top=0.87, bottom=0.14)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])
axC = fig.add_subplot(gs[2])

# ── Panel A: Scatter naive vs fvb — catalog genes highlighted ─────────────────
REP = "GSE123893"
rep = df[df["gse"] == REP].copy()
rep = rep[(rep["log2fc_naive"].between(-4,4)) & (rep["log2fc_fvb"].between(-4,4))]

# Background: all expressed genes
axA.scatter(rep["log2fc_naive"], rep["log2fc_fvb"],
            s=3, alpha=0.14, color=LGRAY, edgecolors="none",
            rasterized=True, zorder=1)

# Catalog overlay colored by n_datasets_4t (non-PyVT datasets with bias≥0.5)
n_col  = {4: RED, 3: "#D04030", 2: ORANGE}
n_size = {4: 42,  3: 32,        2: 25}
cat_rep = (rep[rep["gene_id"].isin(cat["gene_id"])]
           .merge(cat[["gene_id","gene_name","n_datasets_4t","mean_bias"]],
                  on="gene_id", how="left")
           .sort_values("n_datasets_4t"))

for nd in [2, 3, 4]:
    sub = cat_rep[cat_rep["n_datasets_4t"] == nd]
    if len(sub) == 0:
        continue
    axA.scatter(sub["log2fc_naive"], sub["log2fc_fvb"],
                s=n_size[nd], color=n_col[nd], edgecolors="white", lw=0.5,
                alpha=0.92, zorder=4,
                label=f"Catalog: {nd}/4 tissues (n={len(sub)})")

# Labels: top 7 by mean_bias that appear in GSE123893
top7 = cat.nlargest(7, "mean_bias")
for _, row in cat_rep[cat_rep["gene_id"].isin(top7["gene_id"])].iterrows():
    x, y = row["log2fc_naive"], row["log2fc_fvb"]
    axA.annotate(
        row["gene_name"], xy=(x, y),
        xytext=(x + 0.18, y + 0.18),
        fontsize=8, color="#990000", fontweight="bold",
        arrowprops=dict(arrowstyle="-", color="#990000", lw=0.5, shrinkB=3),
        zorder=6
    )

# Diagonal and reference lines
axA.plot([-4,4], [-4,4], "k--", lw=0.8, alpha=0.4, zorder=2)
axA.axhline(0, color=GREEN, lw=0.7, ls=":", alpha=0.45, zorder=2)
axA.axvline(0, color=BLUE,  lw=0.7, ls=":", alpha=0.45, zorder=2)

axA.set_xlim(-4, 4)
axA.set_ylim(-4, 4)
axA.set_xlabel("log₂FC naive (B6 reference)", fontsize=10.5)
axA.set_ylabel("log₂FC FVB reference (g2gtools)", fontsize=10.5)
axA.spines[["top","right"]].set_visible(False)
axA.set_title(
    "Catalog genes systematically above diagonal\n"
    f"(GSE123893, liver; {len(rep):,} expressed genes)",
    fontsize=11, pad=6
)
axA.text(-0.09, 1.04, "A", transform=axA.transAxes, fontsize=15, fontweight="bold")

# Legend
leg = [mpatches.Patch(facecolor=n_col[nd], alpha=0.9,
                       label=f"Catalog: {nd}/4 tissues")
       for nd in [4,3,2] if (cat_rep["n_datasets_4t"]==nd).any()]
leg.append(mpatches.Patch(facecolor=LGRAY, alpha=0.8, label="All expressed genes"))
axA.legend(handles=leg, fontsize=8, loc="upper left",
           frameon=True, framealpha=0.92, edgecolor="#CCCCCC")

axA.text(0.97, 0.04,
         f"Catalog: {len(cat)} genes\n(bias≥0.5 in ≥2 of 4 normal tissues)",
         transform=axA.transAxes, fontsize=8.5, va="bottom", ha="right",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                   edgecolor="#CCCCCC", alpha=0.9))

# ── Panel B: Bias distribution per dataset (violin) ───────────────────────────
tissue_colors = {
    "GSE123875": "#1F78B4", "GSE123893": "#33A02C", "GSE123894": "#FF7F00",
    "GSE135230": "#E31A1C", "GSE175625": "#6A3D9A",
}
bias_per_ds = []
for gse in GSES:
    vals = df[df["gse"] == gse]["bias"].dropna().clip(-3, 3).values
    bias_per_ds.append(vals)

positions = np.arange(len(GSES))
vp = axB.violinplot(bias_per_ds, positions=positions, widths=0.7,
                    showmedians=False, showextrema=False)
for pc, gse in zip(vp["bodies"], GSES):
    pc.set_facecolor(tissue_colors[gse])
    pc.set_alpha(0.75)
    pc.set_edgecolor("white")
    pc.set_linewidth(0.5)

# Median dots
for i, (gse, vals) in enumerate(zip(GSES, bias_per_ds)):
    axB.scatter(i, np.median(vals), s=22, color="white",
                edgecolors=tissue_colors[gse], linewidths=1.3, zorder=5)
    axB.scatter(i, np.percentile(vals, 95), s=8, color=tissue_colors[gse],
                alpha=0.6, zorder=4, marker="^")

axB.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
axB.axhline(0.5, color=RED, lw=0.7, ls=":", alpha=0.5)
axB.text(len(GSES)-0.45, 0.5, "catalog\nthreshold", fontsize=8, color=RED,
         va="bottom", ha="right")

axB.set_xticks(positions)
axB.set_xticklabels([TISSUES[g] for g in GSES], fontsize=8.5, rotation=28, ha="right")
axB.set_ylabel("Bias (log₂FC_fvb − log₂FC_naive)", fontsize=9)
axB.spines[["top","right"]].set_visible(False)
axB.set_title("Per-gene bias distribution\nper tissue", fontsize=10.5, pad=6)
axB.text(-0.24, 1.04, "B", transform=axB.transAxes, fontsize=15, fontweight="bold")
axB.text(len(GSES)-0.45, axB.get_ylim()[0]*0.92 if axB.get_ylim()[0] < 0 else -2.7,
         "* CNVs (PyVT)", fontsize=8, color=GRAY, ha="right", style="italic")

# ── Panel C: Mean log2FC ± SE per arm (secondary context) ────────────────────
summ = df.groupby("gse").agg(
    naive_m=("log2fc_naive","mean"),
    naive_se=("log2fc_naive", lambda x: x.std(ddof=1)/np.sqrt(len(x))),
    wasp_m= ("log2fc_wasp", "mean"),
    wasp_se=("log2fc_wasp",  lambda x: x.std(ddof=1)/np.sqrt(len(x))),
    fvb_m=  ("log2fc_fvb",  "mean"),
    fvb_se= ("log2fc_fvb",   lambda x: x.std(ddof=1)/np.sqrt(len(x))),
).loc[GSES]

x = np.arange(len(GSES))
w = 0.7
axC.errorbar(x - w/2, summ["naive_m"], yerr=summ["naive_se"],
             fmt="o-", color=BLUE,   ms=6, lw=1.4, capsize=3, label="Naive", zorder=4)
axC.errorbar(x,       summ["wasp_m"],  yerr=summ["wasp_se"],
             fmt="s--", color=ORANGE, ms=6, lw=1.2, capsize=3, label="WASP",   zorder=3,
             alpha=0.75)
axC.errorbar(x + w/2, summ["fvb_m"],  yerr=summ["fvb_se"],
             fmt="^-", color=GREEN,   ms=6, lw=1.4, capsize=3, label="FVB ref", zorder=4)

# WASP note
axC.text(0.04, 0.07,
         "WASP Δ from naive: ≤0.005 log₂ (all tissues)\n"
         "Attrition: 0.2–2.3% (vW=2 only; both strains)\n"
         "B6 excess attrition vs FVB (Table S4)",
         transform=axC.transAxes, fontsize=8, va="bottom", color=ORANGE,
         style="italic",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                   edgecolor=ORANGE, alpha=0.85, lw=0.8))

axC.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
axC.set_xticks(x)
axC.set_xticklabels([TISSUES[g] for g in GSES], fontsize=8.5, rotation=28, ha="right")
axC.set_ylabel("Mean log₂FC ± SE", fontsize=9)
axC.spines[["top","right"]].set_visible(False)
axC.legend(fontsize=8, loc="lower left", frameon=True,
           framealpha=0.9, edgecolor="#CCCCCC")
axC.set_title("Dataset-level summary\n(mean ± SE, all 3 arms)", fontsize=10.5, pad=6)
axC.text(-0.27, 1.04, "C", transform=axC.transAxes, fontsize=15, fontweight="bold")

fig.suptitle(
    "Personalised reference alignment stably reduces locus-specific SNP-mediated reference-genome mapping bias;\n"
    "32 catalog genes above diagonal in naive but not FVB-ref mapping"
    " (WASP: negligible genome-wide effect, Δ ≤ 0.005 log₂; attrition 0.2–2.3% [vW=2]; Table S4)",
    fontsize=9.5, y=0.997, color="#333333"
)

for ext in ["png","pdf"]:
    p = f"{OUT}/fig3_strain_comp.{ext}"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {p}")
plt.close(fig)
print("Fig 3 done.")
