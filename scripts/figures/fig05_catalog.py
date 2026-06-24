#!/usr/bin/env python3
"""
Fig 5: Bias catalog
Panel A: Heatmap — top N consistently biased genes × datasets (bias value)
Panel B: Dot plot — mean bias vs n_datasets, sized by mean_bias, colored by tissue breadth
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd
import re, os

PROJ = os.environ.get("FVB_PROJ", str(pathlib.Path(__file__).parent.parent.parent.resolve()))
OUT  = f"{PROJ}/results/figures"
GTF  = f"{PROJ}/ref/b6_gtf/GRCm39.110.gtf"
os.makedirs(OUT, exist_ok=True)

BLUE    = "#2166AC"
RED     = "#CC0000"
GRAY    = "#888888"
ORANGE  = "#E08214"
GREEN   = "#1A9850"

# ── 1. Load data ──────────────────────────────────────────────────────────────
print("Loading bias tables...")
df  = pd.read_csv(f"{PROJ}/results/stage6/all_datasets_bias.tsv", sep="\t")
# catalog_final.tsv = 32 genes (4-tissue criterion; PyVT-dependent genes excluded)
cat_full = pd.read_csv(f"{PROJ}/results/stage6/catalog_final.tsv", sep="\t")
catalog = cat_full.sort_values(["n_datasets_4t","mean_bias"], ascending=[False, False]).copy()
print(f"  Catalog: {len(catalog)} genes (bias≥0.5 in ≥2 of 4 normal tissues)")

# Top 30 for heatmap
top_n = min(30, len(catalog))
top_genes = catalog.head(top_n)

GSES    = ["GSE123875","GSE123893","GSE123894","GSE135230","GSE175625"]
TISSUES = {
    "GSE123875": "Adipose", "GSE123893": "Liver", "GSE123894": "Kidney",
    "GSE135230": "Heart",   "GSE175625": "Mammary\n(PyVT*)",
}

# Build heatmap matrix (raw per-dataset bias values)
hm_data = np.full((top_n, len(GSES)), np.nan)
for i, row in enumerate(top_genes.itertuples()):
    for j, gse in enumerate(GSES):
        sub = df[(df["gene_id"] == row.gene_id) & (df["gse"] == gse)]
        if len(sub) == 1:
            hm_data[i, j] = sub["bias"].values[0]

# Gene labels: "‡" for pvt_discordant (positive 4-tissue mean but negative PyVT bias)
pvt_disc_ids = set(catalog.loc[catalog["pvt_discordant"], "gene_id"]) \
    if "pvt_discordant" in catalog.columns else set()

y_labels = []
for row in top_genes.itertuples():
    suffix = "‡" if row.gene_id in pvt_disc_ids else ""
    y_labels.append(f"{row.gene_name}{suffix}")
x_labels = [TISSUES[g] for g in GSES]

# ── 3. Figure ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 9))
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1.0],
                      wspace=0.40, left=0.16, right=0.97,
                      top=0.88, bottom=0.09)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])

axA.text(-0.40, 1.03, "A", transform=axA.transAxes, fontsize=15, fontweight="bold")
axB.text(-0.22, 1.03, "B", transform=axB.transAxes, fontsize=15, fontweight="bold")

# ── Panel A: heatmap ──────────────────────────────────────────────────────────
vmax = np.nanpercentile(np.abs(hm_data), 95)
vmax = max(vmax, 1.0)
cmap = plt.get_cmap("RdBu_r")
im = axA.imshow(hm_data, aspect="auto", cmap=cmap,
                vmin=-vmax, vmax=vmax, interpolation="nearest")

axA.set_xticks(range(len(GSES)))
axA.set_xticklabels(x_labels, fontsize=9.5)
axA.set_yticks(range(top_n))
axA.set_yticklabels(y_labels, fontsize=8.0)
axA.set_title(f"Top {top_n} of {len(catalog)} consistently biased genes\n"
              f"(‡ = negative PyVT bias; bias≥0.5 in ≥2 of 4 normal tissues)",
              fontsize=11, pad=8)

# Add value text for cells with large bias
for i in range(top_n):
    for j in range(len(GSES)):
        v = hm_data[i, j]
        if np.isnan(v):
            axA.text(j, i, "—", ha="center", va="center",
                     fontsize=8, color="#AAAAAA")
        else:
            tc = "white" if abs(v) > vmax*0.6 else "black"
            axA.text(j, i, f"{v:.1f}", ha="center", va="center",
                     fontsize=8, color=tc, fontweight="bold" if abs(v) >= 1 else "normal")

cbar = fig.colorbar(im, ax=axA, fraction=0.025, pad=0.03)
cbar.set_label("Bias (log₂FC FVB ref − naive)", fontsize=9)
cbar.ax.tick_params(labelsize=8)
axA.text(len(GSES)-0.5, top_n+0.8, "* CNV confounding possible (PyVT tumour)",
         fontsize=8, color=GRAY, ha="right", style="italic")

# ── Panel B: dot plot (mean_bias vs tissue breadth) ───────────────────────────
bt_map = {
    "protein_coding": BLUE,
    "lncRNA": GREEN,
    "lincRNA": GREEN,
    "processed_pseudogene": ORANGE,
    "unprocessed_pseudogene": ORANGE,
    "pseudogene": ORANGE,
}
bt_default = GRAY

# mean_bias = 4-tissue mean (primary); pvt_discordant marks genes where PyVT is negative
disc_mask = catalog["pvt_discordant"].values if "pvt_discordant" in catalog.columns else \
            np.zeros(len(catalog), dtype=bool)

x_jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(catalog))
x_vals = catalog["n_datasets_4t"].values + x_jitter
y_vals = catalog["mean_bias"].values
sizes  = np.clip(y_vals * 18, 20, 200)
colors = [bt_map.get(bt, bt_default) for bt in catalog["biotype"].values]

axB.scatter(x_vals[~disc_mask], y_vals[~disc_mask],
            s=sizes[~disc_mask], c=[c for c, d in zip(colors, disc_mask) if not d],
            alpha=0.78, edgecolors="white", linewidths=0.5, zorder=3)
axB.scatter(x_vals[disc_mask], y_vals[disc_mask],
            s=sizes[disc_mask], c=[c for c, d in zip(colors, disc_mask) if d],
            alpha=0.78, edgecolors="#AA4400", linewidths=1.2, zorder=4,
            marker="D")

# Annotate top 10
top10 = catalog.head(10)
for row in top10.itertuples():
    yv = row.mean_bias
    xi = row.n_datasets_4t + np.random.default_rng(int(hash(row.gene_id)) % 2**31).uniform(-0.08, 0.08)
    lbl = f"{row.gene_name}‡" if row.gene_id in pvt_disc_ids else row.gene_name
    axB.annotate(lbl, xy=(xi, yv), xytext=(6, 0), textcoords="offset points",
                 fontsize=8, va="center",
                 arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.5))

axB.axhline(0.5, color=RED, lw=0.8, ls="--", alpha=0.6, label="Threshold (0.5)")

# Star-mark Rsph3a as most robustly evidenced catalog gene
rsph3a_cat = catalog[catalog["gene_name"] == "Rsph3a"]
if len(rsph3a_cat) == 1:
    rx = float(rsph3a_cat["n_datasets_4t"].values[0])
    ry = float(rsph3a_cat["mean_bias"].values[0])
    axB.scatter(rx, ry, s=220, marker="*", color="gold",
                edgecolors="#CC8800", lw=1.0, zorder=6,
                label="Most robust (Rsph3a)")
    axB.annotate("← most robust", xy=(rx, ry),
                 xytext=(rx + 0.12, ry + 0.04),
                 fontsize=8, color="#CC6600", va="center",
                 style="italic")
axB.set_xlabel("Normal tissues with bias ≥ 0.5 (of 4)", fontsize=10)
axB.set_ylabel("Mean bias (4 normal tissues)\n(log₂FC FVB ref − naive)", fontsize=10)
axB.set_xticks([1, 2, 3, 4])
axB.set_xlim(0.5, 4.5)
axB.spines[["top","right"]].set_visible(False)
axB.set_title(f"Bias catalog overview\n"
              f"({len(catalog)} genes; ‡ = negative PyVT bias)",
              fontsize=11, pad=8)

from matplotlib.lines import Line2D
legend_patches = [
    mpatches.Patch(facecolor=BLUE,   label="Protein-coding", alpha=0.78),
    mpatches.Patch(facecolor=GREEN,  label="lncRNA",         alpha=0.78),
    mpatches.Patch(facecolor=ORANGE, label="Pseudogene",     alpha=0.78),
    mpatches.Patch(facecolor=GRAY,   label="Other",          alpha=0.78),
    Line2D([0],[0], marker="D", color="w", markeredgecolor="#AA4400",
           markerfacecolor=GRAY, markersize=8, label="‡ PyVT discordant"),
    Line2D([0],[0], marker="*", color="w", markerfacecolor="gold",
           markeredgecolor="#CC8800", markersize=10, label="Most robustly evidenced"),
]
axB.legend(handles=legend_patches, fontsize=8, loc="upper left",
           frameon=True, framealpha=0.9, edgecolor="#CCCCCC",
           title="Biotype", title_fontsize=8.5)

fig.suptitle(
    "Bias catalog: consistently biased genes identify loci vulnerable to FVB reference genome mapping errors",
    fontsize=9.5, y=0.995, color="#333333"
)

for ext in ["png","pdf"]:
    p = f"{OUT}/fig5_catalog.{ext}"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {p}")
plt.close(fig)
print("Fig 5 done.")
