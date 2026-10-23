#!/usr/bin/env python3
"""
Supplementary Figure S_multi: Multi-mapping sensitivity analysis
Four panels:
  A: genome-wide Δ distribution (histogram)
  B: catalog bias_unique vs bias_multimap scatter
  C: residual fraction horizontal bar chart for all 32 catalog genes
  D: representative genes per-tissue bars (Fau, Gm55594, Gm8909, Mir6236)
"""
import os
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

PROJ = os.environ.get("FVB_PROJ", str(pathlib.Path(__file__).parent.parent.parent.resolve()))
OUT   = f"{PROJ}/results/paralog_analysis"
FIGS  = f"{PROJ}/results/figures"
STAGE6= f"{PROJ}/results/stage6"
os.makedirs(FIGS, exist_ok=True)

# ── colours ──────────────────────────────────────────────────────────────────
C_MULTI  = "#E07B3A"   # orange: SNP-induced multi-mapping
C_DIRECT = "#4C72B0"   # blue: direct SNP misalignment
C_GREY   = "#AAAAAA"

# ── data ─────────────────────────────────────────────────────────────────────
merged   = pd.read_csv(f"{OUT}/bias_unique_vs_multi.tsv",   sep="\t")
catalog  = pd.read_csv(f"{STAGE6}/catalog_final.tsv",        sep="\t")
mech     = pd.read_csv(f"{OUT}/catalog_mechanism.tsv",       sep="\t")

cat_ids = set(catalog["gene_id"])

# 4T data
data_4t  = merged[merged["tissue"] != "mammary_PyVT"]
cat_4t   = data_4t[data_4t["gene_id"].isin(cat_ids)]

# Per-gene catalog summary
per_gene = cat_4t.groupby("gene_id").agg(
    bias_unique_mean = ("bias_unique","mean"),
    bias_multi_mean  = ("bias_multi","mean"),
    delta_mean       = ("delta_bias","mean"),
).reset_index()
per_gene = per_gene.merge(catalog[["gene_id","gene_name","mean_bias"]], on="gene_id")
per_gene["residual_fraction"] = per_gene["bias_multi_mean"] / per_gene["bias_unique_mean"]
THRESH = 0.3
per_gene["mclass"] = np.where(per_gene["delta_mean"].abs() > THRESH,
                               "SNP-induced multi-mapping", "Direct SNP misalignment")
per_gene = per_gene.sort_values("residual_fraction")

# ── figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 14))
gs  = GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.35)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

# ── Panel A: genome-wide Δ distribution ──────────────────────────────────────
delta_gw = data_4t["delta_bias"].clip(-1, 1)
bins = np.arange(-1.025, 1.05, 0.05)
ax_a.hist(delta_gw, bins=bins, color=C_GREY, alpha=0.7, edgecolor="none", zorder=1)

# catalog gene-dataset Δ (rug at top of histogram)
cat_deltas = cat_4t["delta_bias"].clip(-1, 1)
rug_y = ax_a.get_ylim()[1] * 0.97 if ax_a.get_ylim()[1] > 0 else 200
ax_a.scatter(cat_deltas, np.full(len(cat_deltas), rug_y),
             c=C_MULTI, s=12, alpha=0.6, zorder=3, label="Catalog genes")
ax_a.axvline(0, color="black", lw=1.0, ls="--")
ax_a.set_xlabel("Δ (bias_multi − bias_unique) [log₂]", fontsize=11)
ax_a.set_ylabel("Number of gene-dataset pairs", fontsize=11)
ax_a.set_title("A   Genome-wide Δ distribution", fontsize=12, fontweight="bold", loc="left")
ax_a.set_xlim(-1.05, 1.05)
# Annotate stats
ax_a.text(0.97, 0.95, f"n = {len(data_4t):,} pairs\nmedian Δ = −0.001 log₂\nmean Δ = −0.007 log₂",
          transform=ax_a.transAxes, ha="right", va="top", fontsize=9,
          bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

# ── Panel B: bias_unique vs bias_multimap scatter ─────────────────────────────
max_val = per_gene["bias_unique_mean"].max() * 1.1
ax_b.plot([0, max_val], [0, max_val], "--", color="grey", lw=1.0, zorder=1, label="Δ = 0")
for _, r in per_gene.iterrows():
    col = C_MULTI if r["mclass"] == "SNP-induced multi-mapping" else C_DIRECT
    ax_b.scatter(r["bias_unique_mean"], r["bias_multi_mean"], c=col, s=60, zorder=3, alpha=0.85)

# Label top-10 |Δ|
per_gene_lab = per_gene.copy()
per_gene_lab["abs_delta"] = per_gene_lab["delta_mean"].abs()
for _, r in per_gene_lab.nlargest(10, "abs_delta").iterrows():
    ax_b.annotate(r["gene_name"], (r["bias_unique_mean"], r["bias_multi_mean"]),
                  fontsize=7, ha="left", xytext=(4, 2), textcoords="offset points")

ax_b.set_xlabel("bias_unique (primary pipeline) [log₂]", fontsize=11)
ax_b.set_ylabel("bias_multimap (−M −−fraction) [log₂]", fontsize=11)
ax_b.set_title("B   Catalog genes: unique vs multi bias", fontsize=12, fontweight="bold", loc="left")
ax_b.axhline(0, color="grey", lw=0.5, ls=":")
handles = [mpatches.Patch(color=C_MULTI, label="SNP-induced multi-mapping (n=16)"),
           mpatches.Patch(color=C_DIRECT, label="Direct SNP misalignment (n=16)")]
ax_b.legend(handles=handles, fontsize=9, loc="upper left")

# ── Panel C: residual fraction bar chart ─────────────────────────────────────
y_pos   = np.arange(len(per_gene))
colors  = [C_MULTI if m == "SNP-induced multi-mapping" else C_DIRECT
           for m in per_gene["mclass"]]
rf_clip = per_gene["residual_fraction"].clip(-0.5, 1.7)
ax_c.barh(y_pos, rf_clip, color=colors, edgecolor="none", height=0.7, alpha=0.85)
ax_c.axvline(0, color="black", lw=0.8, ls="-")
ax_c.axvline(1, color="black", lw=1.0, ls="--")
ax_c.set_yticks(y_pos)
ax_c.set_yticklabels([f"{'*' if r['mclass']=='SNP-induced multi-mapping' else ''}{r['gene_name']}"
                       for _, r in per_gene.iterrows()], fontsize=8)
ax_c.set_xlabel("Residual fraction (bias_multi / bias_unique)", fontsize=11)
ax_c.set_title("C   Residual fraction per catalog gene", fontsize=12, fontweight="bold", loc="left")
ax_c.set_xlim(-0.55, 1.75)
ax_c.text(1.02, 0.01, "← multi-mapping\n   dominant", transform=ax_c.transAxes,
          fontsize=8, va="bottom", color=C_MULTI)

# ── Panel D: representative genes per-tissue ─────────────────────────────────
rep_genes = {
    "Fau":     ("ENSMUSG00000038274", "SNP-induced\nmulti-mapping"),
    "Gm55594": ("ENSMUSG00000095562", "SNP-induced\nmulti-mapping"),
    "Gm8909":  ("ENSMUSG00000064194", "Direct SNP\nmisalignment"),
    "Mir6236": ("ENSMUSG00000083773", "Direct SNP\nmisalignment"),
}

# Find correct IDs from catalog
rep_ids = {}
for gname, (guess_id, mclass) in rep_genes.items():
    hit = catalog[catalog["gene_name"] == gname]
    if len(hit):
        rep_ids[gname] = (hit.iloc[0]["gene_id"], mclass)
    else:
        rep_ids[gname] = (guess_id, mclass)

TISSUES = ["adipose","liver","kidney","heart"]
T_LABELS = {"adipose":"Adi", "liver":"Liv", "kidney":"Kid", "heart":"Hrt"}
bar_w = 0.35
gene_gap = 1.6
group_x = {g: i * gene_gap for i, g in enumerate(rep_genes)}

for gname, (gid, mclass) in rep_ids.items():
    sub = merged[(merged["gene_id"] == gid) & (merged["tissue"].isin(TISSUES))]
    sub = sub.set_index("tissue")
    col = C_MULTI if "multi-mapping" in mclass else C_DIRECT
    gx  = group_x[gname]
    for ti, tis in enumerate(TISSUES):
        if tis not in sub.index:
            continue
        xc = gx + (ti - 1.5) * bar_w * 0.55
        bu = sub.loc[tis, "bias_unique"]
        bm = sub.loc[tis, "bias_multi"]
        ax_d.bar(xc,        bu, width=bar_w*0.5, color=col,   alpha=0.85, label=None)
        ax_d.bar(xc + bar_w*0.5, bm, width=bar_w*0.5, color=col,   alpha=0.4,
                 hatch="///", edgecolor=col, label=None)

# Gene name labels
for gname, gx in group_x.items():
    mclass = rep_ids[gname][1]
    col = C_MULTI if "multi-mapping" in mclass else C_DIRECT
    ax_d.text(gx, -0.8, gname, ha="center", va="top", fontsize=9, color=col, fontweight="bold")
    ax_d.text(gx, -1.2, rep_ids[gname][1], ha="center", va="top", fontsize=7, color=col)

ax_d.axhline(0, color="black", lw=0.8)
ax_d.set_xticks([])
ax_d.set_ylabel("Bias [log₂FC_fvb − log₂FC_naive]", fontsize=11)
ax_d.set_title("D   Representative genes: unique vs multi bias", fontsize=12, fontweight="bold", loc="left")
handles_d = [mpatches.Patch(facecolor="grey", alpha=0.85, label="bias_unique (primary)"),
             mpatches.Patch(facecolor="grey", alpha=0.4, hatch="///", edgecolor="grey",
                            label="bias_multimap (−M −−fraction)")]
ax_d.legend(handles=handles_d, fontsize=9, loc="upper right")

# ── save ──────────────────────────────────────────────────────────────────────
for ext in ["png", "pdf"]:
    path = f"{FIGS}/sfig_multimap.{ext}"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    print(f"Saved: {path}")
plt.close()
print("Done.")
