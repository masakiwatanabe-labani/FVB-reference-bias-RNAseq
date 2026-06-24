#!/usr/bin/env python3
"""
Fig 1: Survey of reference genome practice in FVB RNA-seq studies (GEO)
Panel A: 154 studies — B6 reference universal, strain-aware correction = 0
Panel B: Schematic of the reference bias problem
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import os
import pathlib

PROJ = os.environ.get("FVB_PROJ", str(pathlib.Path(__file__).parent.parent.parent.resolve()))
OUT  = f"{PROJ}/results/figures"
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(f"{PROJ}/results/survey/survey_of_practice.tsv", sep="\t")
n_total       = len(df)           # 154
n_b6_stated   = (df["category"] == "B6_stated").sum()       # 1
n_not_stated  = (df["category"] == "not_stated").sum()      # 153
n_strain_aware = 0

# ── Figure layout ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(11, 5))
fig.patch.set_facecolor("white")

gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1], wspace=0.35,
                      left=0.07, right=0.97, top=0.88, bottom=0.12)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])

# ── Panel A: stacked bar ───────────────────────────────────────────────────────
BLUE   = "#2166AC"
LBLUE  = "#92C5DE"
GRAY   = "#CCCCCC"
RED    = "#D6604D"

categories  = ["Reference genome\nused", "Strain-aware\ncorrection applied"]
b6_vals     = [n_total,          0]
none_vals   = [0,                n_total]

x = np.array([0, 1])
w = 0.45

bars_b6   = axA.bar(x, b6_vals,   width=w, color=BLUE,  label="GRCm39/38 (B6)",  zorder=3)
bars_none = axA.bar(x, none_vals, width=w, color=GRAY,  label="None / not applied", zorder=3)

# 数値アノテーション
for bar, val in zip(bars_b6, b6_vals):
    if val > 0:
        axA.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                 f"{val}\n({val/n_total*100:.0f}%)",
                 ha="center", va="bottom", fontsize=11, fontweight="bold", color=BLUE)

axA.text(x[1], 2, f"0\n(0%)",
         ha="center", va="bottom", fontsize=11, fontweight="bold", color="#888888")

axA.set_xticks(x)
axA.set_xticklabels(categories, fontsize=10.5)
axA.set_ylabel("Number of studies", fontsize=11)
axA.set_ylim(0, n_total * 1.22)
axA.set_xlim(-0.5, 1.5)
axA.spines[["top","right"]].set_visible(False)
axA.yaxis.set_tick_params(labelsize=10)

leg = axA.legend(loc="upper right", fontsize=9.5, frameon=True,
                 framealpha=0.9, edgecolor="#CCCCCC")

axA.set_title(f"Survey of {n_total} FVB RNA-seq studies (GEO)", fontsize=11.5, pad=8)
axA.text(-0.12, 1.04, "A", transform=axA.transAxes,
         fontsize=15, fontweight="bold")

# ── Panel B: bias schematic ────────────────────────────────────────────────────
axB.set_xlim(0, 10)
axB.set_ylim(0, 10)
axB.axis("off")
axB.text(-0.05, 1.04, "B", transform=axB.transAxes,
         fontsize=15, fontweight="bold")
axB.set_title("Reference bias in F1 hybrid RNA-seq", fontsize=11.5, pad=8)

def draw_box(ax, x, y, w, h, color, label, fontsize=9, text_color="white"):
    rect = mpatches.FancyBboxPatch((x, y), w, h,
        boxstyle="round,pad=0.15", facecolor=color, edgecolor="none", zorder=3)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, label, ha="center", va="center",
            fontsize=fontsize, color=text_color, fontweight="bold", zorder=4)

def arrow(ax, x1, y1, x2, y2, color="#555555", lw=1.5):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", color=color, lw=lw))

# F1 mouse
draw_box(axB, 0.3, 7.8, 3.6, 1.0, "#4DAF4A", "F1 mouse\n(C57BL/6 × FVB)", fontsize=8.5)

# RNA-seq reads (two alleles)
draw_box(axB, 0.3, 5.6, 1.6, 1.4, BLUE,  "B6 allele\nreads", fontsize=8)
draw_box(axB, 2.3, 5.6, 1.6, 1.4, "#E08214", "FVB allele\nreads", fontsize=8)

arrow(axB, 1.3, 7.8, 1.1, 7.0)
arrow(axB, 2.7, 7.8, 3.1, 7.0)

# Reference
draw_box(axB, 5.5, 5.6, 3.8, 1.4, BLUE, "GRCm39\n(B6 reference)", fontsize=8.5)

# Alignment arrows
arrow(axB, 3.9, 6.3, 5.4, 6.3, color=BLUE, lw=2)      # B6 → ref (easy)
axB.annotate("", xy=(5.4, 6.0), xytext=(3.9, 6.0),
             arrowprops=dict(arrowstyle="->", color="#E08214", lw=2,
                             linestyle="dashed"))

# Labels on arrows
axB.text(4.6, 6.55, "maps well", fontsize=8, color=BLUE, ha="center")
axB.text(4.6, 5.65, "mis-mapped\n/ filtered", fontsize=8, color="#E08214",
         ha="center", style="italic")

# Result boxes
draw_box(axB, 5.5, 3.5, 1.7, 1.4, BLUE,   "B6 reads\ncounted", fontsize=8)
draw_box(axB, 7.5, 3.5, 1.7, 1.4, GRAY,   "FVB reads\nlost", fontsize=8,
         text_color="#555555")

arrow(axB, 7.2, 5.6, 6.3, 4.9)
arrow(axB, 8.2, 5.6, 8.3, 4.9)

# Bias result
draw_box(axB, 5.2, 1.5, 4.2, 1.4, RED,
         "Apparent FVB downregulation\n(reference bias artifact)", fontsize=8)
arrow(axB, 7.2, 3.5, 7.2, 2.9, color=RED, lw=2)

# FVB allele fraction annotation
axB.text(0.5, 4.6, "Expected FVB allele fraction: 0.50\n"
         "Observed (naive B6 mapping): ~0.33",
         fontsize=8.5, color="#333333",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF3CD",
                   edgecolor="#DDAA00", linewidth=1))

fig.suptitle("", fontsize=1)

for ext in ["png", "pdf"]:
    path = f"{OUT}/fig1_survey.{ext}"
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {path}")

plt.close(fig)
print("Fig 1 done.")
