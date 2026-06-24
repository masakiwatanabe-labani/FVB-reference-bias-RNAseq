#!/usr/bin/env python3
"""
Fig 6: Worked example — Fau (ENSMUSG00000038274, chr19)
Protein-coding gene biased in all 5 datasets (mean bias = +1.25 log2FC)

Panel A: Gene model (Fau-201) + FVB SNP positions
Panel B: Per-sample normalized counts (CPM) per arm — GSE135230 heart (highest bias tissue)
Panel C: log2FC per arm across all 5 datasets (the "bias" pattern)
Panel D: This gene on the SNP-density × bias scatter (context vs genome-wide)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import pysam, re, os

PROJ = os.environ.get("FVB_PROJ", str(pathlib.Path(__file__).parent.parent.parent.resolve()))
OUT    = f"{PROJ}/results/figures"
GTF    = f"{PROJ}/ref/b6_gtf/GRCm39.110.gtf"
VCF    = f"{PROJ}/ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf.gz"
COUNTS = f"{PROJ}/results/strain_comp"
os.makedirs(OUT, exist_ok=True)

BLUE   = "#2166AC"
ORANGE = "#E08214"
GREEN  = "#1A9850"
GRAY   = "#888888"
RED    = "#CC0000"
LGRAY  = "#DDDDDD"

GENE_ID   = "ENSMUSG00000038274"
GENE_NAME = "Fau"
CHROM     = "19"
GSTART    = 6107874
GEND      = 6109554
STRAND    = "+"
# Fau-201 canonical transcript exons
EXONS     = [(6107973,6108031),(6108276,6108355),(6108458,6108602),
             (6109174,6109229),(6109379,6109546)]
INTRONS   = [(EXONS[i][1], EXONS[i+1][0]) for i in range(len(EXONS)-1)]

GSES = ["GSE123875","GSE123893","GSE123894","GSE135230","GSE175625"]
TISSUES = {
    "GSE123875": "Adipose", "GSE123893": "Liver", "GSE123894": "Kidney",
    "GSE135230": "Heart",   "GSE175625": "Mammary\n(PyVT*)",
}
FVB_SRRS = {
    "GSE123875": ["SRR8321484","SRR8321485","SRR8321486","SRR8321487","SRR8321488","SRR8321489"],
    "GSE123893": ["SRR8321791","SRR8321792","SRR8321793","SRR8321794","SRR8321795","SRR8321796"],
    "GSE123894": ["SRR8322104","SRR8322105","SRR8322106","SRR8322107","SRR8322108","SRR8322109"],
    "GSE135230": ["SRR9879444","SRR9879445","SRR9879446","SRR9879447","SRR9879448",
                  "SRR9879449","SRR9879450","SRR9879451","SRR9879452","SRR9879453"],
    "GSE175625": ["SRR14664845","SRR14664846","SRR14664847","SRR14664848",
                  "SRR14664849","SRR14664850","SRR14664851"],
}
B6_SRRS = {
    "GSE123875": ["SRR8321439","SRR8321440","SRR8321441","SRR8321442","SRR8321443","SRR8321444"],
    "GSE123893": ["SRR8321732","SRR8321733","SRR8321734","SRR8321735","SRR8321736","SRR8321737"],
    "GSE123894": ["SRR8322030","SRR8322031","SRR8322032","SRR8322033","SRR8322034","SRR8322035"],
    "GSE135230": ["SRR9879454","SRR9879455","SRR9879456","SRR9879457","SRR9879458",
                  "SRR9879459","SRR9879460","SRR9879461","SRR9879462","SRR9879463","SRR9879464","SRR9879465"],
    "GSE175625": ["SRR14664793","SRR14664794","SRR14664795","SRR14664796","SRR14664797","SRR14664798"],
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def _rename_srr(df):
    rename = {}
    for col in df.columns:
        m = re.search(r'(SRR\d+)', os.path.basename(col))
        if m: rename[col] = m.group(1)
    df.rename(columns=rename, inplace=True)

def load_counts(path, srrs):
    df = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    df = df.drop(columns=["Chr","Start","End","Strand","Length"], errors="ignore")
    _rename_srr(df)
    avail = [s for s in srrs if s in df.columns]
    sub = df[avail].astype(float)
    lib = sub.sum(axis=0)
    return sub.divide(lib, axis=1) * 1e6

def load_raw(path, srrs):
    """Return raw integer counts (no normalization)."""
    df = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    df = df.drop(columns=["Chr","Start","End","Strand","Length"], errors="ignore")
    _rename_srr(df)
    avail = [s for s in srrs if s in df.columns]
    return df[avail].astype(float)

# ── 1. SNP positions ──────────────────────────────────────────────────────────
print("Fetching SNPs...")
tbx = pysam.TabixFile(VCF)
snp_pos = []
for rec in tbx.fetch(CHROM, GSTART-1, GEND):
    fields = rec.split("\t")
    if len(fields) >= 2:
        snp_pos.append(int(fields[1]))
tbx.close()
print(f"  {len(snp_pos)} SNPs in Fau region")

# ── 2. Per-sample CPM for Fau — GSE135230 heart (highest-bias tissue) ────────
print("Loading counts for GSE135230 (heart)...")
gse_rep = "GSE135230"
arm_cpms = {}

# naive FVB library size (per sample) — shared denominator for fvb arm so that
# unbiased genes land on the diagonal and biased genes deviate upward.
naive_fvb_raw = load_raw(f"{COUNTS}/{gse_rep}.naive.counts.txt", FVB_SRRS[gse_rep])
naive_fvb_lib = naive_fvb_raw.sum(axis=0)   # shape: (n_FVB_samples,)

# naive CPM for B6 (and FVB naive arm)
naive_cpm_all = load_counts(f"{COUNTS}/{gse_rep}.naive.counts.txt",
                             FVB_SRRS[gse_rep] + B6_SRRS[gse_rep])
b6_naive_cpm_fau = (naive_cpm_all.loc[GENE_ID, B6_SRRS[gse_rep]]
                    if GENE_ID in naive_cpm_all.index else pd.Series(dtype=float))

for arm, suffix in [("naive","naive"), ("wasp","wasp"), ("fvb","fvb_b6gtf")]:
    path = f"{COUNTS}/{gse_rep}.{suffix}.counts.txt"
    if arm == "fvb":
        # FVB: load raw counts from FVB-ref BAMs counted with B6 GTF;
        # normalize by naive FVB library size for shared-denominator CPM.
        # B6: same alignment as naive arm → reuse naive B6 CPM.
        fvb_raw_df = load_raw(path, FVB_SRRS[gse_rep])
        avail = fvb_raw_df.columns.intersection(naive_fvb_lib.index)
        fvb_only_cpm = fvb_raw_df[avail].divide(naive_fvb_lib[avail], axis=1) * 1e6
        fvb_vals = (fvb_only_cpm.loc[GENE_ID] if GENE_ID in fvb_only_cpm.index
                    else pd.Series(dtype=float))
        arm_cpms[arm] = pd.concat([fvb_vals, b6_naive_cpm_fau])
    else:
        cpm = load_counts(path, FVB_SRRS[gse_rep] + B6_SRRS[gse_rep])
        arm_cpms[arm] = (cpm.loc[GENE_ID] if GENE_ID in cpm.index
                         else pd.Series(dtype=float))

# ── 3. log2FC per arm per dataset from all_datasets_bias.tsv ──────────────────
print("Loading bias table...")
bias_df = pd.read_csv(f"{PROJ}/results/stage6/all_datasets_bias.tsv", sep="\t")
fau_bias = bias_df[bias_df["gene_id"] == GENE_ID].set_index("gse")

# ── 4. Gene-level SNP density vs genome-wide for context ──────────────────────
meta_path = f"{PROJ}/results/stage6/gene_meta_with_predictors.tsv"
if os.path.exists(meta_path):
    meta = pd.read_csv(meta_path, sep="\t")
    meta_ok = meta.dropna(subset=["snp_density","bias"])
else:
    meta_ok = None

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 9))
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(2, 2, wspace=0.36, hspace=0.42,
                      left=0.08, right=0.97, top=0.88, bottom=0.09)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

panel_labels = ["A","B","C","D"]
for ax, lbl in zip([axA,axB,axC,axD], panel_labels):
    ax.text(-0.14, 1.05, lbl, transform=ax.transAxes,
            fontsize=15, fontweight="bold")

# ── Panel A: Gene model + SNPs ────────────────────────────────────────────────
axA.set_xlim(GSTART - 50, GEND + 50)
axA.set_ylim(-1.8, 1.5)
axA.axis("off")

# Intron line
axA.plot([GSTART, GEND], [0, 0], color="#555555", lw=1.5, zorder=1)
# Strand arrow
for pos in np.linspace(GSTART+200, GEND-200, 6):
    axA.annotate("", xy=(pos+40, 0), xytext=(pos, 0),
                 arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.6))

# Exons
for i,(es,ee) in enumerate(EXONS):
    rect = mpatches.FancyBboxPatch(
        (es, -0.35), ee-es, 0.7,
        boxstyle="square,pad=0", facecolor=BLUE, edgecolor="white", lw=0.5, zorder=3
    )
    axA.add_patch(rect)
    # exon number
    axA.text((es+ee)/2, 0, str(i+1), ha="center", va="center",
             fontsize=8, color="white", fontweight="bold", zorder=4)

# SNP ticks
for pos in snp_pos:
    axA.plot([pos, pos], [-0.6, -0.95], color=RED, lw=0.8, alpha=0.8, zorder=2)
axA.text((GSTART+GEND)/2, -1.2, f"{len(snp_pos)} FVB SNPs",
         ha="center", va="top", fontsize=8.5, color=RED)

# Scale bar 500 bp
sb_x = GEND - 400
axA.plot([sb_x, sb_x+500], [-1.55, -1.55], color="black", lw=2)
axA.text(sb_x+250, -1.72, "500 bp", ha="center", fontsize=8)

axA.set_title(f"Fau gene structure (chr19:{GSTART:,}–{GEND:,}, {STRAND} strand)\n"
              f"Canonical transcript Fau-201  |  5 exons  |  {GEND-GSTART+1:,} bp",
              fontsize=9.5, pad=6)

# ── Panel B: Per-sample CPM (GSE123875 adipose) ───────────────────────────────
arm_order  = ["naive", "wasp", "fvb"]
arm_colors = {"naive": BLUE, "wasp": ORANGE, "fvb": GREEN}
arm_label  = {"naive": "Naive", "wasp": "WASP", "fvb": "FVB ref"}
x_pos = {"naive": 0, "wasp": 3, "fvb": 6}

for arm in arm_order:
    vals = arm_cpms.get(arm, pd.Series(dtype=float))
    if len(vals) == 0:
        continue
    fvb_ids = [s for s in FVB_SRRS[gse_rep] if s in vals.index]
    b6_ids  = [s for s in B6_SRRS[gse_rep]  if s in vals.index]
    # log₂(CPM+1) so differences are legible on a common scale
    fvb_log = np.log2(vals[fvb_ids].values + 1)
    b6_log  = np.log2(vals[b6_ids].values  + 1)
    xp = x_pos[arm]
    col = arm_colors[arm]

    # FVB dots + mean line
    jx_fvb = xp + np.random.default_rng(1).uniform(-0.22, 0.22, len(fvb_log))
    axB.scatter(jx_fvb, fvb_log, s=40, color=col, alpha=0.85,
                edgecolors="white", lw=0.5, zorder=4)
    axB.plot([xp-0.35, xp+0.35], [fvb_log.mean()]*2,
             color=col, lw=2.0, zorder=5)

    # B6 dots + mean line (hollow)
    jx_b6 = xp + 1.0 + np.random.default_rng(2).uniform(-0.22, 0.22, len(b6_log))
    axB.scatter(jx_b6, b6_log, s=40, facecolors="white",
                edgecolors=col, lw=1.2, zorder=4)
    axB.plot([xp+0.65, xp+1.35], [b6_log.mean()]*2,
             color=col, lw=2.0, ls="--", zorder=5)

# Arm labels as x-ticks
axB.set_xticks([0.5, 3.5, 6.5])
axB.set_xticklabels(
    [f"$\\bf{{Naive}}$", f"$\\bf{{WASP}}$", f"$\\bf{{FVB\\ ref}}$"],
    fontsize=9,
    color="black",
)
for tick, col in zip(axB.get_xticklabels(), [BLUE, ORANGE, GREEN]):
    tick.set_color(col)

# Legend for strain
fvb_patch = Line2D([0],[0], marker='o', color='w', markerfacecolor=GRAY,
                   markersize=8, label="FVB (filled)")
b6_patch  = Line2D([0],[0], marker='o', color='w', markerfacecolor='w',
                   markeredgecolor=GRAY, markersize=8, label="B6 (open)")
axB.legend(handles=[fvb_patch, b6_patch], fontsize=8, loc="upper right",
           frameon=True, framealpha=0.9)

axB.set_xlim(-0.8, 8.3)
axB.set_ylabel("log₂(CPM + 1)", fontsize=10)
axB.spines[["top","right"]].set_visible(False)
axB.set_title("Per-sample Fau expression — GSE135230 (heart)\n"
              "unique-only pipeline; shared naive library size; bar = group mean",
              fontsize=10, pad=6)
all_log_vals = [np.log2(v.values + 1) for v in arm_cpms.values() if len(v) > 0]
ymax_log = max(v.max() for v in all_log_vals)
axB.set_ylim(0, ymax_log * 1.15)

# ── Panel C: log2FC per arm × dataset ─────────────────────────────────────────
x     = np.arange(len(GSES))
width = 0.25
rects_n = axC.bar(x - width, [fau_bias.loc[g,"log2fc_naive"] if g in fau_bias.index else np.nan
                               for g in GSES],
                  width, color=BLUE,   alpha=0.8, label="Naive", edgecolor="white")
rects_w = axC.bar(x,         [fau_bias.loc[g,"log2fc_wasp"]  if g in fau_bias.index else np.nan
                               for g in GSES],
                  width, color=ORANGE, alpha=0.8, label="WASP",  edgecolor="white")
rects_f = axC.bar(x + width, [fau_bias.loc[g,"log2fc_fvb"]   if g in fau_bias.index else np.nan
                               for g in GSES],
                  width, color=GREEN,  alpha=0.8, label="FVB ref", edgecolor="white")

axC.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
axC.set_xticks(x)
axC.set_xticklabels([TISSUES[g] for g in GSES], fontsize=9)
axC.set_ylabel("log₂FC (FVB / B6)", fontsize=10)
axC.spines[["top","right"]].set_visible(False)
axC.legend(fontsize=8.5, loc="lower right", frameon=True,
           framealpha=0.9, edgecolor="#CCCCCC")
axC.set_title(f"Fau log₂FC per arm across 5 datasets (unique-only)\n"
              f"(4-tissue mean bias = +{fau_bias.loc[fau_bias['tissue']!='mammary_PyVT','bias'].mean():.2f}; heart most affected)",
              fontsize=10, pad=6)
axC.text(len(GSES)-0.5+0.1, axC.get_ylim()[0]*0.95, "* CNVs (PyVT)",
         fontsize=8, color=GRAY, ha="right", style="italic")

# ── Panel D: SNP density vs bias — genome-wide context ────────────────────────
if meta_ok is not None:
    plot_d = meta_ok[meta_ok["snp_density"] > 0].copy()
    xd = plot_d["snp_density"].clip(upper=plot_d["snp_density"].quantile(0.99))
    yd = plot_d["bias"]
    axD.scatter(xd, yd, s=3, alpha=0.10, color=GRAY,
                edgecolors="none", rasterized=True, label="All genes")
    # Fau point
    fau_row = meta_ok[meta_ok["gene_id"] == GENE_ID]
    if len(fau_row) == 1:
        fx = float(fau_row["snp_density"].values[0])
        fy = float(fau_row["bias"].values[0])
        axD.scatter(fx, fy, s=120, color=RED, zorder=6,
                    edgecolors="white", lw=1.2, label="Fau")
        axD.annotate(f"Fau\n({fx:.1f} SNPs/kb\nbias={fy:.2f})",
                     xy=(fx, fy), xytext=(fx+1.5, fy+0.3),
                     fontsize=8, color=RED,
                     arrowprops=dict(arrowstyle="-", color=RED, lw=0.8))

    # Rsph3a point (most robustly evidenced catalog gene)
    rsph3a_row = meta_ok[meta_ok["gene_id"] == "ENSMUSG00000073471"]
    if len(rsph3a_row) == 1:
        rx = float(rsph3a_row["snp_density"].values[0])
        ry = float(rsph3a_row["bias"].values[0])
        axD.scatter(rx, ry, s=120, color=GREEN, zorder=6,
                    edgecolors="white", lw=1.2, marker="*", label="Rsph3a (most robust)")
        axD.annotate(f"Rsph3a\n({rx:.1f} SNPs/kb\nbias={ry:.2f})",
                     xy=(rx, ry), xytext=(rx+1.2, ry-0.45),
                     fontsize=8, color=GREEN,
                     arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.8))

    # Trend line (binned)
    bins = np.percentile(xd, np.linspace(0,100,21))
    bins = np.unique(bins)
    bidx = np.digitize(xd, bins)
    bx2, by2 = [], []
    for b in range(1, len(bins)):
        sel = yd[bidx == b]
        if len(sel) > 20:
            bx2.append(bins[b-1]+(bins[b]-bins[b-1])/2)
            by2.append(sel.mean())
    axD.plot(bx2, by2, "o-", color=RED, ms=4, lw=1.5, alpha=0.7, label="Binned mean")

    axD.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    axD.set_xlabel("SNP density (SNPs/kb)", fontsize=10)
    axD.set_ylabel("Mean bias across datasets", fontsize=10)
    axD.spines[["top","right"]].set_visible(False)
    axD.legend(fontsize=8.5, loc="upper left", frameon=True,
               framealpha=0.9, edgecolor="#CCCCCC")
    axD.set_title("Fau and Rsph3a in genome-wide context\n(both lie among SNP-dense biased loci)",
                  fontsize=10, pad=6)
else:
    axD.text(0.5, 0.5, "gene_meta_with_predictors.tsv\nnot found",
             transform=axD.transAxes, ha="center", va="center", color=GRAY)
    axD.axis("off")

n_biased = (fau_bias["bias"] >= 0.5).sum()
max_bias = fau_bias["bias"].max()
fig.suptitle(
    f"Worked examples: Fau and Rsph3a — SNP-mediated reference-genome mapping bias\n"
    f"Fau (chr19): heart bias = +{max_bias:.2f} log₂FC under unique-only alignment (23 FVB SNPs; "
    f"SNP-induced multi-mapping class); Rsph3a: most robustly evidenced catalog gene (n=4 non-tumour tissues)",
    fontsize=9.5, y=0.995, color="#333333"
)

for ext in ["png","pdf"]:
    p = f"{OUT}/fig6_worked_example.{ext}"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {p}")
plt.close(fig)
print("Fig 6 done.")
