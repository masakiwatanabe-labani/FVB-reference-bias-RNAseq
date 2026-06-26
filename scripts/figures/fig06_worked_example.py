#!/usr/bin/env python3
"""
Fig 6: Worked examples — Rsph3a (robust/reproducible) and Fau (extreme/multi-mapping)

Panel A: Rsph3a log2FC per arm × 4 non-tumour tissues  ← left column
Panel B: Fau per-sample CPM (heart) — most impactful    ← right column
Panel C: SNP density × bias genome-wide context (both)  ← left column
Panel D: Fau log2FC per arm × 5 datasets                ← right column
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import pysam, re, os

PROJ = os.environ.get("FVB_PROJ", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
OUT    = f"{PROJ}/results/figures"
VCF    = f"{PROJ}/ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf.gz"
COUNTS = f"{PROJ}/results/strain_comp"
os.makedirs(OUT, exist_ok=True)

BLUE   = "#2166AC"
ORANGE = "#E08214"
GREEN  = "#1A9850"
GRAY   = "#888888"
RED    = "#CC0000"

# ── Gene identifiers ──────────────────────────────────────────────────────────
RSPH3A_ID   = "ENSMUSG00000073471"
RSPH3A_NAME = "Rsph3a"

FAU_ID   = "ENSMUSG00000038274"
FAU_NAME = "Fau"
FAU_GSE  = "GSE135230"   # heart — highest-bias tissue

GSES = ["GSE123875","GSE123893","GSE123894","GSE135230","GSE175625"]
GSES_4T = ["GSE123875","GSE123893","GSE123894","GSE135230"]
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
    df = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    df = df.drop(columns=["Chr","Start","End","Strand","Length"], errors="ignore")
    _rename_srr(df)
    avail = [s for s in srrs if s in df.columns]
    return df[avail].astype(float)

# ── Load bias table ───────────────────────────────────────────────────────────
print("Loading bias table...")
bias_df = pd.read_csv(f"{PROJ}/results/stage6/all_datasets_bias.tsv", sep="\t")
rsph3a_bias = bias_df[bias_df["gene_id"] == RSPH3A_ID].set_index("gse")
fau_bias    = bias_df[bias_df["gene_id"] == FAU_ID].set_index("gse")

# ── Load Fau CPM for heart (Panel B) ──────────────────────────────────────────
print("Loading Fau CPM (GSE135230 heart)...")
naive_fvb_raw = load_raw(f"{COUNTS}/{FAU_GSE}.naive.counts.txt", FVB_SRRS[FAU_GSE])
naive_fvb_lib = naive_fvb_raw.sum(axis=0)

naive_cpm_all = load_counts(f"{COUNTS}/{FAU_GSE}.naive.counts.txt",
                             FVB_SRRS[FAU_GSE] + B6_SRRS[FAU_GSE])
b6_naive_cpm_fau = (naive_cpm_all.loc[FAU_ID, B6_SRRS[FAU_GSE]]
                    if FAU_ID in naive_cpm_all.index else pd.Series(dtype=float))

arm_cpms = {}
for arm, suffix in [("naive","naive"), ("wasp","wasp"), ("fvb","fvb_b6gtf")]:
    path = f"{COUNTS}/{FAU_GSE}.{suffix}.counts.txt"
    if arm == "fvb":
        fvb_raw_df = load_raw(path, FVB_SRRS[FAU_GSE])
        avail = fvb_raw_df.columns.intersection(naive_fvb_lib.index)
        fvb_only_cpm = fvb_raw_df[avail].divide(naive_fvb_lib[avail], axis=1) * 1e6
        fvb_vals = (fvb_only_cpm.loc[FAU_ID] if FAU_ID in fvb_only_cpm.index
                    else pd.Series(dtype=float))
        arm_cpms[arm] = pd.concat([fvb_vals, b6_naive_cpm_fau])
    else:
        cpm = load_counts(path, FVB_SRRS[FAU_GSE] + B6_SRRS[FAU_GSE])
        arm_cpms[arm] = (cpm.loc[FAU_ID] if FAU_ID in cpm.index
                         else pd.Series(dtype=float))

# ── Load SNP density context ───────────────────────────────────────────────────
meta_path = f"{PROJ}/results/stage6/gene_meta_with_predictors.tsv"
meta_ok = None
if os.path.exists(meta_path):
    meta = pd.read_csv(meta_path, sep="\t")
    meta_ok = meta.dropna(subset=["snp_density","bias"])

# ── Figure layout: 2×2 (left=Rsph3a, right=Fau) ──────────────────────────────
fig = plt.figure(figsize=(13, 9))
fig.patch.set_facecolor("white")
gs = fig.add_gridspec(2, 2, wspace=0.36, hspace=0.46,
                      left=0.08, right=0.97, top=0.88, bottom=0.09)
axA = fig.add_subplot(gs[0, 0])   # Rsph3a 4-tissue log2FC
axB = fig.add_subplot(gs[0, 1])   # Fau heart CPM
axC = fig.add_subplot(gs[1, 0])   # SNP density context
axD = fig.add_subplot(gs[1, 1])   # Fau 5-dataset log2FC

for ax, lbl in zip([axA,axB,axC,axD], ["A","B","C","D"]):
    ax.text(-0.14, 1.05, lbl, transform=ax.transAxes, fontsize=15, fontweight="bold")

# ── Panel A: Rsph3a log2FC per arm × 4 non-tumour tissues ────────────────────
x = np.arange(len(GSES_4T))
width = 0.25
arms_spec = [("log2fc_naive", BLUE,   "Naive"),
             ("log2fc_wasp",  ORANGE, "WASP"),
             ("log2fc_fvb",   GREEN,  "FVB ref")]
for i, (col_name, col, label) in enumerate(arms_spec):
    vals = [rsph3a_bias.loc[g, col_name] if g in rsph3a_bias.index else np.nan
            for g in GSES_4T]
    axA.bar(x + (i-1)*width, vals, width, color=col, alpha=0.85,
            label=label, edgecolor="white", zorder=3)

axA.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
axA.axhline(0.5, color=RED, lw=0.7, ls=":", alpha=0.6)
axA.text(len(GSES_4T)-0.1, 0.52, "catalog\nthreshold", fontsize=8.5,
         color=RED, ha="right", va="bottom")
axA.set_xticks(x)
axA.set_xticklabels([TISSUES[g] for g in GSES_4T], fontsize=10.5)
axA.set_ylabel("log₂FC (FVB / B6)", fontsize=11)
axA.set_ylim(-0.2, axA.get_ylim()[1])
axA.spines[["top","right"]].set_visible(False)
axA.legend(fontsize=9.5, loc="upper center",
           bbox_to_anchor=(0.5, -0.14), bbox_transform=axA.transAxes,
           frameon=True, framealpha=0.9, edgecolor="#CCCCCC", ncol=3)
rsph3a_mean = rsph3a_bias.loc[GSES_4T, "bias"].mean() if set(GSES_4T) <= set(rsph3a_bias.index) else float("nan")
axA.set_title(f"Rsph3a — 4-tissue mean bias = +{rsph3a_mean:.2f} log₂FC\n"
              f"Most robustly evidenced catalog gene (4/4 non-tumour tissues)",
              fontsize=11, pad=6)

# ── Panel B: Fau per-sample CPM (heart) ───────────────────────────────────────
arm_order  = ["naive", "wasp", "fvb"]
arm_colors = {"naive": BLUE, "wasp": ORANGE, "fvb": GREEN}
x_pos      = {"naive": 0, "wasp": 3, "fvb": 6}

for arm in arm_order:
    vals = arm_cpms.get(arm, pd.Series(dtype=float))
    if len(vals) == 0:
        continue
    fvb_ids = [s for s in FVB_SRRS[FAU_GSE] if s in vals.index]
    b6_ids  = [s for s in B6_SRRS[FAU_GSE]  if s in vals.index]
    fvb_log = np.log2(vals[fvb_ids].values + 1)
    b6_log  = np.log2(vals[b6_ids].values  + 1)
    xp  = x_pos[arm]
    col = arm_colors[arm]
    jx_fvb = xp + np.random.default_rng(1).uniform(-0.22, 0.22, len(fvb_log))
    axB.scatter(jx_fvb, fvb_log, s=40, color=col, alpha=0.85,
                edgecolors="white", lw=0.5, zorder=4)
    axB.plot([xp-0.35, xp+0.35], [fvb_log.mean()]*2, color=col, lw=2.0, zorder=5)
    jx_b6 = xp + 1.0 + np.random.default_rng(2).uniform(-0.22, 0.22, len(b6_log))
    axB.scatter(jx_b6, b6_log, s=40, facecolors="white",
                edgecolors=col, lw=1.2, zorder=4)
    axB.plot([xp+0.65, xp+1.35], [b6_log.mean()]*2, color=col, lw=2.0, ls="--", zorder=5)

axB.set_xticks([0.5, 3.5, 6.5])
axB.set_xticklabels(["Naive", "WASP", "FVB ref"], fontsize=10)
for tick, col in zip(axB.get_xticklabels(), [BLUE, ORANGE, GREEN]):
    tick.set_color(col)
fvb_patch = Line2D([0],[0], marker='o', color='w', markerfacecolor=GRAY,
                   markersize=8, label="FVB (filled)")
b6_patch  = Line2D([0],[0], marker='o', color='w', markerfacecolor='w',
                   markeredgecolor=GRAY, markersize=8, label="B6 (open)")
axB.legend(handles=[fvb_patch, b6_patch], fontsize=9, loc="upper center",
           bbox_to_anchor=(0.5, -0.14), bbox_transform=axB.transAxes,
           frameon=True, framealpha=0.9, edgecolor="#CCCCCC", ncol=2)
axB.set_xlim(-0.8, 8.3)
axB.set_ylabel("log₂(CPM + 1)", fontsize=11)
axB.spines[["top","right"]].set_visible(False)
fau_heart_bias = fau_bias.loc[FAU_GSE, "bias"] if FAU_GSE in fau_bias.index else float("nan")
axB.set_title(f"Fau — GSE135230 heart, per-sample expression\n"
              f"Heart bias = +{fau_heart_bias:.2f} log₂FC (unique-only); bar = group mean",
              fontsize=11, pad=6)
all_log_vals = [np.log2(v.values + 1) for v in arm_cpms.values() if len(v) > 0]
axB.set_ylim(0, max(v.max() for v in all_log_vals) * 1.15)

# ── Panel C: SNP density × bias — genome-wide context ─────────────────────────
if meta_ok is not None:
    plot_d = meta_ok[meta_ok["snp_density"] > 0].copy()
    sd_cap = plot_d["snp_density"].quantile(0.99)
    xd = plot_d["snp_density"].clip(upper=sd_cap)
    yd = plot_d["bias"]
    axC.scatter(xd, yd, s=3, alpha=0.10, color=GRAY,
                edgecolors="none", rasterized=True, label="All genes")

    # Binned trend
    bins = np.unique(np.percentile(xd, np.linspace(0,100,21)))
    bidx = np.digitize(xd, bins)
    bx2, by2 = [], []
    for b in range(1, len(bins)):
        sel = yd[bidx == b]
        if len(sel) > 20:
            bx2.append(bins[b-1]+(bins[b]-bins[b-1])/2)
            by2.append(sel.mean())
    axC.plot(bx2, by2, "o-", color="#666666", ms=4, lw=1.5, alpha=0.7, label="Binned mean")

    # Rsph3a (green star)
    rsph3a_row = meta_ok[meta_ok["gene_id"] == RSPH3A_ID]
    if len(rsph3a_row) == 1:
        rx = min(float(rsph3a_row["snp_density"].values[0]), sd_cap)
        ry = float(rsph3a_row["bias"].values[0])
        axC.scatter(rx, ry, s=180, color=GREEN, zorder=7,
                    edgecolors="white", lw=1.2, marker="*", label="Rsph3a (4/4 tissues)")
        axC.annotate("Rsph3a",
                     xy=(rx, ry), xytext=(rx+1.8, ry+0.55),
                     fontsize=11, color=GREEN, fontweight="bold", zorder=7,
                     arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.8),
                     bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                               edgecolor=GREEN, alpha=0.85, lw=0.8))

    # Fau (red circle)
    fau_row = meta_ok[meta_ok["gene_id"] == FAU_ID]
    if len(fau_row) == 1:
        fx = min(float(fau_row["snp_density"].values[0]), sd_cap)
        fy = float(fau_row["bias"].values[0])
        axC.scatter(fx, fy, s=100, color=RED, zorder=6,
                    edgecolors="white", lw=1.2, label="Fau (heart +3.67)")
        axC.annotate("Fau",
                     xy=(fx, fy), xytext=(fx-4.0, fy+0.55),
                     fontsize=11, color=RED, fontweight="bold",
                     arrowprops=dict(arrowstyle="-", color=RED, lw=0.8),
                     bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                               edgecolor=RED, alpha=0.85, lw=0.8))

    axC.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
    axC.axhline(0.5, color=RED, lw=0.7, ls=":", alpha=0.5)
    axC.set_xlabel("SNP density (SNPs/kb)", fontsize=11)
    axC.set_ylabel("Mean bias across datasets (log₂FC)", fontsize=11)
    axC.spines[["top","right"]].set_visible(False)
    axC.legend(fontsize=9.5, loc="lower right", frameon=True,
               framealpha=0.9, edgecolor="#CCCCCC")
    axC.set_title("Genome-wide context: SNP density × bias\n"
                  "Both Rsph3a and Fau lie among SNP-dense loci",
                  fontsize=11, pad=6)
else:
    axC.text(0.5, 0.5, "gene_meta_with_predictors.tsv\nnot found",
             transform=axC.transAxes, ha="center", va="center", color=GRAY)
    axC.axis("off")

# ── Panel D: Fau log2FC per arm × 5 datasets ─────────────────────────────────
x2    = np.arange(len(GSES))
for i, (col_name, col, label) in enumerate(arms_spec):
    vals = [fau_bias.loc[g, col_name] if g in fau_bias.index else np.nan for g in GSES]
    axD.bar(x2 + (i-1)*width, vals, width, color=col, alpha=0.85,
            label=label, edgecolor="white", zorder=3)

axD.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
axD.set_xticks(x2)
axD.set_xticklabels([TISSUES[g] for g in GSES], fontsize=10)
axD.set_ylabel("log₂FC (FVB / B6)", fontsize=11)
axD.set_ylim(min(-0.7, axD.get_ylim()[0]), axD.get_ylim()[1])
axD.spines[["top","right"]].set_visible(False)
axD.legend(fontsize=9.5, loc="upper center",
           bbox_to_anchor=(0.5, -0.14), bbox_transform=axD.transAxes,
           frameon=True, framealpha=0.9, edgecolor="#CCCCCC", ncol=3)
fau_4t_mean = fau_bias.loc[GSES_4T, "bias"].mean() if set(GSES_4T) <= set(fau_bias.index) else float("nan")
axD.set_title(f"Fau log₂FC per arm — all 5 datasets (unique-only)\n"
              f"4-tissue mean bias = +{fau_4t_mean:.2f} log₂FC; heart most affected",
              fontsize=11, pad=6)
axD.text(len(GSES)-0.5, axD.get_ylim()[0]*0.95 if axD.get_ylim()[0] < 0 else -0.2,
         "* CNVs (PyVT)", fontsize=9, color=GRAY, ha="right", style="italic")

# ── Suptitle ──────────────────────────────────────────────────────────────────
fig.suptitle(
    "Worked examples of SNP-mediated reference-genome mapping bias\n"
    "Left: Rsph3a — most robustly evidenced catalog gene (4/4 non-tumour tissues; methodologically clean)   "
    "Right: Fau — largest bias in catalog (heart +{:.2f} log₂FC; SNP-induced multi-mapping class)".format(fau_heart_bias),
    fontsize=10.5, y=0.995, color="#333333"
)

for ext in ["png","pdf"]:
    p = f"{OUT}/fig6_worked_example.{ext}"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {p}")
plt.close(fig)
print("Fig 6 done.")
