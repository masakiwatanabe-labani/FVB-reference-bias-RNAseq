#!/usr/bin/env python3
"""
Fig 4: Predictors of reference bias
Panel A: Bias vs SNP density (SNPs/kb) — per-gene, all datasets pooled
Panel B: Bias vs gene length (log10 bp)
Panel C: Bias vs expression level (log2 mean CPM)
Panel D: Bias by gene biotype (protein_coding vs lncRNA vs other)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import pysam
import re, os
from scipy import stats

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

# ── 1. Load cached predictor table or compute from scratch ───────────────────
META_PATH = f"{PROJ}/results/stage6/gene_meta_with_predictors.tsv"

if os.path.exists(META_PATH):
    print(f"Loading cached predictor table: {META_PATH}")
    meta = pd.read_csv(META_PATH, sep="\t")
    print(f"  {len(meta):,} genes loaded")
    # ensure derived columns
    if "log10_length" not in meta.columns:
        meta["log10_length"]  = np.log10(meta["length"].clip(lower=1))
    if "log2_mean_cpm" not in meta.columns:
        meta["log2_mean_cpm"] = np.log2(meta["mean_cpm"].clip(lower=0.01))
    if "snp_density_kb" not in meta.columns:
        meta["snp_density_kb"] = meta["snp_density"]

else:
    # ── 1b. Load bias table ───────────────────────────────────────────────────
    print("Loading bias table...")
    df = pd.read_csv(f"{PROJ}/results/stage6/all_datasets_bias.tsv", sep="\t")
    df = df.dropna(subset=["bias"])

    # ── 2. Parse B6 GTF for gene attributes ──────────────────────────────────
    print("Parsing GTF...")
    genes = {}
    with open(GTF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom  = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            attrs  = parts[8]
            gid = re.search(r'gene_id "([^"]+)"', attrs)
            gbt = re.search(r'gene_biotype "([^"]+)"', attrs)
            gnm = re.search(r'gene_name "([^"]+)"', attrs)
            if not gid:
                continue
            gid = gid.group(1)
            genes[gid] = {
                "chrom":   chrom, "start": start, "end": end,
                "length":  end - start + 1,
                "biotype": gbt.group(1) if gbt else "unknown",
                "name":    gnm.group(1) if gnm else gid,
            }
    print(f"  {len(genes):,} genes parsed from GTF")
    gene_df = pd.DataFrame.from_dict(genes, orient="index")
    gene_df.index.name = "gene_id"
    gene_df = gene_df.reset_index()

    # ── 3. SNP density per gene via tabix ────────────────────────────────────
    print("Counting SNPs per gene (tabix)...")
    tbx = pysam.TabixFile(VCF)
    vcf_chroms = set(tbx.contigs)
    snp_counts = {}
    n = len(gene_df)
    for i, row in gene_df.iterrows():
        if i % 5000 == 0:
            print(f"  {i}/{n} ...", flush=True)
        chrom = row["chrom"]
        if chrom not in vcf_chroms:
            chrom2 = chrom.replace("chr","") if chrom.startswith("chr") else f"chr{chrom}"
            if chrom2 not in vcf_chroms:
                snp_counts[row["gene_id"]] = 0
                continue
            chrom = chrom2
        try:
            cnt = sum(1 for _ in tbx.fetch(chrom, row["start"]-1, row["end"]))
        except ValueError:
            cnt = 0
        snp_counts[row["gene_id"]] = cnt
    tbx.close()
    gene_df["snp_count"]   = gene_df["gene_id"].map(snp_counts)
    gene_df["snp_density"] = gene_df["snp_count"] / (gene_df["length"] / 1000.0)
    print(f"  SNP density computed. Mean: {gene_df['snp_density'].mean():.2f} SNPs/kb")

    # ── 4. Expression level ───────────────────────────────────────────────────
    print("Computing mean expression per gene across all datasets...")
    GSES = ["GSE123875","GSE123893","GSE123894","GSE135230","GSE175625"]
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
    expr_list = []
    for gse in GSES:
        path = f"{COUNTS}/{gse}.naive.counts.txt"
        raw = pd.read_csv(path, sep="\t", comment="#", index_col=0)
        raw = raw.drop(columns=["Chr","Start","End","Strand","Length"], errors="ignore")
        rename = {}
        for col in raw.columns:
            m = re.search(r'(SRR\d+)', os.path.basename(col))
            if m: rename[col] = m.group(1)
        raw.rename(columns=rename, inplace=True)
        srrs = [s for s in FVB_SRRS[gse] + B6_SRRS[gse] if s in raw.columns]
        sub = raw[srrs]
        lib = sub.sum(axis=0)
        cpm = sub.divide(lib, axis=1) * 1e6
        expr_list.append(cpm.mean(axis=1).rename(gse))
    gene_mean_cpm = pd.concat(expr_list, axis=1).mean(axis=1)
    gene_mean_cpm.index.name = "gene_id"
    gene_mean_cpm = gene_mean_cpm.reset_index(name="mean_cpm")
    print(f"  Expression computed for {len(gene_mean_cpm):,} genes")

    # ── 5. Merge ──────────────────────────────────────────────────────────────
    bias_gene = df.groupby("gene_id")["bias"].mean().reset_index()
    meta = gene_df.merge(bias_gene, on="gene_id", how="inner")
    meta = meta.merge(gene_mean_cpm, on="gene_id", how="left")
    meta = meta[meta["bias"].notna()].copy()
    meta["log10_length"]   = np.log10(meta["length"].clip(lower=1))
    meta["log2_mean_cpm"]  = np.log2(meta["mean_cpm"].clip(lower=0.01))
    meta["log2_snp_dens"]  = np.log2(meta["snp_density"].clip(lower=0.01))
    meta["snp_density_kb"] = meta["snp_density"]
    print(f"  Merged meta: {len(meta):,} genes with bias + predictors")
    meta.to_csv(META_PATH, sep="\t", index=False)

# ── 6. Biotype grouping ───────────────────────────────────────────────────────
def classify_biotype(bt):
    if bt == "protein_coding":
        return "protein_coding"
    elif bt in ("lncRNA","lincRNA","antisense","sense_intronic","sense_overlapping"):
        return "lncRNA"
    elif bt in ("processed_pseudogene","unprocessed_pseudogene","pseudogene",
                "transcribed_processed_pseudogene","transcribed_unprocessed_pseudogene"):
        return "pseudogene"
    else:
        return "other"

meta["biotype_group"] = meta["biotype"].apply(classify_biotype)

# ── 7. Figure ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
fig.patch.set_facecolor("white")
fig.subplots_adjust(wspace=0.32, hspace=0.38, left=0.08, right=0.97,
                    top=0.90, bottom=0.09)

ax_labels = ["A","B","C","D"]
for ax, lbl in zip(axes.flat, ax_labels):
    ax.text(-0.14, 1.04, lbl, transform=ax.transAxes,
            fontsize=15, fontweight="bold")

# ── Panel A: bias vs SNP density ─────────────────────────────────────────────
axA = axes[0,0]
plotA = meta[meta["snp_density_kb"] > 0].copy()
x_sd = plotA["snp_density_kb"].clip(upper=plotA["snp_density_kb"].quantile(0.995))
axA.scatter(x_sd, plotA["bias"], s=3, alpha=0.12, color=BLUE,
            edgecolors="none", rasterized=True)
# binned trend
bins = np.percentile(x_sd, np.linspace(0,100,21))
bins = np.unique(bins)
bin_idx = np.digitize(x_sd, bins)
bx, by, be = [], [], []
for b in range(1, len(bins)):
    sel = plotA["bias"][bin_idx == b]
    if len(sel) > 10:
        bx.append(bins[b-1] + (bins[b]-bins[b-1])/2)
        by.append(sel.mean())
        be.append(sel.sem())
axA.errorbar(bx, by, yerr=be, fmt="o-", color=RED, ms=5, lw=1.5,
             capsize=3, zorder=5, label="Binned mean ± SE")
axA.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)

# Highlight Rsph3a (most robust catalog gene) and Fau (most extreme)
sd_cap = plotA["snp_density_kb"].quantile(0.995)
for gid, gname, col, offset in [
    ("ENSMUSG00000073471", "Rsph3a", GREEN, (0.5, 0.10)),
    ("ENSMUSG00000038274", "Fau",    RED,   (0.5, 0.10)),
]:
    row = meta[meta["gene_id"] == gid]
    if len(row) == 1:
        xpt = min(float(row["snp_density_kb"].values[0]), sd_cap)
        ypt = float(row["bias"].values[0])
        axA.scatter(xpt, ypt, s=90, color=col, zorder=7,
                    edgecolors="white", lw=1.0, marker="*")
        axA.annotate(gname, xy=(xpt, ypt),
                     xytext=(xpt + offset[0], ypt + offset[1]),
                     fontsize=8, color=col, zorder=7,
                     arrowprops=dict(arrowstyle="-", color=col, lw=0.6))

r, p = stats.spearmanr(x_sd, plotA["bias"])
axA.text(0.97, 0.97, f"ρ = {r:.3f}\np = {p:.1e}",
         transform=axA.transAxes, fontsize=8.5, va="top", ha="right",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                   edgecolor="#CCCCCC", alpha=0.9))
axA.set_xlabel("SNP density (SNPs/kb)", fontsize=10)
axA.set_ylabel("Bias (log₂FC FVB ref − naive)", fontsize=10)
axA.set_title("Bias vs SNP density\n(SNP-dense loci enriched among biased genes)", fontsize=10.5, pad=6)
axA.spines[["top","right"]].set_visible(False)
axA.legend(fontsize=8, loc="lower right")

# ── Panel B: bias vs gene length ─────────────────────────────────────────────
axB = axes[0,1]
axB.scatter(meta["log10_length"], meta["bias"], s=3, alpha=0.12,
            color=ORANGE, edgecolors="none", rasterized=True)
bins_l = np.percentile(meta["log10_length"], np.linspace(0,100,21))
bins_l = np.unique(bins_l)
bidx_l = np.digitize(meta["log10_length"], bins_l)
bxl, byl, bel = [], [], []
for b in range(1, len(bins_l)):
    sel = meta["bias"][bidx_l == b]
    if len(sel) > 10:
        bxl.append(bins_l[b-1] + (bins_l[b]-bins_l[b-1])/2)
        byl.append(sel.mean())
        bel.append(sel.sem())
axB.errorbar(bxl, byl, yerr=bel, fmt="o-", color=RED, ms=5, lw=1.5, capsize=3, zorder=5)
axB.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
r2, p2 = stats.spearmanr(meta["log10_length"], meta["bias"])
axB.text(0.97, 0.97, f"ρ = {r2:.3f}\np = {p2:.1e}",
         transform=axB.transAxes, fontsize=8.5, va="top", ha="right",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                   edgecolor="#CCCCCC", alpha=0.9))
xticks = [3,4,5,6]
axB.set_xticks(xticks)
axB.set_xticklabels([f"10³","10⁴","10⁵","10⁶"], fontsize=9)
axB.set_xlabel("Gene length (bp)", fontsize=10)
axB.set_ylabel("Bias (log₂FC FVB ref − naive)", fontsize=10)
axB.set_title("Bias vs gene length", fontsize=11, pad=6)
axB.spines[["top","right"]].set_visible(False)

# ── Panel C: bias vs expression ───────────────────────────────────────────────
axC = axes[1,0]
plotC = meta[meta["mean_cpm"] > 0.01].copy()
axC.scatter(plotC["log2_mean_cpm"], plotC["bias"], s=3, alpha=0.12,
            color=GREEN, edgecolors="none", rasterized=True)
bins_e = np.percentile(plotC["log2_mean_cpm"], np.linspace(0,100,21))
bins_e = np.unique(bins_e)
bidx_e = np.digitize(plotC["log2_mean_cpm"], bins_e)
bxe, bye, bee = [], [], []
for b in range(1, len(bins_e)):
    sel = plotC["bias"][bidx_e == b]
    if len(sel) > 10:
        bxe.append(bins_e[b-1] + (bins_e[b]-bins_e[b-1])/2)
        bye.append(sel.mean())
        bee.append(sel.sem())
axC.errorbar(bxe, bye, yerr=bee, fmt="o-", color=RED, ms=5, lw=1.5, capsize=3, zorder=5)
axC.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
r3, p3 = stats.spearmanr(plotC["log2_mean_cpm"], plotC["bias"])
axC.text(0.97, 0.97, f"ρ = {r3:.3f}\np = {p3:.1e}",
         transform=axC.transAxes, fontsize=8.5, va="top", ha="right",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                   edgecolor="#CCCCCC", alpha=0.9))
axC.set_xlabel("Expression level (log₂ mean CPM)", fontsize=10)
axC.set_ylabel("Bias (log₂FC FVB ref − naive)", fontsize=10)
axC.set_title("Bias vs expression level", fontsize=11, pad=6)
axC.spines[["top","right"]].set_visible(False)

# ── Panel D: bias by biotype ──────────────────────────────────────────────────
axD = axes[1,1]
biotype_order  = ["protein_coding","lncRNA","pseudogene","other"]
biotype_colors = [BLUE, GREEN, ORANGE, GRAY]
biotype_labels = ["Protein-coding","lncRNA","Pseudogene","Other"]

btype_data = [meta[meta["biotype_group"]==bt]["bias"].values for bt in biotype_order]
btype_ns   = [len(d) for d in btype_data]

vp = axD.violinplot(btype_data, positions=range(len(biotype_order)),
                    widths=0.6, showmedians=False, showextrema=False)
for pc, col in zip(vp["bodies"], biotype_colors):
    pc.set_facecolor(col); pc.set_alpha(0.7); pc.set_edgecolor("white")

# median markers
for i, vals in enumerate(btype_data):
    axD.scatter(i, np.median(vals), s=25, color="white",
                edgecolors=biotype_colors[i], linewidths=1.5, zorder=5)

axD.axhline(0, color="black", lw=0.8, ls="--", alpha=0.5)
axD.set_xticks(range(len(biotype_order)))
axD.set_xticklabels([f"{biotype_labels[i]}\n(n={btype_ns[i]:,})"
                     for i in range(len(biotype_order))], fontsize=8.5)
axD.set_ylabel("Bias (log₂FC FVB ref − naive)", fontsize=10)
axD.set_title("Bias by gene biotype", fontsize=11, pad=6)
axD.spines[["top","right"]].set_visible(False)

# Kruskal-Wallis
kw_stat, kw_p = stats.kruskal(*btype_data)
axD.text(0.97, 0.97, f"Kruskal-Wallis\np = {kw_p:.1e}",
         transform=axD.transAxes, fontsize=8.5, va="top", ha="right",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                   edgecolor="#CCCCCC", alpha=0.9))

fig.suptitle(
    "Multi-factor predictors of reference bias: SNP density, gene biotype, and expression level",
    fontsize=10.5, y=0.995, color="#333333"
)

for ext in ["png","pdf"]:
    p = f"{OUT}/fig4_predictors.{ext}"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {p}")
plt.close(fig)
print("Fig 4 done.")
