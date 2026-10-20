#!/usr/bin/env python3
"""Shared helpers for the v2 figures (revision_v2/scripts/fig*.py, sfig*.py).
Inputs come only from analysis_v2/results (never legacy/): bias/, catalog/, predictors/, counts/, f1_ase/, survey/, gene_meta.tsv.
Environment: RES_DIR  -> analysis results dir (default analysis_v2/results; per-tissue dry runs point it to scratch)
             FIG_DIR  -> output dir (default analysis_v2/figures)
             DESIGN   -> edgeR design key, "primary" (default; ~ diet + strain where two diets exist) or "water"
Every figure is written as vector PDF + 600-dpi LZW TIFF, and its plotted numbers as <name>_data.tsv next to it."""
import os, matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"]=42; matplotlib.rcParams["ps.fonttype"]=42; matplotlib.rcParams["font.size"]=9
import matplotlib.pyplot as plt, numpy as np, pandas as pd
PROJ="/usr/local/jupyter/FVB_bias"
FINAL=f"{PROJ}/analysis_v2/results"
RES=os.environ.get("RES_DIR", FINAL)
FIG=os.environ.get("FIG_DIR", f"{PROJ}/analysis_v2/figures"); os.makedirs(FIG, exist_ok=True)
DESIGN=os.environ.get("DESIGN", "primary")
BIAS_THR=0.5; MIN_TISSUES=2; DELTA_THR=0.3          # fixed catalog / mechanism criteria (params.yaml, rule 0.4-8)

BLUE="#2166AC"; LBLUE="#92C5DE"; ORANGE="#E08214"; LORANGE="#FDB863"; GREEN="#1A9850"; GRAY="#888888"; LGRAY="#CCCCCC"; RED="#CC0000"; PURPLE="#6A3D9A"
TISSUE_ORDER=["adipose","hypothalamus","liver","skin"]
TISSUE_LABEL={"adipose":"Adipose","hypothalamus":"Hypothalamus","liver":"Liver","skin":"Skin","mammary_tumour_PyVT":"Mammary tumour\n(PyVT)"}
TISSUE_SHORT={"adipose":"Adi","hypothalamus":"Hyp","liver":"Liv","skin":"Skn","mammary_tumour_PyVT":"PyVT"}
TISSUE_COLOR={"adipose":"#1F78B4","hypothalamus":"#6A3D9A","liver":"#33A02C","skin":"#FF7F00","mammary_tumour_PyVT":"#999999"}
ARM_ORDER=["naive","wasp","fvb"]
ARM_COLOR={"naive":BLUE,"wasp":ORANGE,"fvb":GREEN}
ARM_LABEL={"naive":"Naive (B6 reference)","wasp":"WASP-filtered","fvb":"FVB/NJ reference"}
ARM_SHORT={"naive":"Naive","wasp":"WASP","fvb":"FVB ref"}

def load_bias(design=DESIGN):
    """long table (one row per gene x tissue) and the tissue list in display order"""
    long=pd.read_csv(f"{RES}/bias/all_tissues.{design}.tsv", sep="\t")
    tissues=[t for t in TISSUE_ORDER if t in set(long["tissue"])]
    return long, tissues

def load_catalog(design=DESIGN):
    p=f"{RES}/catalog/catalog.{design}.tsv"
    cat=pd.read_csv(p, sep="\t", index_col=0) if os.path.exists(p) else pd.DataFrame()
    return cat

def load_predictors(design=DESIGN):
    return pd.read_csv(f"{RES}/predictors/gene_predictors.{design}.tsv", sep="\t", index_col=0)

def gene_meta():
    return pd.read_csv(f"{FINAL}/gene_meta.tsv", sep="\t", index_col=0)

def load_counts(tissue, arm, mode="unique"):
    return pd.read_csv(f"{RES}/counts/{tissue}.{arm}.{mode}.tsv", sep="\t", index_col=0)

def load_samples(tissue):
    ss=pd.read_csv(f"{RES}/counts/{tissue}.samples.tsv", sep="\t")
    ss["strain_short"]=ss["strain"].map(lambda s: "FVB" if s.startswith("FVB") else ("B6" if s.startswith("C57") else s))
    return ss

def n_animals(tissue):
    ss=load_samples(tissue); return int((ss.strain_short=="FVB").sum()), int((ss.strain_short=="B6").sum())

def panel(ax, letter, x=-0.12, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom", ha="right")

def despine(ax):
    ax.spines[["top","right"]].set_visible(False)

def save(fig, name, data=None):
    """PDF (vector; rasterized point clouds keep the file small) + 600-dpi TIFF, plus the plotted data table"""
    fig.savefig(f"{FIG}/{name}.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(f"{FIG}/{name}.tiff", dpi=600, bbox_inches="tight", facecolor="white", pil_kwargs={"compression":"tiff_lzw"})
    if data is not None:
        if isinstance(data, dict):
            with open(f"{FIG}/{name}_data.tsv","w") as o:
                for k,v in data.items():
                    o.write(f"## {k}\n"); v.to_csv(o, sep="\t"); o.write("\n")
        else: data.to_csv(f"{FIG}/{name}_data.tsv", sep="\t")
    plt.close(fig); print(f"saved {FIG}/{name}.{{pdf,tiff}}")

def label_points(ax, xs, ys, names, color="#333333", fontsize=7.5, dx=4, dy=2):
    for x,y,n in zip(xs,ys,names):
        ax.annotate(n, (x,y), xytext=(dx,dy), textcoords="offset points", fontsize=fontsize, color=color,
                    arrowprops=dict(arrowstyle="-", color=color, lw=0.4, shrinkB=2))

def clip_for_display(ax, y, lo=-3, hi=3, q=(0.001,0.999)):
    """clip extreme y values to [min(lo,q0.1%), max(hi,q99.9%)] for the scatter and note how many were clipped"""
    y=pd.Series(np.asarray(y,dtype=float)); a=min(lo,y.quantile(q[0])); b=max(hi,y.quantile(q[1])); n=int(((y<a)|(y>b)).sum())
    if n: ax.text(0.01,0.01,f"{n} gene(s) beyond the axis range shown at the limit",transform=ax.transAxes,fontsize=7,color=GRAY,va="bottom")
    ax.set_ylim(a-0.05*(b-a), b+0.05*(b-a)); return y.clip(a,b).values
