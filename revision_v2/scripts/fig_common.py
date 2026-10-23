#!/usr/bin/env python3
"""Shared helpers for the v2 figures (revision_v2/scripts/fig*.py, sfig*.py).
Inputs come only from analysis_v2/results (never legacy/): bias/, catalog/, predictors/, counts/, f1_ase/, survey/, gene_meta.tsv.
Environment: RES_DIR  -> analysis results dir (default analysis_v2/results; per-tissue dry runs point it to scratch)
             FIG_DIR  -> output dir (default analysis_v2/figures)
             DESIGN   -> edgeR design key, "primary" (default; ~ diet + strain where two diets exist) or "water"
Every figure is written as vector PDF + 600-dpi LZW TIFF, and its plotted numbers as <name>_data.tsv next to it."""
import os, matplotlib, matplotlib.legend
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"]=42; matplotlib.rcParams["ps.fonttype"]=42
# Phase N-3: figures are drawn at their final print size (full width 170 mm, single column 85 mm, height <= 225 mm) with all text >= 7 pt
matplotlib.rcParams.update({"font.size":7,"axes.titlesize":8,"axes.labelsize":7.5,"xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":7,
                            "legend.title_fontsize":7,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"lines.linewidth":1.0})
MM=1/25.4; FULL_W=170*MM; HALF_W=85*MM; MAX_H=225*MM; MIN_PT=7.0
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
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")

def despine(ax):
    ax.spines[["top","right"]].set_visible(False)

def overlap_check(fig, name):
    """Phase O-3: (i) no two visible text artists overlap (rendered bounding boxes); (ii) no text drawn with a background box, and no legend frame, overlaps a
    bar or patch of the same axes. Allowed exceptions: an axes title over its own axes area (not checked: titles are only compared with
    other texts), colour-bar labels and colour-bar tick labels (axes created by fig.colorbar), and texts inside the same legend
    (laid out by matplotlib). Overlap = intersection wider and taller than 0.5 px at 100 dpi. Result logged to
    revision_v2/freeze/figure_overlap_check.tsv."""
    import matplotlib.text as mtext, matplotlib.patches as mp
    fig.set_dpi(100); fig.canvas.draw(); rend=fig.canvas.get_renderer()
    cbar_axes={a for a in fig.axes if getattr(a, "_colorbar", None) is not None or a.get_label()=="<colorbar>"}
    leg_of={}
    for a in fig.axes + [fig]:
        for lg in ([a.get_legend()] if hasattr(a, "get_legend") and a.get_legend() else []) + (list(fig.legends) if a is fig else []):
            for t in lg.get_texts() + [lg.get_title()]: leg_of[id(t)]=id(lg)
    T=[]
    for t in fig.findobj(mtext.Text):
        if not t.get_visible() or not t.get_text().strip(): continue
        if t.axes in cbar_axes: continue
        bb=mtext.Text.get_window_extent(t, rend)   # text only (an Annotation's own extent includes its leader line)
        if bb.width <= 0 or bb.height <= 0: continue
        T.append((t, bb))
    out=[]
    for i in range(len(T)):
        for j in range(i+1, len(T)):
            (a, ba), (b, bb)=T[i], T[j]
            if leg_of.get(id(a)) is not None and leg_of.get(id(a))==leg_of.get(id(b)): continue
            w=min(ba.x1, bb.x1)-max(ba.x0, bb.x0); h=min(ba.y1, bb.y1)-max(ba.y0, bb.y0)
            if w > 0.5 and h > 0.5: out.append(f"text/text: '{a.get_text()[:30]}' x '{b.get_text()[:30]}'")
    for t, _ in T:
        bp=t.get_bbox_patch()
        if bp is None or t.axes is None: continue
        bpb=bp.get_window_extent(rend)
        for p in t.axes.patches:
            if p is bp or p is t.axes.patch or not p.get_visible(): continue
            pb=p.get_window_extent(rend); w=min(bpb.x1, pb.x1)-max(bpb.x0, pb.x0); h=min(bpb.y1, pb.y1)-max(bpb.y0, pb.y0)
            if w > 0.5 and h > 0.5: out.append(f"boxed text over patch: '{t.get_text()[:30]}'"); break
    for ax in fig.axes:   # legend frames are boxed text too: a legend may not cover a bar or patch of its own axes
        lg=ax.get_legend()
        if lg is None or not lg.get_visible(): continue
        lb=lg.get_window_extent(rend)
        for p in ax.patches:
            if p is ax.patch or not p.get_visible(): continue
            pb=p.get_window_extent(rend); w=min(lb.x1, pb.x1)-max(lb.x0, pb.x0); h=min(lb.y1, pb.y1)-max(lb.y0, pb.y0)
            if w > 0.5 and h > 0.5: out.append(f"legend over patch: '{(lg.get_texts()[0].get_text() if lg.get_texts() else '')[:30]}'"); break
    log=f"{PROJ}/revision_v2/freeze/figure_overlap_check.tsv"; rows={}
    if os.path.exists(log): rows={l.split("\t")[0]:l for l in open(log).read().splitlines()[1:]}
    rows[name]=f"{name}\t{len(T)}\t{len(out)}\t{' | '.join(out) if out else 'none'}"
    open(log,"w").write("figure\ttexts_checked\toverlaps\tdetail\n"+"\n".join(rows[k] for k in sorted(rows))+"\n")
    return out

def save(fig, name, data=None):
    """PDF (vector) + 600-dpi RGB TIFF on white (no alpha channel, Phase N-2c), plus the plotted data table. The figure is saved at its
    exact size (no tight cropping), so its size is the print size; the size and the smallest text are checked (Phase N-3) and logged to
    revision_v2/freeze/figure_sizes.tsv."""
    import io, matplotlib.text as mtext
    from PIL import Image
    w_mm, h_mm = (x*25.4 for x in fig.get_size_inches())
    fig.canvas.draw()
    sizes=[t.get_fontsize() for t in fig.findobj(mtext.Text) if t.get_visible() and t.get_text().strip() and t.get_window_extent().width > 0]
    ov=overlap_check(fig, name)
    if os.environ.get("OVERLAP_REPORT_ONLY"): print(f"[overlap-check] {name}: {len(ov)}"); [print("   ", x) for x in ov]
    else: assert not ov, f"{name}: {len(ov)} overlap(s): " + "; ".join(ov[:12])
    fmin=min(sizes); assert w_mm <= 170.5 and h_mm <= 225.5, f"{name}: {w_mm:.0f} x {h_mm:.0f} mm exceeds 170 x 225 mm"
    assert fmin >= MIN_PT - 1e-9, f"{name}: smallest text {fmin} pt < {MIN_PT} pt"
    fig.savefig(f"{FIG}/{name}.pdf", facecolor="white")
    buf=io.BytesIO(); fig.savefig(buf, format="png", dpi=600, facecolor="white"); buf.seek(0)
    im=Image.open(buf); bg=Image.new("RGB", im.size, (255,255,255)); bg.paste(im, mask=im.split()[3] if im.mode=="RGBA" else None)
    bg.save(f"{FIG}/{name}.tiff", compression="tiff_lzw", dpi=(600,600))
    log=f"{PROJ}/revision_v2/freeze/figure_sizes.tsv"; rows={}
    if os.path.exists(log): rows={l.split("\t")[0]:l for l in open(log).read().splitlines()[1:]}
    rows[name]=f"{name}\t{w_mm:.1f}\t{h_mm:.1f}\t{fmin:.1f}\t{max(sizes):.1f}"
    open(log,"w").write("figure\twidth_mm\theight_mm\tmin_text_pt\tmax_text_pt\n"+"\n".join(rows[k] for k in sorted(rows))+"\n")
    if data is not None:
        if isinstance(data, dict):
            with open(f"{FIG}/{name}_data.tsv","w") as o:
                for k,v in data.items():
                    o.write(f"## {k}\n"); v.to_csv(o, sep="\t"); o.write("\n")
        else: data.to_csv(f"{FIG}/{name}_data.tsv", sep="\t")
    plt.close(fig); print(f"saved {FIG}/{name}.{{pdf,tiff}}  {w_mm:.0f} x {h_mm:.0f} mm, smallest text {fmin:.1f} pt")

def stack_labels(ax, xs, ys, names, x_text, color="#333333", fontsize=7, min_gap_pt=9.5, ha="left", style=None, ylim=None, box=True):
    """Phase O: place labels in a vertical stack at data x = x_text, each as close as possible to its point's y, at least min_gap_pt
    apart (1-D repel inside ylim), joined to the point by a leader line; white text background. Returns the text artists."""
    lo, hi = ylim if ylim else ax.get_ylim()
    pts_per_unit = ax.get_window_extent().height / ((hi - lo) if hi != lo else 1) * 72 / ax.figure.dpi
    gap = min_gap_pt / pts_per_unit
    order = np.argsort(ys); y = np.array(ys, float)[order]
    pos = y.copy()
    for k in range(1, len(pos)): pos[k] = max(pos[k], pos[k-1] + gap)
    over = pos[-1] - (hi - gap/2) if len(pos) else 0
    if over > 0: pos -= over
    for k in range(len(pos)-2, -1, -1): pos[k] = min(pos[k], pos[k+1] - gap)
    pos = np.clip(pos, lo + gap/2, hi - gap/2)
    out = []
    cols = color if isinstance(color, (list, tuple, np.ndarray)) else [color]*len(xs)
    for k, i in enumerate(order):
        out.append(ax.annotate(names[i], (xs[i], ys[i]), xytext=(x_text, pos[k]), textcoords="data", fontsize=fontsize, color=cols[i], ha=ha, va="center",
                               style=style, arrowprops=dict(arrowstyle="-", color=cols[i], lw=0.5, shrinkA=0, shrinkB=2),
                               bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.9) if box else None, zorder=9))
    return out

def label_points(ax, xs, ys, names, color="#333333", fontsize=7, dx=4, dy=2):
    for x,y,n in zip(xs,ys,names):
        ax.annotate(n, (x,y), xytext=(dx,dy), textcoords="offset points", fontsize=fontsize, color=color,
                    arrowprops=dict(arrowstyle="-", color=color, lw=0.4, shrinkB=2))

def clip_for_display(ax, y, lo=-3, hi=3, q=(0.001,0.999)):
    """clip extreme y values to [min(lo,q0.1%), max(hi,q99.9%)] for the scatter and note how many were clipped"""
    y=pd.Series(np.asarray(y,dtype=float)); a=min(lo,y.quantile(q[0])); b=max(hi,y.quantile(q[1])); n=int(((y<a)|(y>b)).sum())
    if n: ax.text(0.01,0.01,f"{n} gene(s) beyond the axis range shown at the limit",transform=ax.transAxes,fontsize=7,color=GRAY,va="bottom")
    ax.set_ylim(a-0.05*(b-a), b+0.05*(b-a)); return y.clip(a,b).values
