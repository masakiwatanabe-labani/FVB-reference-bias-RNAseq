#!/usr/bin/env python3
"""Figure 3. Pure-strain FVB vs B6 comparison, three alignment arms (analysis_v2/results/bias/all_tissues.<design>.tsv).
(A) per-gene log2FC naive vs FVB-reference for the representative tissue — chosen by rule: the tissue with the most animals
    (ties: TISSUE_ORDER) — catalog genes coloured by the number of tissues meeting the criterion;
(B) per-gene bias (log2FC_fvb − log2FC_naive) per tissue, violins with median; (C) tissue-level mean log2FC ± SE per arm."""
import matplotlib.patches as mpatches
from fig_common import *
long,tissues=load_bias(); cat=load_catalog(); meta=gene_meta()
rep=max(tissues, key=lambda t: (sum(n_animals(t)), -TISSUE_ORDER.index(t)))
fig=plt.figure(figsize=(FULL_W,6.6)); gs=fig.add_gridspec(2,2,wspace=0.42,hspace=0.55,left=0.1,right=0.985,top=0.9,bottom=0.1)
axA=fig.add_subplot(gs[0,0]); axL=fig.add_subplot(gs[0,1]); axL.axis("off"); axB=fig.add_subplot(gs[1,0]); axC=fig.add_subplot(gs[1,1])
# A
r=long[long.tissue==rep].dropna(subset=["logFC_naive","logFC_fvb"]); LIM=4
rr=r[(r.logFC_naive.abs()<=LIM)&(r.logFC_fvb.abs()<=LIM)]
axA.scatter(rr.logFC_naive, rr.logFC_fvb, s=2, alpha=0.35, color=LGRAY, edgecolors="none", rasterized=True, zorder=1)
ncol={2:ORANGE,3:"#D04030",4:RED}; nsize={2:12,3:16,4:20}; leg=[]
if len(cat):
    cr=r[r.gene_id.isin(cat.index)].copy(); cr["n_tissues"]=cat.loc[cr.gene_id,"n_tissues"].values; cr["gene_name"]=cat.loc[cr.gene_id,"gene_name"].values
    for k in sorted(cr.n_tissues.unique()):
        s=cr[cr.n_tissues==k]; axA.scatter(s.logFC_naive.clip(-LIM,LIM), s.logFC_fvb.clip(-LIM,LIM), s=nsize.get(k,42), color=ncol.get(k,RED), edgecolors="white", lw=0.5, zorder=4)
        leg.append(mpatches.Patch(facecolor=ncol.get(k,RED), label=f"catalog, bias ≥ {BIAS_THR} in {k} of {len(tissues)} tissues\n(n = {len(s)} tested in {TISSUE_LABEL[rep].lower()})"))   # N-2b
    top=cat.sort_values(["n_tissues","mean_bias"],ascending=False).head(7); lab=cr[cr.gene_id.isin(top.index)]
    label_points(axA, lab.logFC_naive.clip(-LIM,LIM), lab.logFC_fvb.clip(-LIM,LIM), lab.gene_name, color="#990000", fontsize=7.5, dx=6, dy=6)
leg.append(mpatches.Patch(facecolor=LGRAY, label="all tested genes"))
axA.plot([-LIM,LIM],[-LIM,LIM],"k--",lw=0.8,alpha=0.4); axA.axhline(0,color=GREEN,lw=0.7,ls=":",alpha=0.5); axA.axvline(0,color=BLUE,lw=0.7,ls=":",alpha=0.5)
axA.set_xlim(-LIM,LIM); axA.set_ylim(-LIM,LIM); axA.set_xlabel("log₂FC FVB/B6, naive (B6 reference)", fontsize=8); axA.set_ylabel("log₂FC FVB/B6, FVB/NJ reference", fontsize=8); despine(axA)
nf,nb=n_animals(rep)
axA.set_title(f"{TISSUE_LABEL[rep]} ({nf} FVB/NJ vs {nb} C57BL/6J;\n{len(r):,} genes tested; design: {DESIGN})", fontsize=8, pad=6); panel(axA,"A",x=-0.2)
axL.legend(handles=leg, fontsize=7, loc="upper left", frameon=True, edgecolor=LGRAY, title="Panel A")
axL.text(0.0,0.12,f"catalog: {len(cat)} genes\n(bias ≥ {BIAS_THR} log₂ in ≥ {MIN_TISSUES} of {len(tissues)} tissues);\ncatalog genes tested in {TISSUE_LABEL[rep].lower()}: {len(cr) if len(cat) else 0} of {len(cat)}",transform=axL.transAxes,fontsize=7,va="bottom",ha="left",bbox=dict(boxstyle="round,pad=0.3",facecolor="white",edgecolor=LGRAY))
# B
vals=[long[long.tissue==t].bias.dropna().values for t in tissues]; vp=axB.violinplot([np.clip(v,-3,3) for v in vals],positions=range(len(tissues)),widths=0.7,showmedians=False,showextrema=False)
for pc,t in zip(vp["bodies"],tissues): pc.set_facecolor(TISSUE_COLOR[t]); pc.set_alpha(0.75); pc.set_edgecolor("white")
statsB=[]
for i,(t,v) in enumerate(zip(tissues,vals)):
    axB.scatter(i,np.median(v),s=14,color="white",edgecolors=TISSUE_COLOR[t],linewidths=1.3,zorder=5); axB.scatter(i,np.percentile(v,95),s=10,color=TISSUE_COLOR[t],marker="^",zorder=4)
    statsB.append({"tissue":t,"n_genes":len(v),"median":np.median(v),"mean":v.mean(),"q95":np.percentile(v,95),"frac_ge_0.5":(v>=0.5).mean(),"frac_le_-0.5":(v<=-0.5).mean()})
axB.axhline(0,color="black",lw=0.8,ls="--",alpha=0.5); axB.axhline(BIAS_THR,color=RED,lw=0.7,ls=":",alpha=0.6); axB.text(len(tissues)-0.5,BIAS_THR+0.03,"catalog\nthreshold",fontsize=7.5,color=RED,ha="right",va="bottom")
axB.set_xticks(range(len(tissues))); axB.set_xticklabels([TISSUE_LABEL[t] for t in tissues],fontsize=7.5,rotation=28,ha="right"); axB.set_ylim(-1.5,2); despine(axB)
axB.set_ylabel("Bias (log₂FC FVB ref − log₂FC naive)", fontsize=8); axB.set_title("Per-gene bias per tissue\n(circle: median; triangle: 95th percentile)", fontsize=8, pad=6); panel(axB,"B",x=-0.2)
# C
rows=[]
for t in tissues:
    d=long[long.tissue==t]
    for arm in ARM_ORDER:
        v=d[f"logFC_{arm}"].dropna(); rows.append({"tissue":t,"arm":arm,"mean":v.mean(),"se":v.std(ddof=1)/np.sqrt(len(v)),"n":len(v)})
sC=pd.DataFrame(rows); x=np.arange(len(tissues)); off={"naive":-0.22,"wasp":0,"fvb":0.22}; fmt={"naive":"o-","wasp":"s--","fvb":"^-"}
for arm in ARM_ORDER:
    s=sC[sC.arm==arm].set_index("tissue").loc[tissues]; axC.errorbar(x+off[arm],s["mean"],yerr=s["se"],fmt=fmt[arm],color=ARM_COLOR[arm],ms=4,lw=1.0,capsize=2,label=ARM_SHORT[arm],alpha=0.85 if arm!="wasp" else 0.7)
axC.axhline(0,color="black",lw=0.8,ls="--",alpha=0.5); axC.set_xticks(x); axC.set_xticklabels([TISSUE_LABEL[t] for t in tissues],fontsize=7.5,rotation=28,ha="right")
axC.set_ylabel("Mean log₂FC FVB/B6 over genes ± SE", fontsize=8); despine(axC); axC.legend(fontsize=7.5,loc="best",frameon=True,edgecolor=LGRAY)
axC.set_title("Genome-wide mean per arm\n(no global shift between arms)", fontsize=8, pad=6); panel(axC,"C",x=-0.22)
save(fig,"Fig3_strain_comparison",{"panelA_representative":pd.DataFrame({"tissue":[rep],"n_FVB":[nf],"n_B6":[nb],"n_genes":[len(r)]}),"panelB":pd.DataFrame(statsB).set_index("tissue"),"panelC":sC.set_index(["tissue","arm"])})
