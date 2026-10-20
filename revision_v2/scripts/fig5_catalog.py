#!/usr/bin/env python3
"""Figure 5. The bias catalog (analysis_v2/results/catalog/catalog.<design>.tsv; bias >= 0.5 log2 in >= 2 tissues, gene present in the
FVB/NJ GTF). (A) heatmap of per-tissue bias for the top 30 catalog genes (sorted by number of qualifying tissues, then mean bias);
(B) all catalog genes: tissues meeting the criterion vs mean bias, coloured by biotype, genes qualifying in every tissue labelled."""
import matplotlib.patches as mpatches
from fig_common import *
long,tissues=load_bias(); cat=load_catalog()
fig=plt.figure(figsize=(12,9)); gs=fig.add_gridspec(1,2,width_ratios=[1.6,1],wspace=0.4,left=0.16,right=0.97,top=0.9,bottom=0.1)
axA=fig.add_subplot(gs[0]); axB=fig.add_subplot(gs[1]); panel(axA,"A",x=-0.38,y=1.03); panel(axB,"B",x=-0.2,y=1.03)
if len(cat)==0:
    for ax in (axA,axB): ax.text(0.5,0.5,f"no catalog genes yet\n({len(tissues)} tissue(s) available; criterion needs ≥ {MIN_TISSUES})",ha="center",va="center",transform=ax.transAxes,color=GRAY); ax.axis("off")
    save(fig,"Fig5_catalog"); raise SystemExit
cat=cat.sort_values(["n_tissues","mean_bias"],ascending=False); top=cat.head(30); hm=top[tissues].values.astype(float)
vmax=max(np.nanpercentile(np.abs(hm),95),1.0); im=axA.imshow(hm,aspect="auto",cmap="RdBu_r",vmin=-vmax,vmax=vmax,interpolation="nearest")
axA.set_xticks(range(len(tissues))); axA.set_xticklabels([TISSUE_LABEL[t] for t in tissues],fontsize=10); axA.set_yticks(range(len(top))); axA.set_yticklabels(top.gene_name,fontsize=8.5,style="italic")
for i in range(len(top)):
    for j in range(len(tissues)):
        v=hm[i,j]
        if np.isnan(v): axA.text(j,i,"—",ha="center",va="center",fontsize=8,color="#AAAAAA")
        else: axA.text(j,i,f"{v:.1f}",ha="center",va="center",fontsize=8,color="white" if abs(v)>vmax*0.6 else "black",fontweight="bold" if v>=BIAS_THR else "normal")
cb=fig.colorbar(im,ax=axA,fraction=0.025,pad=0.03); cb.set_label("Bias (log₂FC FVB ref − log₂FC naive)",fontsize=9.5); cb.ax.tick_params(labelsize=8)
axA.set_title(f"Top {len(top)} of {len(cat)} catalog genes\n(— : not tested in that tissue; bold: ≥ {BIAS_THR})",fontsize=11,pad=8)
# B
def bgroup(b): return b if b in ("protein_coding","lncRNA") else ("pseudogene" if "pseudogene" in str(b) else "other")
col={"protein_coding":BLUE,"lncRNA":GREEN,"pseudogene":ORANGE,"other":GRAY}; lab={"protein_coding":"Protein-coding","lncRNA":"lncRNA","pseudogene":"Pseudogene","other":"Other"}
rng=np.random.default_rng(42); xj=cat.n_tissues.values+rng.uniform(-0.12,0.12,len(cat)); grp=cat.biotype.map(bgroup)
axB.scatter(xj,cat.mean_bias,s=np.clip(cat.mean_bias*18,20,200),c=[col[g] for g in grp],alpha=0.8,edgecolors="white",linewidths=0.5,zorder=3)
# labels: the 8 highest-mean-bias genes qualifying in every tissue plus any gene with mean bias >= 3 (all 46 four-tissue labels overlapped; polish 2026-10-20)
NLAB=8; labdf=pd.concat([cat[cat.n_tissues==len(tissues)].sort_values("mean_bias",ascending=False).head(NLAB),cat[cat.mean_bias>=3]]); labdf=labdf[~labdf.index.duplicated()]   # ("lab" is the biotype label dict below)
if len(labdf):
    pos={g:i for i,g in enumerate(cat.index)}; side=[7 if k%2==0 else -7 for k in range(len(labdf))]
    for (g,r),dx in zip(labdf.iterrows(),side): axB.annotate(r.gene_name,(xj[pos[g]],r.mean_bias),xytext=(dx,2),textcoords="offset points",fontsize=8,style="italic",ha="left" if dx>0 else "right",arrowprops=dict(arrowstyle="-",color="#333333",lw=0.4,shrinkB=2))
axB.axhline(BIAS_THR,color=RED,lw=0.8,ls="--",alpha=0.6); axB.set_xticks(range(MIN_TISSUES,len(tissues)+1)); axB.set_xlim(MIN_TISSUES-0.5,len(tissues)+0.5)
axB.set_xlabel(f"Tissues with bias ≥ {BIAS_THR} (of {len(tissues)})",fontsize=10.5); axB.set_ylabel("Mean bias across tested tissues (log₂)",fontsize=10.5); despine(axB)
axB.set_title(f"Catalog overview ({len(cat)} genes; point size ∝ mean bias;\ntop {NLAB} four-tissue genes and genes with mean bias ≥ 3 labelled)",fontsize=10.5,pad=8)
axB.legend(handles=[mpatches.Patch(facecolor=col[g],label=f"{lab[g]} (n = {int((grp==g).sum())})") for g in col if (grp==g).any()],fontsize=8.5,loc="upper center",bbox_to_anchor=(0.5,-0.1),ncol=2,frameon=True,edgecolor=LGRAY,title="Biotype",title_fontsize=9)
save(fig,"Fig5_catalog",cat)
