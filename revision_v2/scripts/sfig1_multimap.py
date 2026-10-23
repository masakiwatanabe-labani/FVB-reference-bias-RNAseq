#!/usr/bin/env python3
"""Supplementary Figure S2 (renumbered 2026-10-20; cited after the identity figure). Multi-mapping sensitivity (bias_multi from featureCounts -M --fraction vs unique-only bias; both in
analysis_v2/results/bias/all_tissues.<design>.tsv). (A) genome-wide Δ = bias_multi − bias_unique over gene × tissue pairs;
(B) catalog genes: mean unique vs mean multi bias, coloured by mechanism class (|Δ| > 0.3: SNP-induced multi-mapping);
(C) residual fraction (bias_multi / bias_unique) per catalog gene; (D) per-tissue unique vs multi bias for four catalog genes chosen
by rule: the two with the largest |Δ| and the two with residual fraction closest to 1."""
import matplotlib.patches as mpatches
from fig_common import *
long,tissues=load_bias(); cat=load_catalog(); C_MULTI="#E07B3A"; C_DIRECT="#4C72B0"
pairs=long.dropna(subset=["delta_multi_minus_unique"])
fig=plt.figure(figsize=(FULL_W,8.3)); gs=fig.add_gridspec(2,2,hspace=0.5,wspace=0.42,left=0.1,right=0.985,top=0.93,bottom=0.13)
_gc=gs[1,0].subgridspec(1,2,width_ratios=[1.55,1],wspace=0.05)   # Phase O-1d: the residual-fraction lists sit right of the bars, never on them
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(_gc[0]); axCt=fig.add_subplot(_gc[1]); axCt.axis("off"); axD=fig.add_subplot(gs[1,1])
for ax,l in zip((axA,axB,axC,axD),"ABCD"): panel(ax,l,x=-0.2)
d=pairs.delta_multi_minus_unique; out=(d.abs()>1).mean()
axA.hist(d.clip(-1,1),bins=np.arange(-1.025,1.05,0.05),color="#AAAAAA",alpha=0.8,edgecolor="none")
if len(cat):
    cd=pairs[pairs.gene_id.isin(cat.index)].delta_multi_minus_unique.clip(-1,1); axA.scatter(cd,np.full(len(cd),axA.get_ylim()[1]*0.45),c=C_MULTI,s=5,alpha=0.7,zorder=3,label=f"catalog gene\n× tissue pairs\n(n = {len(cd)})"); axA.legend(fontsize=7,loc="lower left",bbox_to_anchor=(0.0,0.1),handletextpad=0.3,borderpad=0.3)   # points at 78 % height so that the text box (top right) does not cover them (polish 2026-10-20)
axA.axvline(0,color="black",lw=1,ls="--"); axA.set_xlim(-1.05,1.05); despine(axA)
axA.text(0.98,0.98,f"n = {len(pairs):,}\ngene × tissue pairs\n({len(tissues)} tissues)\nmedian Δ = {d.median():+.3f}\nmean Δ = {d.mean():+.3f} log₂\n{out*100:.1f}% of pairs\nbeyond ±1 (clipped)",transform=axA.transAxes,ha="right",va="top",fontsize=7,bbox=dict(boxstyle="round,pad=0.3",facecolor="white",edgecolor=LGRAY))
axA.set_xlabel("Δ = bias(multi-mapping, fractional)\n− bias(unique-only) [log₂]",fontsize=7.5); axA.set_ylabel("Gene × tissue pairs",fontsize=8); axA.set_title("Genome-wide effect of including\nmulti-mapping reads",fontsize=8,pad=6)
if len(cat)==0:
    for ax in (axB,axC,axD): ax.text(0.5,0.5,"no catalog genes yet",ha="center",va="center",transform=ax.transAxes,color=GRAY); ax.axis("off")
    save(fig,"SFig2_multimapping",pd.DataFrame({"n_pairs":[len(pairs)],"median_delta":[d.median()],"mean_delta":[d.mean()]})); raise SystemExit
cat=cat.copy(); cat["cls"]=np.where(cat.delta_multi_minus_unique.abs()>DELTA_THR,"multi","direct"); col=cat.cls.map({"multi":C_MULTI,"direct":C_DIRECT})
# B
m=max(cat.mean_bias.max(),cat.mean_bias_multi.max())*1.1; axB.plot([0,m],[0,m],"--",color="grey",lw=1,label="Δ = 0"); axB.axhline(0,color="grey",lw=0.5,ls=":")
axB.scatter(cat.mean_bias,cat.mean_bias_multi,c=col,s=16,alpha=0.85,zorder=3); top=cat.reindex(cat.delta_multi_minus_unique.abs().sort_values(ascending=False).index).head(10)
axB.set_xlim(-0.7,m*1.42); axB.set_ylim(min(-1.6,cat.mean_bias_multi.min()-0.3),m)
stack_labels(axB,list(top.mean_bias),list(top.mean_bias_multi),list(top.gene_name),x_text=m*1.05,fontsize=7,min_gap_pt=9.5)   # Phase O-1c: stacked right of the points with leader lines
axB.legend(handles=[mpatches.Patch(color=C_MULTI,label=f"SNP-induced multi-mapping,\n|Δ| > {DELTA_THR} (n = {(cat.cls=='multi').sum()})"),mpatches.Patch(color=C_DIRECT,label=f"direct SNP misalignment,\n|Δ| ≤ {DELTA_THR} (n = {(cat.cls=='direct').sum()})")],fontsize=7,loc="upper left")
axB.set_xlabel("Mean bias, unique-only (log₂)",fontsize=8); axB.set_ylabel("Mean bias, multi-mapping included (log₂)",fontsize=7.5); axB.set_title("Catalog genes: unique vs\nmulti-mapping bias",fontsize=8,pad=6); despine(axB)
# C
cs=cat.sort_values("residual_fraction"); y=np.arange(len(cs)); axC.barh(y,cs.residual_fraction.clip(-0.5,1.7),color=cs.cls.map({"multi":C_MULTI,"direct":C_DIRECT}),height=0.7,alpha=0.85)
axC.axvline(0,color="black",lw=0.8); axC.axvline(1,color="black",lw=1,ls="--"); axC.set_xlim(-0.55,1.75)
if len(cs)<=40: axC.set_yticks(y); axC.set_yticklabels(cs.gene_name,fontsize=7.5,style="italic")
else:   # 138 names do not fit: no tick labels; the 5 lowest and 5 highest residual fractions are named, all values are in Supplementary Table S2 (polish 2026-10-20)
    axC.set_yticks([]); axC.set_ylabel("Catalog genes, sorted by residual fraction",fontsize=8)
    hi=cs.tail(5).iloc[::-1]; lo=cs.head(5)
    axCt.text(0.02,1.0,"highest residual fraction:\n"+"\n".join(f"{n} ({v:.2f})" for n,v in zip(hi.gene_name,hi.residual_fraction)),ha="left",va="top",fontsize=7,style="italic",color="#333333",transform=axCt.transAxes)
    axCt.text(0.02,0.0,"lowest residual fraction:\n"+"\n".join(f"{n} ({v:.2f})" for n,v in zip(lo.gene_name,lo.residual_fraction)),ha="left",va="bottom",fontsize=7,style="italic",color="#333333",transform=axCt.transAxes)
axC.set_xlabel("Residual fraction = bias(multi) / bias(unique)",fontsize=7.5); axC.set_title(f"Residual fraction per catalog gene\n(n = {len(cs)}; values clipped to [−0.5, 1.7])",fontsize=8,pad=6,loc="left"); despine(axC)
# D — four genes by rule
pick=list(cat.reindex(cat.delta_multi_minus_unique.abs().sort_values(ascending=False).index).head(2).index)
pick+=[g for g in (cat.residual_fraction-1).abs().sort_values().index if g not in pick][:2]
import matplotlib.transforms as mtransforms
# Phase O-1e: one slot of width w per tissue (unique and multi bars side by side, no overlap between tissues)
w=0.42; gap=2.3; rows=[]; xt=[]; xl=[]; blend=mtransforms.blended_transform_factory(axD.transData,axD.transAxes)
for gi,g in enumerate(pick):
    sub=long[long.gene_id==g].set_index("tissue"); c=C_MULTI if cat.loc[g,"cls"]=="multi" else C_DIRECT; gx=gi*gap
    for ti,t in enumerate(tissues):
        if t not in sub.index or np.isnan(sub.loc[t,"bias"]): continue
        xc=gx+(ti-(len(tissues)-1)/2)*w; axD.bar(xc-w*0.22,sub.loc[t,"bias"],width=w*0.42,color=c,alpha=0.85); axD.bar(xc+w*0.22,sub.loc[t,"bias_multi"],width=w*0.42,color=c,alpha=0.4,hatch="///",edgecolor=c)
        xt.append(xc); xl.append(t[0].upper()); rows.append({"gene_id":g,"gene_name":cat.loc[g,"gene_name"],"tissue":t,"bias_unique":sub.loc[t,"bias"],"bias_multi":sub.loc[t,"bias_multi"]})
    # gene name and class below the axis in axes coordinates (data-coordinate labels collided with the tissue ticks and negative bars; polish 2026-10-20)
    axD.text(gx,-0.13 if gi%2==0 else -0.19,cat.loc[g,"gene_name"],transform=blend,ha="center",va="top",fontsize=7,color=c,fontweight="bold",style="italic")   # class shown by colour (legend)
axD.axhline(0,color="black",lw=0.8); axD.set_xticks(xt); axD.set_xticklabels(xl,fontsize=7,color="#555555")   # Phase O-1e: single letters, explained in the legend; axD.tick_params(axis="x",length=2,pad=1)
axD.set_ylim(bottom=min(0,axD.get_ylim()[0])-0.2); axD.set_ylabel("Bias (log₂FC FVB ref − log₂FC naive)",fontsize=7.5); despine(axD)
axD.legend(handles=[mpatches.Patch(facecolor="grey",alpha=0.85,label="unique-only (primary)"),mpatches.Patch(facecolor="grey",alpha=0.4,hatch="///",edgecolor="grey",label="multi-mapping included\n(-M --fraction)"),mpatches.Patch(color=C_MULTI,label="SNP-induced multi-mapping"),mpatches.Patch(color=C_DIRECT,label="direct SNP misalignment")],fontsize=7,loc="upper right",title="A, adipose; H, hypothalamus;\nL, liver; S, skin",title_fontsize=7)
axD.set_title("Representative catalog genes per tissue\n(two largest |Δ| of the catalog;\ntwo with residual fraction closest to 1)",fontsize=8,pad=6)   # Phase N-1a: selection over the whole catalog
save(fig,"SFig2_multimapping",{"genome_wide":pd.DataFrame({"n_pairs":[len(pairs)],"median_delta":[d.median()],"mean_delta":[d.mean()],"frac_abs_gt_1":[out]}),"catalog":cat,"panelD":pd.DataFrame(rows).set_index(["gene_id","tissue"])})
