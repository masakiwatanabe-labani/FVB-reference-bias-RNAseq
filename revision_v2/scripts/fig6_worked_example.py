#!/usr/bin/env python3
"""Figure 6. Worked examples chosen by rule (worked_examples.py): "robust" = catalog gene qualifying in the most tissues (ties: mean
bias); "extreme" = catalog gene of the other mechanism class (multi-mapping vs direct misalignment) with the largest single-tissue bias. (A) robust gene: log2FC per arm per tissue; (B) extreme gene:
per-sample log2(CPM + 1) in its highest-bias tissue, three arms (FVB-arm counts divided by the naive-arm library size of the same
sample, shared denominator); (C) SNP density vs mean bias, both genes highlighted; (D) extreme gene: log2FC per arm per tissue."""
from matplotlib.lines import Line2D
from fig_common import *
from worked_examples import pick_examples
long,tissues=load_bias(); pr=load_predictors(); meta=gene_meta(); ex=pick_examples()
G1,G2,T2=ex["robust"],ex["extreme"],ex["extreme_tissue"]; N1,N2=ex["robust_name"],ex["extreme_name"]
b1=long[long.gene_id==G1].set_index("tissue"); b2=long[long.gene_id==G2].set_index("tissue")
fig=plt.figure(figsize=(13,9)); gs=fig.add_gridspec(2,2,wspace=0.36,hspace=0.5,left=0.08,right=0.97,top=0.9,bottom=0.09)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])
for ax,l in zip((axA,axB,axC,axD),"ABCD"): panel(ax,l,x=-0.13,y=1.05)
def arm_bars(ax,b,name):
    x=np.arange(len(tissues)); w=0.25
    for i,arm in enumerate(ARM_ORDER):
        v=[b.loc[t,f"logFC_{arm}"] if t in b.index else np.nan for t in tissues]; ax.bar(x+(i-1)*w,v,w,color=ARM_COLOR[arm],alpha=0.85,label=ARM_SHORT[arm],edgecolor="white",zorder=3)
    ax.axhline(0,color="black",lw=0.8,ls="--",alpha=0.5); ax.set_xticks(x); ax.set_xticklabels([TISSUE_LABEL[t] for t in tissues],fontsize=9.5); ax.set_ylabel("log₂FC FVB/B6",fontsize=10.5); despine(ax)
    ax.legend(fontsize=9,loc="upper center",bbox_to_anchor=(0.5,-0.12),ncol=3,frameon=True,edgecolor=LGRAY)
    for t in tissues:
        if t in b.index and not np.isnan(b.loc[t,"bias"]): ax.text(tissues.index(t),max(b.loc[t,["logFC_naive","logFC_wasp","logFC_fvb"]].max(),0)+0.05,f"bias {b.loc[t,'bias']:+.2f}",ha="center",va="bottom",fontsize=8,color="#333333")
arm_bars(axA,b1,N1); arm_bars(axD,b2,N2)
axA.set_title(f"{N1} — log₂FC per arm; qualifies in {int(b1.bias.ge(BIAS_THR).sum())} of {len(tissues)} tissues\n(mean bias {b1.bias.mean():+.2f} log₂)",fontsize=10.5,pad=6,style="italic")
axD.set_title(f"{N2} — log₂FC per arm; largest single-tissue bias\n({TISSUE_LABEL[T2].replace(chr(10),' ')} {b2.loc[T2,'bias']:+.2f} log₂)",fontsize=10.5,pad=6,style="italic")
# B per-sample CPM
ss=load_samples(T2).set_index("sample_id"); nv=load_counts(T2,"naive"); lib=nv.sum(); cp={}
for arm in ARM_ORDER:
    c=load_counts(T2,arm); cp[arm]=(c.loc[G2]/lib[c.columns]*1e6) if G2 in c.index else pd.Series(dtype=float)   # all arms over the naive library size of the same sample
xp={"naive":0,"wasp":3,"fvb":6}; rows=[]
for arm in ARM_ORDER:
    v=np.log2(cp[arm]+1)
    for k,(strain,dx,fill) in enumerate((("FVB",0,True),("B6",1,False))):
        ids=[s for s in v.index if ss.loc[s,"strain_short"]==strain]; y=v[ids].values; jx=xp[arm]+dx+np.random.default_rng(k+1).uniform(-0.22,0.22,len(y))
        axB.scatter(jx,y,s=40,color=ARM_COLOR[arm] if fill else "white",edgecolors=ARM_COLOR[arm],lw=1.2,zorder=4); axB.plot([xp[arm]+dx-0.35,xp[arm]+dx+0.35],[y.mean()]*2,color=ARM_COLOR[arm],lw=2,ls="-" if fill else "--",zorder=5)
        for s,yy in zip(ids,y): rows.append({"arm":arm,"strain":strain,"sample_id":s,"log2_cpm_plus1":yy})
axB.set_xticks([0.5,3.5,6.5]); axB.set_xticklabels([ARM_SHORT[a] for a in ARM_ORDER],fontsize=10)
for tk,a in zip(axB.get_xticklabels(),ARM_ORDER): tk.set_color(ARM_COLOR[a])
nf,nb=n_animals(T2)
axB.legend(handles=[Line2D([0],[0],marker="o",color="w",markerfacecolor=GRAY,markersize=8,label=f"FVB/NJ (filled, n = {nf})"),Line2D([0],[0],marker="o",color="w",markerfacecolor="w",markeredgecolor=GRAY,markersize=8,label=f"C57BL/6J (open, n = {nb})")],fontsize=9,loc="upper center",bbox_to_anchor=(0.5,-0.12),ncol=2,frameon=True,edgecolor=LGRAY)
axB.set_xlim(-0.8,8.3); axB.set_ylabel("log₂(CPM + 1), naive-arm library size",fontsize=10.5); despine(axB)
axB.set_title(f"{N2} — per-animal expression, {TISSUE_LABEL[T2].replace(chr(10),' ')}\n(bar: group mean; bias {b2.loc[T2,'bias']:+.2f} log₂)",fontsize=10.5,pad=6,style="italic")
# C context
a=pr[pr.snp_density_per_kb.notna()]; cap=a.snp_density_per_kb.quantile(0.99); xs=a.snp_density_per_kb.clip(upper=cap)
yc=clip_for_display(axC,a.mean_bias); axC.scatter(xs,yc,s=3,alpha=0.1,color=GRAY,edgecolors="none",rasterized=True,label="all tested genes")
edges=np.unique(np.percentile(xs,np.linspace(0,100,21))); idx=np.digitize(xs,edges); bx=[];by=[]
for b in range(1,len(edges)):
    s=a.mean_bias[idx==b]
    if len(s)>20: bx.append((edges[b-1]+edges[b])/2); by.append(s.mean())
axC.plot(bx,by,"o-",color="#666666",ms=4,lw=1.4,alpha=0.7,label="binned mean")
for g,name,col,mk,dxy in (((G1,N1,GREEN,"*",(1.5,0.5)),) if G1==G2 else ((G1,N1,GREEN,"*",(1.5,0.5)),(G2,N2,RED,"o",(1.5,-0.6)))):
    if g in a.index:
        x0=min(a.loc[g,"snp_density_per_kb"],cap); y0=float(np.clip(a.loc[g,"mean_bias"],*axC.get_ylim())); lo,hi=axC.get_ylim(); dxy=(dxy[0], -abs(dxy[1])-0.3) if y0>lo+0.8*(hi-lo) else dxy; axC.scatter(x0,y0,s=170 if mk=="*" else 90,color=col,marker=mk,edgecolors="white",lw=1.2,zorder=7,label=f"{name} ({a.loc[g,'snp_density_per_kb']:.1f} SNPs/kb)")
        axC.annotate(name,(x0,y0),xytext=(x0+dxy[0],y0+dxy[1]),fontsize=10,color=col,fontweight="bold",style="italic",arrowprops=dict(arrowstyle="-",color=col,lw=0.8),bbox=dict(boxstyle="round,pad=0.2",facecolor="white",edgecolor=col,alpha=0.85,lw=0.8))
axC.axhline(0,color="black",lw=0.8,ls="--",alpha=0.5); axC.axhline(BIAS_THR,color=RED,lw=0.7,ls=":",alpha=0.5)
axC.set_xlabel("FVB/NJ SNP density (SNPs/kb)",fontsize=10.5); axC.set_ylabel("Mean bias across tissues (log₂)",fontsize=10.5); despine(axC); axC.legend(fontsize=9,loc="lower right",frameon=True,edgecolor=LGRAY)
axC.set_title("Genome-wide context: SNP density vs bias",fontsize=10.5,pad=6)
fig.suptitle(("TEST ONLY (no catalog yet) — " if ex["test_only"] else "")+f"Worked examples of SNP-mediated reference-mapping bias (design: {DESIGN})",fontsize=11,y=0.985,color="#333333")
save(fig,"Fig6_worked_examples",{"choice":pd.DataFrame([ex]),"robust_gene":b1,"extreme_gene":b2,"panelB":pd.DataFrame(rows).set_index(["arm","strain","sample_id"])})
