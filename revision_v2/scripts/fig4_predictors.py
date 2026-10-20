#!/usr/bin/env python3
"""Figure 4. Gene-level predictors of bias (analysis_v2/results/predictors/gene_predictors.<design>.tsv: genes tested in >= 1 tissue, mean bias over the tissues tested (old Methods definition),
mean bias across tissues; SNP density / length from gene_meta.tsv; expression = mean logCPM in the naive arm).
(A) bias vs FVB/NJ SNP density with catalog genes as stars and the two worked-example genes (same rule as Fig. 6) labelled;
(B) bias vs gene length; (C) bias vs expression with binned mean ± SE (20 percentile bins); (D) bias by biotype group."""
from scipy import stats
from fig_common import *
from worked_examples import pick_examples
pr=load_predictors(); cat=load_catalog(); ex=pick_examples()
fig,axes=plt.subplots(2,2,figsize=(11,9)); fig.subplots_adjust(wspace=0.32,hspace=0.4,left=0.08,right=0.97,top=0.93,bottom=0.08)
for ax,l in zip(axes.flat,"ABCD"): panel(ax,l,x=-0.13)
BOX=dict(boxstyle="round,pad=0.3",facecolor="white",edgecolor=LGRAY,alpha=0.9)
# A
axA=axes[0,0]; a=pr[pr.snp_density_per_kb.notna()].copy(); cap=a.snp_density_per_kb.quantile(0.995); xs=a.snp_density_per_kb.clip(upper=cap)
ya=clip_for_display(axA,a.mean_bias); axA.scatter(xs,ya,s=3,alpha=0.12,color=BLUE,edgecolors="none",rasterized=True); axA.axhline(0,color="black",lw=0.8,ls="--",alpha=0.5)
q75=a.snp_density_per_kb.quantile(0.75); txt=""
if len(cat):
    c=a[a.index.isin(cat.index)]; axA.scatter(c.snp_density_per_kb.clip(upper=cap),c.mean_bias,s=80,color=GREEN,marker="*",edgecolors="white",lw=0.8,zorder=6,label=f"catalog genes (n = {len(c)})")
    ntop=int((c.snp_density_per_kb>=q75).sum()); txt=f"{ntop} of {len(c)} catalog genes in the top quartile\nof SNP density (≥ {q75:.2f} SNPs/kb)"
    axA.legend(fontsize=9,loc="lower right")
seen=set()
for key,col,dy in (("robust",GREEN,0.15),("extreme",RED,-0.25)):
    g=ex.get(key)
    if g and g in a.index and g not in seen:
        seen.add(g)
        axA.scatter(min(a.loc[g,"snp_density_per_kb"],cap),np.clip(a.loc[g,"mean_bias"],*axA.get_ylim()),s=120,color=col,marker="*",edgecolors="white",lw=1,zorder=8)
        yy=float(np.clip(a.loc[g,"mean_bias"],*axA.get_ylim())); axA.annotate(ex[key+"_name"],(min(a.loc[g,"snp_density_per_kb"],cap),yy),xytext=(min(a.loc[g,"snp_density_per_kb"],cap)+0.6,yy+dy),fontsize=9,color=col,arrowprops=dict(arrowstyle="-",color=col,lw=0.6))
rho,p=stats.spearmanr(a.snp_density_per_kb,a.mean_bias)
axA.text(0.97,0.42,(txt+"\n" if txt else "")+f"Spearman ρ = {rho:.3f}, p = {p:.1e}",transform=axA.transAxes,fontsize=9,va="top",ha="right",bbox=BOX)   # lower half of the panel (empty): the box at the top overlapped the catalog stars and the Rps2 label (polish 2026-10-20)
axA.set_xlabel("FVB/NJ SNP density (SNPs per kb of gene body)",fontsize=10.5); axA.set_ylabel("Mean bias across tissues (log₂)",fontsize=10.5); axA.set_title("Bias vs SNP density",fontsize=11.5,pad=6); despine(axA)
# B
axB=axes[0,1]; axB.scatter(pr.log10_length,clip_for_display(axB,pr.mean_bias),s=3,alpha=0.12,color=ORANGE,edgecolors="none",rasterized=True); axB.axhline(0,color="black",lw=0.8,ls="--",alpha=0.5)
rb,pb=stats.spearmanr(pr.log10_length,pr.mean_bias,nan_policy="omit"); axB.text(0.97,0.97,f"Spearman ρ = {rb:.3f}\np = {pb:.1e}",transform=axB.transAxes,fontsize=9,va="top",ha="right",bbox=BOX)
axB.set_xticks([3,4,5,6]); axB.set_xticklabels(["10³","10⁴","10⁵","10⁶"]); axB.set_xlabel("Gene length (bp)",fontsize=10.5); axB.set_ylabel("Mean bias across tissues (log₂)",fontsize=10.5); axB.set_title("Bias vs gene length",fontsize=11.5,pad=6); despine(axB)
# C
axC=axes[1,0]; c=pr.dropna(subset=["mean_logCPM_naive"]); axC.scatter(c.mean_logCPM_naive,clip_for_display(axC,c.mean_bias),s=3,alpha=0.12,color=GREEN,edgecolors="none",rasterized=True)
edges=np.unique(np.percentile(c.mean_logCPM_naive,np.linspace(0,100,21))); idx=np.digitize(c.mean_logCPM_naive,edges); bx=[];by=[];be=[]
for b in range(1,len(edges)):
    s=c.mean_bias[idx==b]
    if len(s)>10: bx.append((edges[b-1]+edges[b])/2); by.append(s.mean()); be.append(s.sem())
axC.errorbar(bx,by,yerr=be,fmt="o-",color=RED,ms=4,lw=1.4,capsize=3,zorder=5); axC.axhline(0,color="black",lw=0.8,ls="--",alpha=0.5)
rc,pc=stats.spearmanr(c.mean_logCPM_naive,c.mean_bias); axC.text(0.97,0.97,f"Spearman ρ = {rc:.3f}\np = {pc:.1e}",transform=axC.transAxes,fontsize=9,va="top",ha="right",bbox=BOX)
axC.set_xlabel("Expression (mean log₂CPM, naive arm)",fontsize=10.5); axC.set_ylabel("Mean bias across tissues (log₂)",fontsize=10.5); axC.set_title("Bias vs expression level (red: binned mean ± SE)",fontsize=11.5,pad=6); despine(axC)
# D
axD=axes[1,1]; order=[b for b in ("protein_coding","lncRNA","pseudogene","miRNA","other") if b in set(pr.biotype_group)]; cols={"protein_coding":BLUE,"lncRNA":GREEN,"pseudogene":ORANGE,"miRNA":PURPLE,"other":GRAY}
labels={"protein_coding":"Protein-coding","lncRNA":"lncRNA","pseudogene":"Pseudogene","miRNA":"miRNA","other":"Other"}
data=[pr[pr.biotype_group==b].mean_bias.dropna().values for b in order]
vp=axD.violinplot([np.clip(d,-3,2) for d in data],positions=range(len(order)),widths=0.6,showmedians=False,showextrema=False)
for pc_,b in zip(vp["bodies"],order): pc_.set_facecolor(cols[b]); pc_.set_alpha(0.7); pc_.set_edgecolor("white")
for i,d in enumerate(data): axD.scatter(i,np.median(d),s=25,color="white",edgecolors=cols[order[i]],linewidths=1.5,zorder=5)
axD.axhline(0,color="black",lw=0.8,ls="--",alpha=0.5); axD.set_xticks(range(len(order))); axD.set_xticklabels([f"{labels[b]}\n(n = {len(d):,})" for b,d in zip(order,data)],fontsize=9)
kw=stats.kruskal(*[d for d in data if len(d)>1]); axD.text(0.97,0.97,f"Kruskal–Wallis p = {kw.pvalue:.1e}",transform=axD.transAxes,fontsize=9,va="top",ha="right",bbox=BOX)
axD.set_ylim(-3,2); axD.set_ylabel("Mean bias across tissues (log₂)",fontsize=10.5); axD.set_title("Bias by gene biotype",fontsize=11.5,pad=6); despine(axD)
fig.suptitle(f"Gene-level predictors of SNP-mediated reference-mapping bias ({len(pr):,} genes tested in ≥ 1 tissue, mean bias over tested tissues; design: {DESIGN})",fontsize=11,y=0.985,color="#333333")
summ=pd.DataFrame([{"predictor":"snp_density_per_kb","rho":rho,"p":p,"n":len(a)},{"predictor":"log10_length","rho":rb,"p":pb,"n":int(pr.log10_length.notna().sum())},{"predictor":"mean_logCPM_naive","rho":rc,"p":pc,"n":len(c)}]).set_index("predictor")
bt=pd.DataFrame({"biotype_group":order,"n":[len(d) for d in data],"median_bias":[np.median(d) for d in data],"mean_bias":[d.mean() for d in data]}).set_index("biotype_group")
save(fig,"Fig4_predictors",{"spearman":summ,"biotype":bt,"catalog_top_quartile_snp_density":pd.DataFrame({"note":[txt.replace("\n"," ")]})})
