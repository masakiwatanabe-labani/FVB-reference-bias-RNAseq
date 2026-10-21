#!/usr/bin/env python3
"""Figure 2. F1 hybrid allele-specific expression, Phase H-1 recomputation (analysis_v2/results/f1_ase/v2/; f1_ase_rerun.sh +
f1_ase_summary_v2.py: GATK ASEReadCounter on the re-aligned F1 BAMs at the 203,510 FVB/NJ-homozygous exonic identity-QC sites,
genotype field set to 0/1 so that ASEReadCounter evaluates them).
(A) per-sample mean FVB allele fraction, naive vs WASP, both F1 datasets (f1_ase_per_sample_v2.tsv; sites total >= 20, ref >= 3, alt >= 3);
(B) site-level naive vs WASP scatter for one representative GSE293440 animal — chosen by rule: the animal whose naive site count is the
    median of the dataset (ties: first by SRR), so no hand-picking; (C) per-site histogram, all GSE293440 animals pooled;
(D) per-animal mean FVB allele fraction at the sites inside the 138 catalog genes vs the other evaluated sites, naive and WASP
    (f1_sites_catalog_vs_other.tsv; H-1d item 7).
Significance stars in (A): one-sample t-test of the per-sample means against 0.50 (as in f1_ase_summary_v2.tsv); in (D): paired t-test."""
import glob, matplotlib.patches as mpatches
from scipy import stats
from fig_common import *
D=f"{FINAL}/f1_ase/v2"; ps=pd.read_csv(f"{D}/f1_ase_per_sample_v2.tsv", sep="\t"); summ=pd.read_csv(f"{D}/f1_ase_summary_v2.tsv", sep="\t")
cat=pd.read_csv(f"{D}/f1_sites_catalog_vs_other.tsv", sep="\t").merge(ps[["sample_id","GSE","SRR"]].drop_duplicates(), on="sample_id")
GSES=sorted(ps.GSE.unique(), key=lambda g: -ps[ps.GSE==g].SRR.nunique())   # larger dataset first
GCOL={GSES[0]:(BLUE,LBLUE)};
if len(GSES)>1: GCOL[GSES[1]]=(ORANGE,LORANGE)
man=pd.read_csv(f"{PROJ}/revision_v2/manifest/manifest.tsv", sep="\t"); man=man[man.role=="f1_ase"]
def gdesc(g):
    m=man[man.GSE==g]; return f"{g}\n{m.tissue.iloc[0].replace('_',' ')}, n = {m.SRR.nunique()}"
def load_filt(p):
    a=pd.read_csv(p, sep="\t"); a=a[(a.totalCount>=20)&(a.refCount>=3)&(a.altCount>=3)].copy(); a["af"]=a.altCount/a.totalCount; return a[["contig","position","af"]]
def star(p): return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "n.s."

fig=plt.figure(figsize=(17,5.5)); gs=fig.add_gridspec(1,4,width_ratios=[1.5,1.1,1.2,1.25],wspace=0.4,left=0.05,right=0.985,top=0.86,bottom=0.15)
axA=fig.add_subplot(gs[0]); axB=fig.add_subplot(gs[1]); axC=fig.add_subplot(gs[2]); axD=fig.add_subplot(gs[3])
# A
pos={}; x=0.0; rows=[]
for g in GSES:
    for arm in ("naive","wasp"):
        pos[(g,arm)]=x; x+=0.7
    x+=0.3
rng=np.random.default_rng(42)
for (g,arm),xp in pos.items():
    v=ps[(ps.GSE==g)&(ps.arm==arm)].mean_fvb_fraction.values; c=GCOL[g][0 if arm=="naive" else 1]
    axA.scatter(xp+rng.uniform(-0.08,0.08,len(v)), v, s=55, color=c, zorder=5, edgecolors="white", linewidths=0.6)
    m=v.mean(); se=v.std(ddof=1)/np.sqrt(len(v)); axA.errorbar(xp,m,yerr=se,fmt="_",color="black",capsize=5,capthick=1.5,elinewidth=1.5,ms=12,zorder=6)
    p=summ[(summ.GSE==g)&(summ.arm==arm)]["p_vs_0.5"].iloc[0]
    axA.text(xp, m+se+0.004, star(p), ha="center", va="bottom", fontsize=10)
    axA.text(xp+0.12, m, f"{m:.3f}", ha="left", va="center", fontsize=8, color="#333333")
    rows.append({"GSE":g,"arm":arm,"n":len(v),"mean":m,"se":se,"p_vs_0.5":p})
axA.axhline(0.5,color=RED,lw=1.2,ls="--",alpha=0.8); axA.text(x-0.25,0.501,"expected 0.50",fontsize=8.5,color=RED,ha="right",va="bottom")
for g in GSES:
    xc=(pos[(g,"naive")]+pos[(g,"wasp")])/2; axA.text(xc, -0.1, gdesc(g), ha="center", va="top", fontsize=8.5, color=GCOL[g][0], fontweight="bold", transform=axA.get_xaxis_transform())
axA.set_xticks([pos[k] for k in pos]); axA.set_xticklabels([ARM_SHORT[k[1]] for k in pos], fontsize=10); axA.set_xlim(-0.45,x-0.3)
lo=min(ps.mean_fvb_fraction.min(),0.5)-0.02; hi=max(ps.mean_fvb_fraction.max(),0.5)+0.02; axA.set_ylim(lo,hi)
axA.set_ylabel("FVB allele fraction (per-animal mean over filtered sites)", fontsize=10); despine(axA)
axA.set_title("Per-animal FVB allele fraction in F1 hybrids", fontsize=11, pad=6); panel(axA,"A")
# B — representative animal by rule
g0=GSES[0]; nv=ps[(ps.GSE==g0)&(ps.arm=="naive")].sort_values(["n_sites","SRR"]); rep=nv.iloc[(len(nv)-1)//2]; rep_sid=rep.sample_id; rep=rep.SRR
n_b=load_filt(f"{D}/{rep_sid}.naive.ase.tsv"); w_b=load_filt(f"{D}/{rep_sid}.wasp.ase.tsv"); mb=n_b.merge(w_b,on=["contig","position"],suffixes=("_naive","_wasp"))
axB.scatter(mb.af_naive, mb.af_wasp, s=6, alpha=0.25, color=GCOL[g0][0], edgecolors="none", zorder=3, rasterized=True)
axB.plot([0,1],[0,1],"k--",lw=0.8,alpha=0.5); axB.axhline(0.5,color=RED,lw=0.9,ls=":",alpha=0.7); axB.axvline(0.5,color=RED,lw=0.9,ls=":",alpha=0.7)
axB.text(0.05,0.96,f"n = {len(mb):,} sites (both arms)\nnaive mean {mb.af_naive.mean():.3f}\nWASP mean {mb.af_wasp.mean():.3f}",transform=axB.transAxes,fontsize=9,va="top",bbox=dict(boxstyle="round,pad=0.3",facecolor="white",edgecolor=LGRAY))
axB.set_xlim(0,1); axB.set_ylim(0,1); axB.set_xlabel("FVB allele fraction, naive", fontsize=10); axB.set_ylabel("FVB allele fraction, WASP", fontsize=10); despine(axB)
axB.set_title(f"Site level, one {g0} animal\n({rep}; median site count of the dataset)", fontsize=10, pad=6); panel(axB,"B",x=-0.2)
# C — pooled
pool=[]
for _,r in ps[(ps.GSE==g0)&(ps.arm=="naive")].iterrows():
    m=load_filt(f"{D}/{r.sample_id}.naive.ase.tsv").merge(load_filt(f"{D}/{r.sample_id}.wasp.ase.tsv"),on=["contig","position"],suffixes=("_naive","_wasp")); m["SRR"]=r.SRR; pool.append(m)
pool=pd.concat(pool); bins=np.linspace(0,1,41)
axC.hist(pool.af_naive,bins=bins,color=GCOL[g0][0],alpha=0.6,density=True,label=f"naive (mean {pool.af_naive.mean():.3f})",zorder=3)
axC.hist(pool.af_wasp,bins=bins,color=GCOL[g0][1],alpha=0.7,density=True,label=f"WASP (mean {pool.af_wasp.mean():.3f})",zorder=4)
axC.axvline(0.5,color=RED,lw=1.2,ls="--",zorder=5); axC.axvline(pool.af_naive.mean(),color=GCOL[g0][0],lw=1.2,ls=":",zorder=5); axC.axvline(pool.af_wasp.mean(),color=GCOL[g0][1],lw=1.2,ls=":",zorder=5)
axC.set_xlabel("FVB allele fraction per site", fontsize=10); axC.set_ylabel("Density", fontsize=10); despine(axC)
axC.legend(fontsize=8.5,loc="upper center",bbox_to_anchor=(0.5,-0.13),ncol=2,frameon=True,edgecolor=LGRAY)
axC.set_title(f"Per-site distribution, {g0}\n({pool.SRR.nunique()} animals pooled, {len(pool):,} site × animal)", fontsize=10, pad=6); panel(axC,"C",x=-0.2)
# D — sites inside catalog genes vs the other sites, per animal (paired), naive and WASP
xpos={("naive","other"):0.0,("naive","catalog"):0.7,("wasp","other"):1.8,("wasp","catalog"):2.5}; drows=[]
for arm in ("naive","wasp"):
    a=cat[cat.arm==arm]
    for _,r in a.iterrows():
        c=GCOL[r.GSE][0 if arm=="naive" else 1]
        axD.plot([xpos[(arm,"other")],xpos[(arm,"catalog")]],[r.mean_af_other,r.mean_af_catalog],color=c,lw=0.8,alpha=0.6,zorder=3)
        axD.scatter([xpos[(arm,"other")],xpos[(arm,"catalog")]],[r.mean_af_other,r.mean_af_catalog],s=32,color=c,edgecolors="white",linewidths=0.5,zorder=5)
    t,p=stats.ttest_rel(a.mean_af_catalog,a.mean_af_other)
    for k,v in (("other",a.mean_af_other),("catalog",a.mean_af_catalog)):
        axD.errorbar(xpos[(arm,k)],v.mean(),yerr=v.sem(),fmt="_",color="black",capsize=4,capthick=1.3,elinewidth=1.3,ms=11,zorder=6)
        axD.text(xpos[(arm,k)]+(0.1 if k=="catalog" else -0.1),v.mean()-(0.004 if (arm,k)==("wasp","catalog") else 0),f"{v.mean():.3f}",ha="left" if k=="catalog" else "right",va="top" if (arm,k)==("wasp","catalog") else "center",fontsize=8,color="#333333")
    if arm=="naive":   # bracket below the points (the 0.50 line and its label sit just above them)
        yb=min(a.mean_af_other.min(),a.mean_af_catalog.min())-0.012; axD.plot([xpos[(arm,"other")],xpos[(arm,"catalog")]],[yb,yb],color="black",lw=0.8); axD.text((xpos[(arm,"other")]+xpos[(arm,"catalog")])/2,yb-0.003,f"paired t-test {star(p)}",ha="center",va="top",fontsize=8.5)
    else:
        ytop=max(a.mean_af_other.max(),a.mean_af_catalog.max())+0.01; axD.plot([xpos[(arm,"other")],xpos[(arm,"catalog")]],[ytop,ytop],color="black",lw=0.8); axD.text((xpos[(arm,"other")]+xpos[(arm,"catalog")])/2,ytop+0.003,f"paired t-test {star(p)}",ha="center",va="bottom",fontsize=8.5)
    drows.append({"arm":arm,"n_animals":len(a),"mean_af_other":a.mean_af_other.mean(),"mean_af_catalog":a.mean_af_catalog.mean(),"paired_t":t,"paired_p":p,"catalog_sites_min":int(a.n_sites_catalog.min()),"catalog_sites_max":int(a.n_sites_catalog.max())})
axD.axhline(0.5,color=RED,lw=1.2,ls="--",alpha=0.8); axD.text(-0.45,0.5015,"expected 0.50",fontsize=8,color=RED,ha="left",va="bottom")
axD.set_xticks(list(xpos.values())); axD.set_xticklabels(["other\nsites","catalog\ngenes","other\nsites","catalog\ngenes"],fontsize=9); axD.set_xlim(-0.5,3.0)
for arm,xc in (("naive",0.35),("wasp",2.15)): axD.text(xc,-0.16,ARM_SHORT[arm],ha="center",va="top",fontsize=10,fontweight="bold",transform=axD.get_xaxis_transform())
axD.set_ylabel("FVB allele fraction (per-animal mean)", fontsize=10); despine(axD)
axD.set_title("Sites inside the 138 catalog genes\nvs the other evaluated sites", fontsize=10, pad=6); panel(axD,"D",x=-0.2)
save(fig,"Fig2_f1_ase",{"panelA":pd.DataFrame(rows).set_index(["GSE","arm"]),"panelB_sample":pd.DataFrame({"representative_SRR":[rep],"sample_id":[rep_sid],"n_sites":[len(mb)],"mean_naive":[mb.af_naive.mean()],"mean_wasp":[mb.af_wasp.mean()]}),
                         "panelC":pd.DataFrame({"n_site_sample":[len(pool)],"mean_naive":[pool.af_naive.mean()],"mean_wasp":[pool.af_wasp.mean()]}),"panelD":pd.DataFrame(drows).set_index("arm")})
