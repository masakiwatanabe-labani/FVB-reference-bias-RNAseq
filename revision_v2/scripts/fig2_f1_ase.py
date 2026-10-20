#!/usr/bin/env python3
"""Figure 2. F1 hybrid allele-specific expression (reused Phase A results, analysis_v2/results/f1_ase/).
(A) per-sample mean FVB allele fraction, naive vs WASP, both F1 datasets (f1_ase_per_sample.tsv; sites total >= 20, ref >= 3, alt >= 3);
(B) site-level naive vs WASP scatter for one representative GSE293440 sample — chosen by rule: the sample whose naive site count is the
    median of the dataset (ties: first by SRR), so no hand-picking; (C) per-site histogram, all GSE293440 samples pooled.
Significance stars in (A): one-sample t-test of the per-sample means against 0.50 (as in f1_ase_summary.tsv)."""
import glob, matplotlib.patches as mpatches
from scipy import stats
from fig_common import *
D=f"{FINAL}/f1_ase"; ps=pd.read_csv(f"{D}/f1_ase_per_sample.tsv", sep="\t"); summ=pd.read_csv(f"{D}/f1_ase_summary.tsv", sep="\t")
GSES=sorted(ps.GSE.unique(), key=lambda g: -ps[ps.GSE==g].SRR.nunique())   # larger dataset first
GCOL={GSES[0]:(BLUE,LBLUE)};
if len(GSES)>1: GCOL[GSES[1]]=(ORANGE,LORANGE)
man=pd.read_csv(f"{PROJ}/revision_v2/manifest/manifest.tsv", sep="\t"); man=man[man.role=="f1_ase"]
def gdesc(g):
    m=man[man.GSE==g]; return f"{g}\n{m.tissue.iloc[0].replace('_',' ')}, n = {m.SRR.nunique()}"
def load_filt(p):
    a=pd.read_csv(p, sep="\t"); tot=a.refCount+a.altCount; a=a[(tot>=20)&(a.refCount>=3)&(a.altCount>=3)].copy(); a["af"]=a.altCount/(a.refCount+a.altCount); return a[["contig","position","af"]]
def star(p): return "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "n.s."

fig=plt.figure(figsize=(13,5.5)); gs=fig.add_gridspec(1,3,width_ratios=[1.5,1.1,1.2],wspace=0.38,left=0.07,right=0.97,top=0.86,bottom=0.15)
axA=fig.add_subplot(gs[0]); axB=fig.add_subplot(gs[1]); axC=fig.add_subplot(gs[2])
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
    axA.text(xp, m+se+0.012, star(p), ha="center", va="bottom", fontsize=10)
    rows.append({"GSE":g,"arm":arm,"n":len(v),"mean":m,"se":se,"p_vs_0.5":p})
axA.axhline(0.5,color=RED,lw=1.2,ls="--",alpha=0.8); axA.text(x-0.25,0.503,"expected 0.50",fontsize=8.5,color=RED,ha="right",va="bottom")
for g in GSES:
    xc=(pos[(g,"naive")]+pos[(g,"wasp")])/2; axA.text(xc, -0.1, gdesc(g), ha="center", va="top", fontsize=8.5, color=GCOL[g][0], fontweight="bold", transform=axA.get_xaxis_transform())
axA.set_xticks([pos[k] for k in pos]); axA.set_xticklabels([ARM_SHORT[k[1]] for k in pos], fontsize=10); axA.set_xlim(-0.45,x-0.3)
axA.set_ylim(0.2,0.58); axA.set_ylabel("FVB allele fraction (per-sample mean over filtered sites)", fontsize=10); despine(axA)
axA.set_title("Per-sample FVB allele fraction in F1 hybrids", fontsize=11, pad=6); panel(axA,"A")
# B — representative sample by rule
g0=GSES[0]; nv=ps[(ps.GSE==g0)&(ps.arm=="naive")].sort_values(["n_sites","SRR"]); rep=nv.iloc[(len(nv)-1)//2].SRR
n_b=load_filt(f"{D}/{g0}/{rep}.naive.ase.tsv"); w_b=load_filt(f"{D}/{g0}/{rep}.wasp.ase.tsv"); mb=n_b.merge(w_b,on=["contig","position"],suffixes=("_naive","_wasp"))
axB.scatter(mb.af_naive, mb.af_wasp, s=12, alpha=0.4, color=GCOL[g0][0], edgecolors="none", zorder=3, rasterized=True)
axB.plot([0,1],[0,1],"k--",lw=0.8,alpha=0.5); axB.axhline(0.5,color=RED,lw=0.9,ls=":",alpha=0.7); axB.axvline(0.5,color=RED,lw=0.9,ls=":",alpha=0.7)
axB.text(0.05,0.96,f"n = {len(mb)} sites (both arms)\nnaive mean {mb.af_naive.mean():.3f}\nWASP mean {mb.af_wasp.mean():.3f}",transform=axB.transAxes,fontsize=9,va="top",bbox=dict(boxstyle="round,pad=0.3",facecolor="white",edgecolor=LGRAY))
axB.set_xlim(0,1); axB.set_ylim(0,1); axB.set_xlabel("FVB allele fraction, naive", fontsize=10); axB.set_ylabel("FVB allele fraction, WASP", fontsize=10); despine(axB)
axB.set_title(f"Site level, one {g0} sample\n({rep}; median site count of the dataset)", fontsize=10, pad=6); panel(axB,"B",x=-0.2)
# C — pooled
pool=[]
for p in sorted(glob.glob(f"{D}/{g0}/*.naive.ase.tsv")):
    srr=os.path.basename(p).split(".")[0]; m=load_filt(p).merge(load_filt(f"{D}/{g0}/{srr}.wasp.ase.tsv"),on=["contig","position"],suffixes=("_naive","_wasp")); m["SRR"]=srr; pool.append(m)
pool=pd.concat(pool); bins=np.linspace(0,1,41)
axC.hist(pool.af_naive,bins=bins,color=GCOL[g0][0],alpha=0.6,density=True,label=f"naive (mean {pool.af_naive.mean():.3f})",zorder=3)
axC.hist(pool.af_wasp,bins=bins,color=GCOL[g0][1],alpha=0.7,density=True,label=f"WASP (mean {pool.af_wasp.mean():.3f})",zorder=4)
axC.axvline(0.5,color=RED,lw=1.2,ls="--",zorder=5); axC.axvline(pool.af_naive.mean(),color=GCOL[g0][0],lw=1.2,ls=":",zorder=5); axC.axvline(pool.af_wasp.mean(),color=GCOL[g0][1],lw=1.2,ls=":",zorder=5)
axC.set_xlabel("FVB allele fraction per site", fontsize=10); axC.set_ylabel("Density", fontsize=10); despine(axC)
axC.legend(fontsize=8.5,loc="upper center",bbox_to_anchor=(0.5,-0.13),ncol=2,frameon=True,edgecolor=LGRAY)
axC.set_title(f"Per-site distribution, {g0}\n({pool.SRR.nunique()} samples pooled, {len(pool):,} site × sample)", fontsize=10, pad=6); panel(axC,"C",x=-0.2)
save(fig,"Fig2_f1_ase",{"panelA":pd.DataFrame(rows).set_index(["GSE","arm"]),"panelB_sample":pd.DataFrame({"representative_SRR":[rep],"n_sites":[len(mb)],"mean_naive":[mb.af_naive.mean()],"mean_wasp":[mb.af_wasp.mean()]}),
                         "panelC":pd.DataFrame({"n_site_sample":[len(pool)],"mean_naive":[pool.af_naive.mean()],"mean_wasp":[pool.af_wasp.mean()]})})
