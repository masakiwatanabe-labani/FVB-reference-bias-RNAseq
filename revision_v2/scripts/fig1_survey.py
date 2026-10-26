#!/usr/bin/env python3
"""Figure 1. (A) reference-genome practice among GEO FVB RNA-seq series (analysis_v2/results/survey/survey_of_practice_rescreen.tsv:
the series that passed the Phase W-1 re-screen — RNA-Seq library strategy and FVB background per sample, survey/rescreen_classify.py;
before Phase W the panel used survey_of_practice_v2.tsv, 154 series); (B) schematic of SNP-mediated reference-mapping bias. The F1 allele fractions quoted in the
schematic are read from analysis_v2/results/f1_ase/v2/f1_ase_per_sample_v2.tsv (the per-animal table of Fig. 2; mean over the 12 F1
animals per arm), not typed in. Phase M-3 D1: the earlier path read the submitted-version table at the old site set."""
import matplotlib.patches as mpatches
from fig_common import *

df=pd.read_csv(f"{FINAL}/survey/survey_of_practice_rescreen.tsv", sep="\t")
n_total=len(df); n_b6=int((df.category=="stated_B6").sum()); n_nr=int((df.category=="not_reported").sum()); n_sa=int((df.category=="strain_aware").sum())
f1=pd.read_csv(f"{FINAL}/f1_ase/v2/f1_ase_per_sample_v2.tsv", sep="\t")
af_naive=f1[f1.arm=="naive"].mean_fvb_fraction.mean(); af_wasp=f1[f1.arm=="wasp"].mean_fvb_fraction.mean()

fig=plt.figure(figsize=(FULL_W,4.9)); gs=fig.add_gridspec(2,1,height_ratios=[1,2.75],hspace=0.22,left=0.2,right=0.97,top=0.94,bottom=0.01)
axA=fig.add_subplot(gs[0]); axB=fig.add_subplot(gs[1]); axA.set_position([0.21,0.755,0.76,0.2])
cats=["B6 reference reported","Reference not reported","Strain-aware correction"]; vals=[n_b6,n_nr,n_sa]; cols=[BLUE,GRAY,GREEN]
bars=axA.barh(range(3), vals, height=0.55, color=cols, zorder=3); axA.invert_yaxis()
for b,v,c in zip(bars,vals,cols):   # Phase W-1: label offset scaled to the axis (was +2 series at n = 154)
    axA.text(v+n_total*0.013, b.get_y()+b.get_height()/2, f"{v} ({v/n_total*100:.1f}%)", ha="left", va="center", fontsize=7.5, fontweight="bold", color=c if c!=GRAY else "#555555")
axA.set_yticks(range(3)); axA.set_yticklabels(cats, fontsize=7.5); axA.set_xlabel("Number of GEO series", fontsize=7.5)
axA.set_xlim(0,n_total*1.2); despine(axA)
axA.set_title(f"Reference-genome practice, {n_total} FVB RNA-seq series (GEO)", fontsize=8, pad=4); panel(axA,"A",x=-0.235)

axB.set_xlim(0,12); axB.set_ylim(0.2,10); axB.axis("off"); axB.set_position([0.02,0.0,0.96,0.585]); panel(axB,"B",x=0.0)
axB.set_title("SNP-mediated reference-mapping bias", fontsize=8, pad=8)
def box(x,y,w,h,color,label,fs=7.5,tc="white"):
    axB.add_patch(mpatches.FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.12",facecolor=color,edgecolor="none",zorder=3))
    axB.text(x+w/2,y+h/2,label,ha="center",va="center",fontsize=fs,color=tc,fontweight="bold",zorder=4)
def arrow(x1,y1,x2,y2,color="#555555",lw=1.5,ls="-"):
    axB.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="->",color=color,lw=lw,linestyle=ls))
box(0.3,8.5,5.2,1.0,"#4DAF4A","FVB/NJ (or F1 hybrid) RNA-seq reads",fs=7)
box(0.3,6.3,2.4,1.3,BLUE,"reads matching\nthe B6 allele",fs=7); box(3.1,6.3,2.4,1.3,ORANGE,"reads carrying\nthe FVB allele",fs=7)
arrow(1.5,8.5,1.5,7.7); arrow(4.3,8.5,4.3,7.7)
box(8.1,6.3,3.7,1.3,BLUE,"GRCm39\n(C57BL/6J reference)")
# Phase O-2d: the solid "maps" arrow starts at the blue B6-allele box and runs above the orange box; the dashed arrow stays on the orange box
axB.plot([2.35,2.35,6.8,6.8],[7.6,8.05,8.05,7.25],color=BLUE,lw=2,solid_capstyle="butt",zorder=2); arrow(6.8,7.3,8.0,7.3,color=BLUE,lw=2); arrow(5.6,6.6,8.0,6.6,color=ORANGE,lw=2,ls="dashed")
axB.text(6.15,8.12,"maps",fontsize=7,color=BLUE,ha="center",va="bottom"); axB.text(6.8,6.45,"mismatch penalty:\nunmapped, mis-placed\nor multi-mapped",fontsize=7,color=ORANGE,ha="center",va="top",style="italic")
box(7.7,4.1,1.9,1.2,BLUE,"B6-allele\nreads counted",fs=7); box(9.9,4.1,1.9,1.2,LGRAY,"FVB-allele\nreads lost",fs=7,tc="#555555")
arrow(9.3,6.3,8.7,5.4); arrow(10.6,6.3,10.8,5.4)
box(6.3,1.9,5.5,1.3,RED,"Apparent FVB under-expression\n(locus-specific artefact)"); arrow(9.8,4.1,9.8,3.3,color=RED,lw=2)
axB.text(0.3,5.0,f"F1 (B6 × FVB) hybrids: FVB allele fraction\nat heterozygous sites, expected 0.50\nnaive B6 mapping {af_naive:.3f}; WASP {af_wasp:.3f} (Fig. 2)",
         fontsize=7,color="#333333",va="top",bbox=dict(boxstyle="round,pad=0.3",facecolor="#FFF3CD",edgecolor="#DDAA00",linewidth=1))
axB.text(0.3,1.9,"Remedy tested here: align to an\nFVB/NJ-personalised reference (g2gtools)\nand compare with the naive and WASP arms",fontsize=7,color="#333333",va="top")
save(fig,"Fig1_survey",pd.DataFrame({"category":["stated_B6","not_reported","strain_aware"],"n":vals,"pct":[v/n_total*100 for v in vals]}).set_index("category"))
