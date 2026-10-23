#!/usr/bin/env python3
"""Supplementary Figure S3 (Phase H-2a). Catalog genes: primary bias (x; mean over tested tissues) vs matched-reference bias (y;
mean over tested tissues; FVB/NJ animals on the FVB/NJ reference vs C57BL/6J animals on B6, minus the naive contrast). Genes
that meet the fixed catalog criterion under BOTH definitions are filled; the two worked-example genes are labelled. Diagonal
(bias = bias_matched, i.e. no C57BL/6J-read loss) and the 0.5 log2 criterion lines. Data: analysis_v2/results/catalog/catalog_decomposition.tsv."""
import matplotlib.patches as mpatches
from fig_common import *
d=pd.read_csv(f"{RES}/catalog/catalog_decomposition.tsv", sep="\t", index_col=0)
fig,ax=plt.subplots(figsize=(HALF_W,4.3)); fig.subplots_adjust(left=0.2,right=0.97,top=0.83,bottom=0.15)
both=d.in_matched_catalog=="Y"; m=max(d.mean_bias.max(), d.bias_matched_mean.max())*1.08; lo=min(0,d.bias_matched_mean.min())-0.3
ax.plot([lo,m],[lo,m],"--",color="grey",lw=1,label="bias = bias_matched\n(no C57BL/6J-read loss)"); ax.axhline(BIAS_THR,color=RED,lw=0.8,ls=":"); ax.axvline(BIAS_THR,color=RED,lw=0.8,ls=":")
ax.scatter(d.mean_bias[~both],d.bias_matched_mean[~both],s=12,facecolors="white",edgecolors=BLUE,lw=0.7,zorder=3,label=f"primary catalog only (n = {int((~both).sum())})")
ax.scatter(d.mean_bias[both],d.bias_matched_mean[both],s=14,color=BLUE,edgecolors="white",lw=0.3,zorder=4,label=f"supported under both definitions (n = {int(both.sum())})")
for g,col,dy,dx in (("ENSMUSG00000028648",GREEN,14,-34),("ENSMUSG00000044533",RED,-16,10)):
    if g in d.index:
        ax.scatter(d.loc[g,"mean_bias"],d.loc[g,"bias_matched_mean"],s=60,marker="*",color=col,edgecolors="white",lw=0.5,zorder=6)
        ax.annotate(d.loc[g,"gene_name"],(d.loc[g,"mean_bias"],d.loc[g,"bias_matched_mean"]),xytext=(dx,dy),textcoords="offset points",fontsize=7,style="italic",color=col,arrowprops=dict(arrowstyle="-",color=col,lw=0.6))
ax.set_xlim(lo,m); ax.set_ylim(lo,m); ax.set_xlabel("Primary bias: log₂FC(FVB ref) − log₂FC(naive),\nmean over tested tissues",fontsize=7)
ax.set_ylabel("Matched-reference bias: log₂FC(matched)\n− log₂FC(naive), mean over tested tissues",fontsize=7); despine(ax)
ax.legend(fontsize=7,loc="upper left",frameon=True,edgecolor=LGRAY,handletextpad=0.4,borderpad=0.3); ax.set_title(f"Decomposition of the catalog bias ({len(d)} genes);\nvertical distance to the diagonal = loss of\nC57BL/6J reads on the FVB/NJ reference",fontsize=7.5,pad=6)
save(fig,"SFig3_bias_decomposition",d)
