#!/usr/bin/env python3
"""Phase H-1d, optional item 7 (report only): in the F1 animals, is the naive-arm FVB allele fraction lower at sites that lie inside
catalog genes (analysis_v2/results/catalog/catalog.primary.tsv; gene bodies from gene_meta.tsv) than at the other evaluated sites?
Per animal and pooled, naive and WASP, filtered sites (total >= 20, ref >= 3, alt >= 3). Also sites inside any tested-but-non-catalog gene
body vs intergenic/other. Output: analysis_v2/results/f1_ase/v2/f1_sites_catalog_vs_other.tsv"""
import pandas as pd, numpy as np, glob, os
from scipy import stats
P="/usr/local/jupyter/FVB_bias"; D=f"{P}/analysis_v2/results/f1_ase/v2"
gm=pd.read_csv(f"{P}/analysis_v2/results/gene_meta.tsv",sep="\t",dtype={"chrom":str}); cat=set(pd.read_csv(f"{P}/analysis_v2/results/catalog/catalog.primary.tsv",sep="\t").gene_id)
def intervals(sub):
    out={}
    for c,g in sub.groupby("chrom"):
        iv=sorted(zip(g.start,g.end)); merged=[]
        for s,e in iv:
            if merged and s<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],e)
            else: merged.append([s,e])
        out[c]=(np.array([m[0] for m in merged]),np.array([m[1] for m in merged]))
    return out
def inside(df,iv):
    r=np.zeros(len(df),bool)
    for c,g in df.groupby("contig"):
        if c not in iv: continue
        s,e=iv[c]; i=np.searchsorted(s,g.position.values,side="right")-1; ok=(i>=0)&(g.position.values<=e[np.clip(i,0,len(e)-1)]); r[g.index]=ok
    return r
IV_CAT=intervals(gm[gm.gene_id.isin(cat)]); IV_GENE=intervals(gm)
rows=[]
for p in sorted(glob.glob(f"{D}/*.naive.ase.tsv")):
    sid=os.path.basename(p).split(".naive")[0]
    for arm in ("naive","wasp"):
        a=pd.read_csv(f"{D}/{sid}.{arm}.ase.tsv",sep="\t",dtype={"contig":str}); a=a[(a.totalCount>=20)&(a.refCount>=3)&(a.altCount>=3)].reset_index(drop=True); a["af"]=a.altCount/a.totalCount
        a["in_catalog_gene"]=inside(a,IV_CAT); a["in_any_gene"]=inside(a,IV_GENE)
        c=a[a.in_catalog_gene]; o=a[~a.in_catalog_gene]; og=a[(~a.in_catalog_gene)&a.in_any_gene]
        t,pv=stats.mannwhitneyu(c.af,o.af,alternative="two-sided") if len(c)>1 else (np.nan,np.nan)
        rows.append({"sample_id":sid,"arm":arm,"n_sites_catalog":len(c),"n_sites_other":len(o),"n_sites_other_genic":len(og),"mean_af_catalog":c.af.mean(),"mean_af_other":o.af.mean(),"mean_af_other_genic":og.af.mean(),
                     "median_af_catalog":c.af.median(),"median_af_other":o.af.median(),"frac_catalog_sites_af_lt_0.5":(c.af<0.5).mean(),"frac_other_sites_af_lt_0.5":(o.af<0.5).mean(),"mannwhitney_p":pv})
df=pd.DataFrame(rows); df.to_csv(f"{D}/f1_sites_catalog_vs_other.tsv",sep="\t",index=False)
pd.set_option("display.width",250); print(df.round(4).to_string(index=False))
for arm in ("naive","wasp"):
    a=df[df.arm==arm]; t,pv=stats.ttest_rel(a.mean_af_catalog,a.mean_af_other)
    print(f"{arm}: per-animal mean AF catalog {a.mean_af_catalog.mean():.4f} vs other {a.mean_af_other.mean():.4f} (other genic {a.mean_af_other_genic.mean():.4f}); paired t = {t:.2f}, p = {pv:.2g}; catalog sites/animal {a.n_sites_catalog.min()}-{a.n_sites_catalog.max()}")
