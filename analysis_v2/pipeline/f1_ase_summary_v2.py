#!/usr/bin/env python3
"""Phase H-1b: per-sample F1 allele-specific expression from the recomputed per-site tables (f1_ase_rerun.sh; GATK ASEReadCounter at
the 203,510 identity-QC sites), with the OLD Methods per-site filter: totalCount >= 20 and refCount >= 3 and altCount >= 3.
FVB allele fraction per site = altCount/totalCount (alt = FVB/NJ allele; B6 reference coordinates). Per-sample mean, per-dataset
means, pooled one-sample t-test vs 0.50 (12 animals), naive and WASP. The same statistics are also computed on the OLD site set
(94,351 MGP-heterozygous positions; f1_oldsites_recount.sh) when those tables exist, to test whether the v1 values (0.333/0.386)
are reproduced by site selection alone.
Outputs: analysis_v2/results/f1_ase/v2/f1_ase_per_sample_v2.tsv, f1_ase_summary_v2.tsv (+ *_oldsites.tsv)."""
import os, glob, numpy as np, pandas as pd
from scipy import stats
PROJ="/usr/local/jupyter/FVB_bias"; OUT=f"{PROJ}/analysis_v2/results/f1_ase/v2"
man=pd.read_csv(f"{PROJ}/revision_v2/manifest/manifest.tsv", sep="\t", dtype=str).set_index("sample_id"); f1=man[(man.role=="f1_ase")&(man["include"]=="Y")]
def summarise(suffix, tag):
    rows=[]
    for sid in f1.index:
        for arm in ("naive","wasp"):
            p=f"{OUT}/{sid}.{arm}{suffix}.ase.tsv"
            if not os.path.exists(p): continue
            a=pd.read_csv(p, sep="\t"); n_raw=len(a); a=a[(a.totalCount>=20)&(a.refCount>=3)&(a.altCount>=3)]; af=a.altCount/a.totalCount
            t,pv=stats.ttest_1samp(af,0.5) if len(af)>1 else (np.nan,np.nan)
            rows.append({"sample_id":sid,"GSE":f1.loc[sid,"GSE"],"SRR":f1.loc[sid,"SRR"],"arm":arm,"n_sites_evaluated":n_raw,"n_sites":len(af),"mean_fvb_fraction":af.mean(),"median_fvb_fraction":af.median(),
                         "frac_sites_lt_0.5":(af<0.5).mean(),"t_vs_0.5_within_sample":t,"p_vs_0.5_within_sample":pv})
    if not rows: print(f"[{tag}] no tables yet"); return
    ps=pd.DataFrame(rows); ps.to_csv(f"{OUT}/f1_ase_per_sample_v2{suffix}.tsv", sep="\t", index=False)
    summ=[]
    for arm in ("naive","wasp"):
        a=ps[ps.arm==arm]
        if len(a)>1:
            t,pv=stats.ttest_1samp(a.mean_fvb_fraction,0.5); summ.append({"GSE":"all","arm":arm,"n_samples":len(a),"mean_fvb_fraction":a.mean_fvb_fraction.mean(),"se":a.mean_fvb_fraction.sem(),"t_vs_0.5":t,"p_vs_0.5":pv,"sites_min":int(a.n_sites.min()),"sites_max":int(a.n_sites.max()),"sites_mean":a.n_sites.mean()})
        for g,gg in a.groupby("GSE"):
            t,pv=stats.ttest_1samp(gg.mean_fvb_fraction,0.5) if len(gg)>1 else (np.nan,np.nan)
            summ.append({"GSE":g,"arm":arm,"n_samples":len(gg),"mean_fvb_fraction":gg.mean_fvb_fraction.mean(),"se":gg.mean_fvb_fraction.sem(),"t_vs_0.5":t,"p_vs_0.5":pv,"sites_min":int(gg.n_sites.min()),"sites_max":int(gg.n_sites.max()),"sites_mean":gg.n_sites.mean()})
    sm=pd.DataFrame(summ); sm.to_csv(f"{OUT}/f1_ase_summary_v2{suffix}.tsv", sep="\t", index=False); print(f"[{tag}] {len(ps)} sample x arm rows"); print(sm.round(4).to_string(index=False))
summarise("", "QC site set (203,510 hom-alt sites)"); summarise(".oldsites", "old site set (94,351 MGP-het sites)")
