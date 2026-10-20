#!/usr/bin/env python3
"""F1 allele-specific expression summary (Methods: sites with total >= 20 and ref,alt >= 3; FVB fraction = alt/total;
per-sample mean; one-sample t-test vs 0.5; group mean and 95% CI for naive and WASP).
Input: analysis_v2/results/f1_ase/<GSE>/<SRR>.{naive,wasp}.ase.tsv ; sample names from the manifest (F1 rows).
Output: analysis_v2/results/f1_ase/f1_ase_per_sample.tsv, f1_ase_summary.tsv"""
import glob, os, csv, numpy as np, pandas as pd
from scipy import stats
PROJ="/usr/local/jupyter/FVB_bias"; D=f"{PROJ}/analysis_v2/results/f1_ase"
man={r["SRR"]:r for r in csv.DictReader(open(f"{PROJ}/revision_v2/manifest/manifest.tsv"), delimiter="\t") if r["role"]=="f1_ase"}
rows=[]
for p in sorted(glob.glob(f"{D}/GSE*/*.ase.tsv")):
    gse=p.split("/")[-2]; srr,arm=os.path.basename(p).split(".")[:2]
    a=pd.read_csv(p, sep="\t"); tot=a.refCount+a.altCount; a=a[(tot>=20)&(a.refCount>=3)&(a.altCount>=3)]
    af=a.altCount/(a.refCount+a.altCount); t,pv=stats.ttest_1samp(af,0.5)
    rows.append({"GSE":gse,"SRR":srr,"sample_id":man.get(srr,{}).get("sample_id",srr),"arm":arm,"n_sites":len(af),"mean_fvb_fraction":af.mean(),"median_fvb_fraction":af.median(),"t_vs_0.5":t,"p_vs_0.5":pv})
ps=pd.DataFrame(rows); ps.to_csv(f"{D}/f1_ase_per_sample.tsv", sep="\t", index=False)
summ=[]
for (gse,arm),g in ps.groupby(["GSE","arm"]):
    m=g.mean_fvb_fraction; ci=stats.t.interval(0.95, len(m)-1, loc=m.mean(), scale=stats.sem(m)) if len(m)>1 else (np.nan,np.nan)
    t,pv=stats.ttest_1samp(m,0.5) if len(m)>1 else (np.nan,np.nan)
    summ.append({"GSE":gse,"arm":arm,"n_samples":len(m),"mean_fvb_fraction":m.mean(),"ci95_low":ci[0],"ci95_high":ci[1],"t_vs_0.5":t,"p_vs_0.5":pv,"mean_n_sites":g.n_sites.mean()})
pd.DataFrame(summ).to_csv(f"{D}/f1_ase_summary.tsv", sep="\t", index=False); print(pd.DataFrame(summ).to_string())
