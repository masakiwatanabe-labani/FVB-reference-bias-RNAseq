#!/usr/bin/env python3
"""Phase H-1b supplement: at the per-site-filtered F1 sites (total >= 20, ref >= 3, alt >= 3 in the naive arm), what fraction of the
reference-allele (B6) reads and of the alternative-allele (FVB) reads does the WASP filter (vW == 1 or untagged) remove, per animal?
Explains the direction of the WASP arm in Fig. 2 (per-site FVB fraction rises above 0.50 because WASP removes B6-allele reads whose
allele-swapped version fails to remap, while FVB-allele reads that failed to align in the naive arm cannot be recovered).
Input: analysis_v2/results/f1_ase/v2/<sid>.{naive,wasp}.ase.tsv. Output: analysis_v2/results/f1_ase/v2/f1_wasp_allele_removal.tsv"""
import pandas as pd, glob, os
D="/usr/local/jupyter/FVB_bias/analysis_v2/results/f1_ase/v2"; rows=[]
for p in sorted(glob.glob(f"{D}/*.naive.ase.tsv")):
    sid=os.path.basename(p).split(".naive")[0]
    n=pd.read_csv(p,sep="\t"); w=pd.read_csv(p.replace(".naive.",".wasp."),sep="\t")
    n=n[(n.totalCount>=20)&(n.refCount>=3)&(n.altCount>=3)]
    m=n.merge(w,on=["contig","position"],suffixes=("_n","_w"),how="left").fillna(0)
    af_n=m.altCount_n/m.totalCount_n; af_w=m.altCount_w/m.totalCount_w.replace(0,float("nan"))
    rows.append({"sample_id":sid,"n_sites":len(m),"ref_reads_naive":int(m.refCount_n.sum()),"alt_reads_naive":int(m.altCount_n.sum()),
                 "frac_ref_reads_removed_by_wasp":(m.refCount_n-m.refCount_w).sum()/m.refCount_n.sum(),"frac_alt_reads_removed_by_wasp":(m.altCount_n-m.altCount_w).sum()/m.altCount_n.sum(),
                 "frac_sites_af_increased":float((af_w>af_n+1e-9).mean()),"frac_sites_af_unchanged":float(((af_w-af_n).abs()<=1e-9).mean()),"frac_sites_af_decreased":float((af_w<af_n-1e-9).mean()),
                 "pooled_af_naive":m.altCount_n.sum()/m.totalCount_n.sum(),"pooled_af_wasp":m.altCount_w.sum()/m.totalCount_w.sum()})
df=pd.DataFrame(rows); df.to_csv(f"{D}/f1_wasp_allele_removal.tsv",sep="\t",index=False); print(df.round(4).to_string(index=False))
