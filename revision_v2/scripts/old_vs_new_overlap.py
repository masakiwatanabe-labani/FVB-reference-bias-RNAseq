#!/usr/bin/env python3
"""Overlap between the frozen v1 catalog and the v2 catalog — for the response to reviewers ONLY (decision 2026-10-20, item 3:
"not in the main text"). This script intentionally reads the frozen v1 table legacy/results/stage6/catalog_final.tsv; it is not
analysis code (rule 0.4-7 concerns analysis_v2/), writes nothing under analysis_v2/, and no v1 number from it is reused in new text.
Output: revision_v2/audit/old_vs_new_catalog_overlap.tsv (one row per v1 gene with its v2 status) and a printed summary."""
import pandas as pd, numpy as np
PROJ="/usr/local/jupyter/FVB_bias"
old=pd.read_csv(f"{PROJ}/legacy/results/stage6/catalog_final.tsv", sep="\t", index_col=0)
new=pd.read_csv(f"{PROJ}/analysis_v2/results/catalog/catalog.primary.tsv", sep="\t", index_col=0)
wide=pd.read_csv(f"{PROJ}/analysis_v2/results/bias/all_tissues.primary.tsv", sep="\t").pivot(index="gene_id",columns="tissue",values="bias")
meta=pd.read_csv(f"{PROJ}/analysis_v2/results/gene_meta.tsv", sep="\t", index_col=0)
rows=[]
for g in old.index:
    tested=g in wide.index; b=wide.loc[g] if tested else None
    rows.append({"gene_id":g,"gene_name_v2":meta.gene_name.get(g,old.loc[g,"gene_name"]),"in_v2_catalog":"Y" if g in new.index else "N",
                 "v2_tested_tissues":int(b.notna().sum()) if tested else 0,"v2_tissues_bias_ge_0.5":int((b>=0.5).sum()) if tested else 0,
                 "v2_max_bias":float(b.max()) if tested else np.nan,"v2_mean_bias":float(b.mean()) if tested else np.nan,"in_fvb_gtf":meta.in_fvb_gtf.get(g,"?")})
d=pd.DataFrame(rows).sort_values(["in_v2_catalog","v2_tissues_bias_ge_0.5"],ascending=False)
d.to_csv(f"{PROJ}/revision_v2/audit/old_vs_new_catalog_overlap.tsv", sep="\t", index=False)
n_old=len(old); n_both=int((d.in_v2_catalog=="Y").sum()); n_tested=int((d.v2_tested_tissues>0).sum()); n_ge1=int((d["v2_tissues_bias_ge_0.5"]>=1).sum())
print(f"v1 catalog genes: {n_old}; in v2 catalog: {n_both}; tested in v2 (>= 1 tissue): {n_tested}; bias >= 0.5 in >= 1 v2 tissue: {n_ge1}; "
      f"not in FVB GTF: {int((d.in_fvb_gtf!='Y').sum())}; v2 catalog size: {len(new)}; v2-only: {len(new)-n_both}")
print(d.to_string(index=False))
