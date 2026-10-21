#!/usr/bin/env python3
"""Phase K-1(a): per-strain mean CPM of Rps2 (ENSMUSG00000044533) and the three processed pseudogenes named in
revision_v2/audit/rps2_indels_pseudogenes.md (Rps2-ps10, Rps2-ps13, Rps2-ps6) in skin (GSE86315), under unique-only and
-M --fraction counting, naive (B6 reference) and FVB-reference arms; shared-denominator normalisation (every arm divided by the
naive-arm unique-only library size of the same animal). Output: analysis_v2/results/worked_examples/rps2_multimap_cpm.tsv"""
import pandas as pd
R="/usr/local/jupyter/FVB_bias/analysis_v2/results"
gm=pd.read_csv(f"{R}/gene_meta.tsv",sep="\t",index_col=0); genes={"ENSMUSG00000044533":"Rps2"}
genes.update({gid:gm.loc[gid,"gene_name"] for gid in gm.index if gm.loc[gid,"gene_name"] in ("Rps2-ps10","Rps2-ps13","Rps2-ps6")})
ss=pd.read_csv(f"{R}/counts/skin.samples.tsv",sep="\t"); fvb=ss[ss.strain.str.startswith("FVB")].sample_id.tolist(); b6=ss[ss.strain.str.startswith("C57")].sample_id.tolist()
lib=pd.read_csv(f"{R}/counts/skin.naive.unique.tsv",sep="\t",index_col=0).sum()
rows=[]
for arm in ("naive","fvb"):
    for mode in ("unique","multi"):
        c=pd.read_csv(f"{R}/counts/skin.{arm}.{mode}.tsv",sep="\t",index_col=0); cpm=c.div(lib,axis=1)*1e6
        for gid,name in genes.items():
            rows.append({"gene_id":gid,"gene":name,"arm":arm,"counting":"unique-only" if mode=="unique" else "-M --fraction","FVBNJ_mean_cpm":cpm.loc[gid,fvb].mean(),"C57BL6J_mean_cpm":cpm.loc[gid,b6].mean(),"FVBNJ_mean_count":c.loc[gid,fvb].mean(),"C57BL6J_mean_count":c.loc[gid,b6].mean(),"n_FVBNJ":len(fvb),"n_C57BL6J":len(b6)})
df=pd.DataFrame(rows); df.to_csv(f"{R}/worked_examples/rps2_multimap_cpm.tsv",sep="\t",index=False)
pd.set_option("display.width",200); print(df.drop(columns=["gene_id","n_FVBNJ","n_C57BL6J"]).round(2).to_string(index=False))
r=df[df.gene=="Rps2"].set_index(["arm","counting"])
for arm in ("naive","fvb"):
    u=r.loc[(arm,"unique-only"),"FVBNJ_mean_cpm"]; m=r.loc[(arm,"-M --fraction"),"FVBNJ_mean_cpm"]; print(f"{arm}: Rps2 FVB/NJ unique {u:.2f} -> multi {m:.2f} CPM; gain {m-u:+.2f} CPM, {m/u:.2f}-fold")
