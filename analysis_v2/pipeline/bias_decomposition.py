#!/usr/bin/env python3
"""Phase H-2a (2026-10-20): decomposition of the primary bias into its two components, for every tested gene x tissue and for the
138 catalog genes.
  bias          = log2FC_FVBref - log2FC_naive                      (primary definition; catalog criterion, unchanged)
  bias_matched  = log2FC_matched - log2FC_naive                     (FVB/NJ animals on the FVB reference vs C57BL/6J on B6: the
                                                                     gain of FVB/NJ reads on their own reference = the
                                                                     underestimation a naive-B6 user actually experiences)
  b6_loss       = bias - bias_matched                                (loss of C57BL/6J reads on the FVB/NJ reference)
Inputs: analysis_v2/results/sensitivity_matched/bias/all_tissues.matched.tsv (followup_analysis.py section 4),
        analysis_v2/results/catalog/catalog.primary.tsv, sensitivity_matched/catalog.matched.tsv, gene_meta.tsv.
Outputs: analysis_v2/results/sensitivity_matched/bias/all_tissues.decomposition.tsv (all genes x tissues),
         analysis_v2/results/catalog/catalog_decomposition.tsv (138 genes: per-tissue and mean components, in_matched_catalog),
         analysis_v2/results/catalog/decomposition_summary.tsv (counts + the two worked-example genes per tissue).
Nothing in the primary catalog changes."""
import pandas as pd, numpy as np
PROJ="/usr/local/jupyter/FVB_bias"; RES=f"{PROJ}/analysis_v2/results"; SM=f"{RES}/sensitivity_matched"
m=pd.read_csv(f"{SM}/bias/all_tissues.matched.tsv", sep="\t"); m["b6_loss"]=m.bias_primary-m.bias_matched
m=m.rename(columns={"bias_primary":"bias"})[["gene_id","tissue","logFC_naive","logFC_matched","bias","bias_matched","b6_loss"]]
m.to_csv(f"{SM}/bias/all_tissues.decomposition.tsv", sep="\t", index=False)
cat=pd.read_csv(f"{RES}/catalog/catalog.primary.tsv", sep="\t", index_col=0); cm=pd.read_csv(f"{SM}/catalog.matched.tsv", sep="\t", index_col=0)
tissues=["adipose","hypothalamus","liver","skin"]
wb=m.pivot(index="gene_id",columns="tissue",values="bias").reindex(cat.index); wm=m.pivot(index="gene_id",columns="tissue",values="bias_matched").reindex(cat.index); wl=m.pivot(index="gene_id",columns="tissue",values="b6_loss").reindex(cat.index)
out=pd.DataFrame({"gene_name":cat.gene_name,"biotype":cat.biotype,"n_tissues":cat.n_tissues,"mean_bias":cat.mean_bias,"mechanism_class":cat.mechanism_class})
for t in tissues: out[f"{t}_bias"]=wb[t]; out[f"{t}_bias_matched"]=wm[t]; out[f"{t}_b6_loss"]=wl[t]
out["bias_matched_mean"]=wm.mean(axis=1); out["b6_loss_mean"]=wl.mean(axis=1); out["n_tissues_matched_ge_0.5"]=(wm>=0.5).sum(axis=1)
out["in_matched_catalog"]=np.where(out.index.isin(cm.index),"Y","N"); out["b6_loss_share_of_bias"]=out.b6_loss_mean/out.mean_bias
out.index.name="gene_id"; out.sort_values("mean_bias",ascending=False).to_csv(f"{RES}/catalog/catalog_decomposition.tsv", sep="\t")
both=out.index[out.in_matched_catalog=="Y"]
rows=[{"item":"primary catalog","value":len(out)},{"item":"matched-reference catalog","value":len(cm)},{"item":"in both (supported under both definitions)","value":len(both)},
      {"item":"primary only","value":int((out.in_matched_catalog=="N").sum())},{"item":"matched only","value":len(cm.index.difference(out.index))},
      {"item":"primary-only genes with b6_loss_mean > bias_matched_mean","value":int(((out.in_matched_catalog=="N")&(out.b6_loss_mean>out.bias_matched_mean)).sum())},
      {"item":"catalog median share of bias due to B6-read loss (b6_loss_mean/mean_bias)","value":round(float(out.b6_loss_share_of_bias.median()),3)},
      {"item":"primary-only median share","value":round(float(out.b6_loss_share_of_bias[out.in_matched_catalog=="N"].median()),3)},
      {"item":"both median share","value":round(float(out.b6_loss_share_of_bias[out.in_matched_catalog=="Y"].median()),3)},
      {"item":"genome-wide median b6_loss (all tested gene x tissue)","value":round(float(m.b6_loss.median()),4)},{"item":"genome-wide median bias_matched","value":round(float(m.bias_matched.median()),4)}]
for g,name in (("ENSMUSG00000028648","Ndufs5"),("ENSMUSG00000044533","Rps2")):
    for t in tissues:
        r=m[(m.gene_id==g)&(m.tissue==t)]
        if len(r): r=r.iloc[0]; rows.append({"item":f"{name} {t}: bias / bias_matched / b6_loss","value":f"{r.bias:+.2f} / {r.bias_matched:+.2f} / {r.b6_loss:+.2f}"})
    rows.append({"item":f"{name} meets the catalog criterion under the matched definition (>= 0.5 in >= 2 tissues)","value":"Y" if g in cm.index else "N"})
    rows.append({"item":f"{name} tissues with bias_matched >= 0.5","value":int(out.loc[g,"n_tissues_matched_ge_0.5"])})
pd.DataFrame(rows).to_csv(f"{RES}/catalog/decomposition_summary.tsv", sep="\t", index=False); print(pd.DataFrame(rows).to_string(index=False))
