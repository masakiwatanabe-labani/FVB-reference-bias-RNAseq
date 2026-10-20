#!/usr/bin/env python3
"""Rule-based choice of the two worked-example genes for Figs 4 and 6 (no gene names are typed anywhere; rule 0.4-8).
  robust  = the catalog gene with the largest number of tissues meeting the criterion, ties broken by mean bias (most reproducible);
  extreme = the catalog gene of the OTHER mechanism class than the robust gene (SNP-induced multi-mapping vs direct SNP
            misalignment, catalog column mechanism_class) with the largest single-tissue bias, and that tissue; if no catalog gene
            of the other class exists, any catalog gene other than the robust one (decision 2026-10-20, item 5: show one example
            of each mechanism).
If the catalog is empty (partial dry runs with < 2 tissues) the same rules are applied to all tested genes and the figure is marked
TEST ONLY. The choice is written to <FIG_DIR>/worked_examples.tsv so the Results text can cite it."""
from fig_common import *
def pick_examples(design=DESIGN):
    long,tissues=load_bias(design); cat=load_catalog(design); meta=gene_meta(); out={"test_only":False}
    if len(cat):
        robust=cat.sort_values(["n_tissues","mean_bias"],ascending=False).index[0]
        sub=long[long.gene_id.isin(cat.index)]
    else:
        out["test_only"]=True; wide=long.pivot(index="gene_id",columns="tissue",values="bias")
        robust=((wide>=BIAS_THR).sum(axis=1)*1000+wide.mean(axis=1).fillna(-9)).sort_values(ascending=False).index[0]; sub=long
    # extreme is chosen among catalog genes OTHER than the robust gene (rule 2026-10-19 22:34), preferring the other mechanism
    # class (rule 2026-10-20: one multi-mapping-type and one direct-misalignment-type example; both recorded in CHANGELOG)
    pool=sub[sub.gene_id!=robust]
    if len(cat) and "mechanism_class" in cat.columns:
        other=pool[pool.gene_id.map(cat["mechanism_class"])!=cat.loc[robust,"mechanism_class"]]
        if len(other): pool=other
    ext=pool.sort_values("bias",ascending=False).iloc[0]
    out.update({"robust":robust,"robust_name":meta.loc[robust,"gene_name"],"extreme":ext.gene_id,"extreme_name":meta.loc[ext.gene_id,"gene_name"],"extreme_tissue":ext.tissue,"extreme_bias":float(ext.bias),
                "robust_mechanism":cat["mechanism_class"].get(robust,"") if len(cat) else "","extreme_mechanism":cat["mechanism_class"].get(ext.gene_id,"") if len(cat) else ""})
    pd.DataFrame([out]).to_csv(f"{FIG}/worked_examples.tsv", sep="\t", index=False)
    return out
if __name__=="__main__": print(pick_examples())
