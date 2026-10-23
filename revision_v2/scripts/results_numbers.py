#!/usr/bin/env python3
"""Every number quoted in the revised text (manuscript/revised_v2/*.md, response letter) is derived here from analysis_v2/results
and analysis_v2/figures and written to manuscript/revised_v2/generated/numbers.tsv (key, value, evidence file) — rule 0.4-3.
Nothing is read from legacy/. Keys are referenced in revision_v2/text_changes.md next to each new sentence."""
import os, glob, collections, re
import numpy as np, pandas as pd
from scipy import stats
PROJ="/usr/local/jupyter/FVB_bias"; RES=f"{PROJ}/analysis_v2/results"; FIG=f"{PROJ}/analysis_v2/figures"; OUT=f"{PROJ}/manuscript/revised_v2/generated/numbers.tsv"
N=[]
def put(key, value, evidence):
    N.append({"key":key,"value":value if isinstance(value,str) else (f"{value:.6g}" if isinstance(value,(float,np.floating)) else str(value)),"evidence":evidence})
rel=lambda p: p.replace(PROJ+"/","")
# ---- datasets / samples ---------------------------------------------------------------------------------------------
man=pd.read_csv(f"{PROJ}/revision_v2/manifest/manifest.tsv", sep="\t", dtype=str); inc=man[man["include"]=="Y"]
sc=inc[inc.role=="strain_comparison"]; put("n_animals_strain_comparison", len(sc), "revision_v2/manifest/manifest.tsv")
put("n_animals_strain_comparison_FVB", int((sc.strain=="FVB/NJ").sum()), "manifest"); put("n_animals_strain_comparison_B6", int((sc.strain=="C57BL/6J").sum()), "manifest")
for t,g in sc.groupby("tissue"): put(f"n_{t}_FVB", int((g.strain=="FVB/NJ").sum()), "manifest"); put(f"n_{t}_B6", int((g.strain=="C57BL/6J").sum()), "manifest")
put("n_pyvt", int((inc.role=="negative_control_candidate").sum()), "manifest"); put("n_f1", int((inc.role=="f1_ase").sum()), "manifest")
put("n_manifest_rows", len(man), "manifest"); put("n_manifest_excluded", int((man["include"]!="Y").sum()), "manifest")
put("n_animals_processed_three_arms", len(sc)+int((inc.role=="negative_control_candidate").sum()), "manifest")
# ---- identity QC ----------------------------------------------------------------------------------------------------
q=pd.read_csv(f"{RES}/qc/sample_identity_verification.tsv", sep="\t"); ev=rel(f"{RES}/qc/sample_identity_verification.tsv")
put("qc_rows", len(q), ev); put("qc_used_rows", int((q.used_in_analysis=="Y").sum()), ev); put("qc_used_match", int(((q.used_in_analysis=="Y")&(q["match"]=="Y")).sum()), ev)
for v,g in q.groupby("qc_verdict"): put(f"qc_verdict_{v.replace(' ','_').replace('/','').replace('(','').replace(')','')}", len(g), ev)
pure=q[q.qc_verdict.isin(["B6 pure","FVB/NJ pure"])&(q.used_in_analysis=="Y")]
put("qc_pure_n_sites_min", int(pure.n_sites.min()), ev); put("qc_pure_n_sites_max", int(pure.n_sites.max()), ev); put("qc_pure_n_sites_median", float(pure.n_sites.median()), ev)
b6=q[q.qc_verdict=="B6 pure"]; fv=q[q.qc_verdict=="FVB/NJ pure"]; f1=q[q.qc_verdict.str.startswith("F1")]
put("qc_B6_mean_af_max", float(b6.mean_alt_fraction.max()), ev); put("qc_B6_frac_lo_min", float(b6.frac_lo.min()), ev)
put("qc_FVB_mean_af_min", float(fv.mean_alt_fraction.min()), ev); put("qc_FVB_frac_hi_min", float(fv.frac_hi.min()), ev)
put("qc_F1_mean_af_min", float(f1.mean_alt_fraction.min()), ev); put("qc_F1_mean_af_max", float(f1.mean_alt_fraction.max()), ev); put("qc_F1_frac_mid_min", float(f1.frac_mid.min()), ev); put("qc_F1_n_sites_min", int(f1.n_sites.min()), ev); put("qc_F1_n_sites_max", int(f1.n_sites.max()), ev)
hx=q[q.used_in_analysis=="N"]; put("qc_failed_n", len(hx), ev); put("qc_failed_GSE", ",".join(sorted(set(hx.GSE))), ev)
put("qc_failed_mean_af_min", float(hx.mean_alt_fraction.min()), ev); put("qc_failed_mean_af_max", float(hx.mean_alt_fraction.max()), ev); put("qc_failed_frac_hi_min", float(hx.frac_hi.min()), ev); put("qc_failed_frac_hi_max", float(hx.frac_hi.max()), ev)
put("qc_failed_frac_mid_min", float(hx.frac_mid.min()), ev); put("qc_failed_frac_mid_max", float(hx.frac_mid.max()), ev); put("qc_failed_n_sites_min", int(hx.n_sites.min()), ev); put("qc_failed_n_sites_max", int(hx.n_sites.max()), ev)
put("qc_failed_verdicts", ";".join(sorted(set(hx.qc_verdict))), ev)
# ---- survey (Fig. 1) --------------------------------------------------------------------------------------------------
s1=pd.read_csv(f"{FIG}/Fig1_survey_data.tsv", sep="\t"); ev=rel(f"{FIG}/Fig1_survey_data.tsv")
for _,r in s1.iterrows(): put(f"survey_{r.category}_n", int(r.n), ev); put(f"survey_{r.category}_pct", float(r.pct), ev)
put("survey_total", int(s1.n.sum()), ev)
# ---- F1 ASE (Fig. 2) -------------------------------------------------------------------------------------------------
# Phase H-1 (2026-10-21): recomputed on the re-aligned F1 BAMs at the 203,510 identity-QC sites (f1_ase_rerun.sh, f1_ase_summary_v2.py);
# the v1 tables in f1_ase/ (94,351 MGP-heterozygous sites, an artefact of a genotype-unfiltered VCF; H-1a) are no longer read.
fa=pd.read_csv(f"{RES}/f1_ase/v2/f1_ase_per_sample_v2.tsv", sep="\t"); ev=rel(f"{RES}/f1_ase/v2/f1_ase_per_sample_v2.tsv")
for arm in ("naive","wasp"):
    a=fa[fa.arm==arm]; m=a.mean_fvb_fraction; t,p=stats.ttest_1samp(m,0.5)
    put(f"f1_{arm}_n_samples", len(a), ev); put(f"f1_{arm}_pooled_mean", float(m.mean()), ev); put(f"f1_{arm}_pooled_se", float(m.sem()), ev); put(f"f1_{arm}_t", float(t), ev); put(f"f1_{arm}_p", float(p), ev)
    put(f"f1_{arm}_per_animal_min", float(m.min()), ev); put(f"f1_{arm}_per_animal_max", float(m.max()), ev)
    put(f"f1_{arm}_sites_min", int(a.n_sites.min()), ev); put(f"f1_{arm}_sites_max", int(a.n_sites.max()), ev); put(f"f1_{arm}_sites_mean", float(a.n_sites.mean()), ev)
    for g,gg in a.groupby("GSE"): put(f"f1_{arm}_{g}_mean", float(gg.mean_fvb_fraction.mean()), ev); put(f"f1_{arm}_{g}_n", len(gg), ev)
put("f1_naive_deficit_pp", float((0.5-fa[fa.arm=="naive"].mean_fvb_fraction.mean())*100), ev)
put("f1_wasp_excess_pp", float((fa[fa.arm=="wasp"].mean_fvb_fraction.mean()-0.5)*100), ev)
put("f1_site_set_n", 203510, "analysis_v2/qc/fvbnj_homalt_exonic.hetGT.vcf.gz")
wr=pd.read_csv(f"{RES}/f1_ase/v2/f1_wasp_allele_removal.tsv", sep="\t"); ev=rel(f"{RES}/f1_ase/v2/f1_wasp_allele_removal.tsv")
put("f1_wasp_ref_removed_pct_min", float(wr.frac_ref_reads_removed_by_wasp.min()*100), ev); put("f1_wasp_ref_removed_pct_max", float(wr.frac_ref_reads_removed_by_wasp.max()*100), ev)
put("f1_wasp_alt_removed_pct_min", float(wr.frac_alt_reads_removed_by_wasp.min()*100), ev); put("f1_wasp_alt_removed_pct_max", float(wr.frac_alt_reads_removed_by_wasp.max()*100), ev)
put("f1_wasp_sites_af_increased_pct_min", float(wr.frac_sites_af_increased.min()*100), ev); put("f1_wasp_sites_af_increased_pct_max", float(wr.frac_sites_af_increased.max()*100), ev)
put("f1_wasp_sites_af_decreased_pct_max", float(wr.frac_sites_af_decreased.max()*100), ev)
# old site set (for the response letter only): the same new alignments evaluated at the 94,351 v1 sites reproduce the v1 values
fo=pd.read_csv(f"{RES}/f1_ase/v2/f1_ase_per_sample_v2.oldsites.tsv", sep="\t"); ev=rel(f"{RES}/f1_ase/v2/f1_ase_per_sample_v2.oldsites.tsv")
for arm in ("naive","wasp"):
    a=fo[fo.arm==arm]; m=a.mean_fvb_fraction; put(f"f1_oldsites_{arm}_pooled_mean", float(m.mean()), ev); put(f"f1_oldsites_{arm}_sites_min", int(a.n_sites.min()), ev); put(f"f1_oldsites_{arm}_sites_max", int(a.n_sites.max()), ev)
put("f1_oldsites_n_vcf", 94351, "revision_v2/audit/f1_old_site_genotypes.md")
# F1 sites inside catalog genes vs the other evaluated sites (H-1d item 7, adopted 2026-10-21)
fc=pd.read_csv(f"{RES}/f1_ase/v2/f1_sites_catalog_vs_other.tsv", sep="\t"); ev=rel(f"{RES}/f1_ase/v2/f1_sites_catalog_vs_other.tsv")
for arm in ("naive","wasp"):
    a=fc[fc.arm==arm]; t,p=stats.ttest_rel(a.mean_af_catalog, a.mean_af_other)
    put(f"f1_catalog_sites_{arm}_mean_af", float(a.mean_af_catalog.mean()), ev); put(f"f1_other_sites_{arm}_mean_af", float(a.mean_af_other.mean()), ev)
    put(f"f1_catalog_vs_other_{arm}_t", float(t), ev); put(f"f1_catalog_vs_other_{arm}_p", float(p), ev)
    put(f"f1_catalog_sites_{arm}_n_min", int(a.n_sites_catalog.min()), ev); put(f"f1_catalog_sites_{arm}_n_max", int(a.n_sites_catalog.max()), ev)
put("f1_catalog_sites_naive_deficit_pp", float((0.5-fc[fc.arm=="naive"].mean_af_catalog.mean())*100), ev)
# Rps2 fact check (revision_v2/audit/rps2_indels_pseudogenes.md; counts table rps2_pseudogene_counts_skin.tsv)
ev="revision_v2/audit/rps2_indels_pseudogenes.md"; put("rps2_indels_gene_body", 12, ev); put("rps2_indels_exonic", 4, ev); put("rps2_snps_exonic", 13, ev); put("rps2_n_processed_pseudogenes", 12, ev)
# ---- genome-wide bias per tissue (Fig. 3) ----------------------------------------------------------------------------
long=pd.read_csv(f"{RES}/bias/all_tissues.primary.tsv", sep="\t"); ev=rel(f"{RES}/bias/all_tissues.primary.tsv"); tissues=["adipose","hypothalamus","liver","skin"]
wide=long.pivot(index="gene_id",columns="tissue",values="bias"); put("genes_tested_union", len(wide), ev); put("genes_tested_all4", int(wide.notna().all(axis=1).sum()), ev)
for t in tissues:
    b=long[long.tissue==t]; put(f"gw_{t}_n_genes", len(b), ev); put(f"gw_{t}_median_bias", float(b.bias.median()), ev); put(f"gw_{t}_mean_bias", float(b.bias.mean()), ev)
    put(f"gw_{t}_n_ge_0.5", int((b.bias>=0.5).sum()), ev); put(f"gw_{t}_pct_ge_0.5", float((b.bias>=0.5).mean()*100), ev); put(f"gw_{t}_n_le_-0.5", int((b.bias<=-0.5).sum()), ev)
    put(f"gw_{t}_median_bias_wasp", float(b.bias_wasp.median()), ev); put(f"gw_{t}_mean_logFC_naive", float(b.logFC_naive.mean()), ev); put(f"gw_{t}_mean_logFC_wasp", float(b.logFC_wasp.mean()), ev); put(f"gw_{t}_mean_logFC_fvb", float(b.logFC_fvb.mean()), ev)
    put(f"gw_{t}_wasp_minus_naive_mean_logFC", float(b.logFC_wasp.mean()-b.logFC_naive.mean()), ev); put(f"gw_{t}_fvb_minus_naive_mean_logFC", float(b.logFC_fvb.mean()-b.logFC_naive.mean()), ev)
put("gw_median_bias_absmax", float(max(abs(long[long.tissue==t].bias.median()) for t in tissues)), ev)
put("gw_mean_bias_min", float(min(long[long.tissue==t].bias.mean() for t in tissues)), ev); put("gw_mean_bias_max", float(max(long[long.tissue==t].bias.mean() for t in tissues)), ev)
put("gw_pct_ge_0.5_min", float(min((long[long.tissue==t].bias>=0.5).mean()*100 for t in tissues)), ev); put("gw_pct_ge_0.5_max", float(max((long[long.tissue==t].bias>=0.5).mean()*100 for t in tissues)), ev)
put("gw_wasp_minus_naive_absmax", float(max(abs(long[long.tissue==t].logFC_wasp.mean()-long[long.tissue==t].logFC_naive.mean()) for t in tissues)), ev)
put("gw_fvb_minus_naive_min", float(min(long[long.tissue==t].logFC_fvb.mean()-long[long.tissue==t].logFC_naive.mean() for t in tissues)), ev); put("gw_fvb_minus_naive_max", float(max(long[long.tissue==t].logFC_fvb.mean()-long[long.tissue==t].logFC_naive.mean() for t in tissues)), ev)
sens=pd.read_csv(f"{RES}/catalog/catalog_sensitivity.primary.tsv", sep="\t"); ev2=rel(f"{RES}/catalog/catalog_sensitivity.primary.tsv")
for _,r in sens.iterrows(): put(f"sens_thr{r.threshold}_min{int(r.min_tissues)}", int(r.n_genes_in_fvb_gtf), ev2)
n1=int(sens[(sens.threshold==0.5)&(sens.min_tissues==1)].n_genes_in_fvb_gtf.iloc[0]); put("genes_ge_0.5_any_tissue", n1, ev2); put("genes_ge_0.5_any_tissue_pct", n1/len(wide)*100, ev2)
sw=pd.read_csv(f"{RES}/catalog/catalog_sensitivity.water.tsv", sep="\t"); ev3=rel(f"{RES}/catalog/catalog_sensitivity.water.tsv")
for _,r in sw.iterrows(): put(f"sens_water_thr{r.threshold}_min{int(r.min_tissues)}", int(r.n_genes_in_fvb_gtf), ev3)
# biotype Kruskal–Wallis test exactly as drawn in Fig. 4D (fig4_predictors.py: biotype_group order, mean_bias dropna, groups with > 1 gene)
prk=pd.read_csv(f"{RES}/predictors/gene_predictors.primary.tsv", sep="\t"); ev=rel(f"{RES}/predictors/gene_predictors.primary.tsv")
kwd=[prk[prk.biotype_group==b].mean_bias.dropna().values for b in ("protein_coding","lncRNA","pseudogene","miRNA","other") if b in set(prk.biotype_group)]
put("pred_biotype_kruskal_p", float(stats.kruskal(*[d for d in kwd if len(d)>1]).pvalue), ev)
# ---- bias decomposition (Phase H-2; bias = bias_matched + b6_loss) -----------------------------------------------------
ds=pd.read_csv(f"{RES}/catalog/decomposition_summary.tsv", sep="\t"); ev=rel(f"{RES}/catalog/decomposition_summary.tsv"); dsv=dict(zip(ds.item, ds.value))
put("decomp_catalog_primary", int(dsv["primary catalog"]), ev); put("decomp_catalog_matched", int(dsv["matched-reference catalog"]), ev); put("decomp_both", int(dsv["in both (supported under both definitions)"]), ev)
put("decomp_primary_only", int(dsv["primary only"]), ev); put("decomp_matched_only", int(dsv["matched only"]), ev); put("decomp_primary_only_b6loss_dominant", int(dsv["primary-only genes with b6_loss_mean > bias_matched_mean"]), ev)
put("decomp_median_b6loss_share_catalog", float(dsv["catalog median share of bias due to B6-read loss (b6_loss_mean/mean_bias)"]), ev); put("decomp_median_b6loss_share_primary_only", float(dsv["primary-only median share"]), ev); put("decomp_median_b6loss_share_both", float(dsv["both median share"]), ev)
put("decomp_gw_median_b6_loss", float(dsv["genome-wide median b6_loss (all tested gene x tissue)"]), ev); put("decomp_gw_median_bias_matched", float(dsv["genome-wide median bias_matched"]), ev)
dc=pd.read_csv(f"{RES}/catalog/catalog_decomposition.tsv", sep="\t").set_index("gene_name"); ev=rel(f"{RES}/catalog/catalog_decomposition.tsv")
for g in ("Ndufs5","Rps2"):
    r=dc.loc[g]; put(f"decomp_{g}_n_tissues_matched_ge_0.5", int(r["n_tissues_matched_ge_0.5"]), ev); put(f"decomp_{g}_in_matched_catalog", str(r["in_matched_catalog"]), ev)
    for t in tissues:
        if pd.notna(r[f"{t}_bias_matched"]): put(f"decomp_{g}_{t}_bias_matched", float(r[f"{t}_bias_matched"]), ev); put(f"decomp_{g}_{t}_b6_loss", float(r[f"{t}_b6_loss"]), ev)
    put(f"decomp_{g}_bias_matched_mean", float(r["bias_matched_mean"]), ev); put(f"decomp_{g}_b6_loss_mean", float(r["b6_loss_mean"]), ev)
# representative scatter tissue
f3=open(f"{FIG}/Fig3_strain_comparison_data.tsv").read().split("## panelA_representative")[1].split("\n")[2].split("\t"); put("fig3A_tissue", f3[1], rel(f"{FIG}/Fig3_strain_comparison_data.tsv")); put("fig3A_n_genes", f3[4], rel(f"{FIG}/Fig3_strain_comparison_data.tsv"))
# ---- WASP read attrition and mapping rates (mapping_stats) ------------------------------------------------------------
ms=pd.read_csv(f"{RES}/qc/mapping_stats_per_sample.tsv", sep="\t"); ev=rel(f"{RES}/qc/mapping_stats_per_sample.tsv")
pct=lambda s: s.str.rstrip("%").astype(float)
nv=ms[ms.arm=="naive"].set_index("sample_id"); wp=ms[ms.arm=="wasp"].set_index("sample_id"); fb=ms[ms.arm=="fvb"].set_index("sample_id")
att=(1-wp.fc_assigned_unique.reindex(nv.index)/nv.fc_assigned_unique)*100; nv=nv.assign(wasp_removed_pct=att, uniq=pct(nv.uniq_pct), fvb_uniq=pct(fb.uniq_pct.reindex(nv.index)))
nv["strain_short"]=nv.strain.map(lambda s:"FVB" if s.startswith("FVB") else "B6")
rows=[]
for (t,s),g in nv.groupby(["tissue","strain_short"]):
    rows.append({"tissue":t,"strain":s,"n":len(g),"wasp_removed_mean":g.wasp_removed_pct.mean(),"wasp_removed_sd":g.wasp_removed_pct.std(),"naive_uniq_mean":g.uniq.mean(),"fvb_uniq_mean":g.fvb_uniq.mean()})
    put(f"wasp_removed_{t}_{s}_mean_pct", g.wasp_removed_pct.mean(), ev); put(f"uniq_naive_{t}_{s}_mean_pct", g.uniq.mean(), ev); put(f"uniq_fvb_{t}_{s}_mean_pct", g.fvb_uniq.mean(), ev)
att_tab=pd.DataFrame(rows); att_tab.to_csv(f"{PROJ}/manuscript/revised_v2/generated/wasp_attrition_by_group.tsv", sep="\t", index=False)
st=att_tab[att_tab.tissue!="mammary_tumour_PyVT"]; put("wasp_removed_min_pct", float(st.wasp_removed_mean.min()), ev); put("wasp_removed_max_pct", float(st.wasp_removed_mean.max()), ev)
for t in tissues:
    d=st[st.tissue==t].set_index("strain"); put(f"wasp_removed_{t}_FVB_minus_B6_pp", float(d.loc["FVB","wasp_removed_mean"]-d.loc["B6","wasp_removed_mean"]), ev)
    put(f"uniq_gap_naive_{t}_FVB_minus_B6_pp", float(d.loc["FVB","naive_uniq_mean"]-d.loc["B6","naive_uniq_mean"]), ev); put(f"uniq_gain_fvbarm_{t}_FVB_pp", float(d.loc["FVB","fvb_uniq_mean"]-d.loc["FVB","naive_uniq_mean"]), ev)
put("wasp_removed_B6_more_than_FVB_tissues", int(sum(1 for t in tissues if st[(st.tissue==t)&(st.strain=="B6")].wasp_removed_mean.iloc[0]>st[(st.tissue==t)&(st.strain=="FVB")].wasp_removed_mean.iloc[0])), ev)
# ---- predictors (Fig. 4) ---------------------------------------------------------------------------------------------
pr=pd.read_csv(f"{RES}/predictors/gene_predictors.primary.tsv", sep="\t", index_col=0); ev=rel(f"{RES}/predictors/gene_predictors.primary.tsv")
cat=pd.read_csv(f"{RES}/catalog/catalog.primary.tsv", sep="\t", index_col=0); evc=rel(f"{RES}/catalog/catalog.primary.tsv")
for c in ("snp_density_per_kb","log10_length","mean_logCPM_naive"):
    rho,p=stats.spearmanr(pr[c],pr.mean_bias,nan_policy="omit"); put(f"pred_{c}_rho", float(rho), ev); put(f"pred_{c}_p", float(p), ev); put(f"pred_{c}_n", int(pr[c].notna().sum()), ev)
all4=pr.loc[wide.index[wide.notna().all(axis=1)]]; rho4,p4=stats.spearmanr(all4.snp_density_per_kb,all4.mean_bias,nan_policy="omit"); put("pred_snp_rho_all4tissue_genes", float(rho4), ev); put("pred_snp_p_all4tissue_genes", float(p4), ev); put("pred_n_all4tissue_genes", len(all4), ev)
for b,g in pr.groupby("biotype_group"): put(f"pred_biotype_{b}_n", len(g), ev); put(f"pred_biotype_{b}_mean", float(g.mean_bias.mean()), ev); put(f"pred_biotype_{b}_median", float(g.mean_bias.median()), ev)
q75=pr.snp_density_per_kb.quantile(0.75); q90=pr.snp_density_per_kb.quantile(0.90); q10=pr.snp_density_per_kb.quantile(0.10)
incat=pr.index.isin(cat.index); put("pred_q75_snp_density", float(q75), ev); put("pred_catalog_in_top_quartile", int((pr.snp_density_per_kb[incat]>=q75).sum()), ev); put("pred_catalog_n_with_density", int(incat.sum()), ev)
put("pred_catalog_median_snp_density", float(pr.snp_density_per_kb[incat].median()), ev); put("pred_noncatalog_median_snp_density", float(pr.snp_density_per_kb[~incat].median()), ev)
u,pu=stats.mannwhitneyu(pr.snp_density_per_kb[incat].dropna(),pr.snp_density_per_kb[~incat].dropna(),alternative="two-sided"); put("pred_catalog_vs_noncatalog_snp_density_MWU_p", float(pu), ev)
top=pr.snp_density_per_kb>=q90; a=int((top&incat).sum()); b_=int((~top&incat).sum()); c_=int((top&~incat).sum()); d_=int((~top&~incat).sum()); OR=(a*d_)/(b_*c_); se=np.sqrt(1/a+1/b_+1/c_+1/d_)
put("pred_q90_snp_density", float(q90), ev); put("pred_top_decile_catalog_n", a, ev); put("pred_top_decile_OR", float(OR), ev); put("pred_top_decile_OR_ci_lo", float(np.exp(np.log(OR)-1.96*se)), ev); put("pred_top_decile_OR_ci_hi", float(np.exp(np.log(OR)+1.96*se)), ev)
put("pred_top_decile_fisher_p", float(stats.fisher_exact([[a,b_],[c_,d_]])[1]), ev)
u2,p2=stats.mannwhitneyu(pr.mean_bias[top],pr.mean_bias[pr.snp_density_per_kb<=q10],alternative="two-sided"); put("pred_top_vs_bottom_decile_bias_MWU_p", float(p2), ev); put("pred_top_decile_mean_bias", float(pr.mean_bias[top].mean()), ev); put("pred_bottom_decile_mean_bias", float(pr.mean_bias[pr.snp_density_per_kb<=q10].mean()), ev)
put("pred_median_snp_density_all_tested", float(pr.snp_density_per_kb.median()), ev)
cc=pr.dropna(subset=["mean_logCPM_naive"]); edges=np.unique(np.percentile(cc.mean_logCPM_naive,np.linspace(0,100,21))); idx=np.digitize(cc.mean_logCPM_naive,edges)
bm=[cc.mean_bias[idx==b].mean() for b in range(1,len(edges)) if (idx==b).sum()>10]; put("pred_expression_binned_mean_absmax", float(max(abs(x) for x in bm)), ev)
put("pred_biotype_median_absmax", float(pr.groupby("biotype_group").mean_bias.median().abs().max()), ev)
kw=stats.kruskal(*[g.mean_bias.dropna().values for _,g in pr.groupby("biotype_group")]); put("pred_biotype_kruskal_p", float(kw.pvalue), ev)
# ---- catalog (Fig. 5) -------------------------------------------------------------------------------------------------
put("catalog_n", len(cat), evc); put("catalog_pct_of_tested", len(cat)/len(wide)*100, evc)
for k in (4,3,2): put(f"catalog_n_tissues_{k}", int((cat.n_tissues==k).sum()), evc)
bt=cat.biotype.map(lambda b: b if b in ("protein_coding","lncRNA") else ("pseudogene" if "pseudogene" in b else "other"))
for b,n in bt.value_counts().items(): put(f"catalog_biotype_{b}", int(n), evc)
put("catalog_other_biotypes", ";".join(sorted(set(cat.biotype[bt=="other"]))), evc)
for t in tissues: put(f"catalog_biased_in_{t}", int((cat[t]>=0.5).sum()), evc)
put("catalog_tested_in_all4", int((cat.n_tissues_tested==4).sum()), evc); put("catalog_mean_bias_lt_0.5", int((cat.mean_bias<0.5).sum()), evc)
put("catalog_mean_bias_max", float(cat.mean_bias.max()), evc); put("catalog_mean_bias_median", float(cat.mean_bias.median()), evc)
put("catalog_single_tissue_bias_max", float(cat[tissues].max().max()), evc); put("catalog_single_tissue_bias_max_tissue", cat[tissues].max().idxmax(), evc); put("catalog_single_tissue_bias_max_gene", cat[tissues].max(axis=1).idxmax(), evc)
cw=pd.read_csv(f"{RES}/catalog/catalog.water.tsv", sep="\t", index_col=0); evw=rel(f"{RES}/catalog/catalog.water.tsv"); put("catalog_water_n", len(cw), evw); put("catalog_water_shared", len(cw.index.intersection(cat.index)), evw)
# mechanism (SFig 1)
mc=cat.mechanism_class; mm=cat[mc=="SNP-induced multi-mapping"]; md=cat[mc=="direct SNP misalignment"]
put("mech_multi_n", len(mm), evc); put("mech_direct_n", len(md), evc); put("mech_multi_mean_absdelta", float(mm.delta_multi_minus_unique.abs().mean()), evc); put("mech_direct_mean_absdelta", float(md.delta_multi_minus_unique.abs().mean()), evc)
put("mech_multi_residual_le_0.30", int((mm.residual_fraction<=0.30).sum()), evc); put("mech_multi_residual_min", float(mm.residual_fraction.min()), evc); put("mech_multi_residual_max", float(mm.residual_fraction.max()), evc)
put("mech_direct_residual_ge_0.80", int((md.residual_fraction>=0.80).sum()), evc); put("mech_direct_residual_min", float(md.residual_fraction.min()), evc); put("mech_direct_residual_max", float(md.residual_fraction.max()), evc)
put("mech_multi_biotype_pseudogene", int(bt[mc=="SNP-induced multi-mapping"].eq("pseudogene").sum()), evc); put("mech_direct_biotype_pseudogene", int(bt[mc=="direct SNP misalignment"].eq("pseudogene").sum()), evc)
sf=open(f"{FIG}/SFig2_multimapping_data.tsv").read().split("## genome_wide")[1].split("\n")[2].split("\t"); evs=rel(f"{FIG}/SFig2_multimapping_data.tsv")
put("mm_n_pairs", int(sf[1]), evs); put("mm_median_delta", float(sf[2]), evs); put("mm_mean_delta", float(sf[3]), evs); put("mm_pct_abs_gt_1", float(sf[4])*100, evs)
pairs=long.dropna(subset=["delta_multi_minus_unique"]); put("mm_pct_abs_gt_0.5", float((pairs.delta_multi_minus_unique.abs()>0.5).mean()*100), ev)
# SFig1 panel D picks
sfd=pd.read_csv(f"{FIG}/SFig2_multimapping_data.tsv", sep="\t", skiprows=lambda i: False) if False else None
txt=open(f"{FIG}/SFig2_multimapping_data.tsv").read().split("## panelD")[1].strip().split("\n"); pd_genes=[]
for l in txt[1:]:
    g=l.split("\t")[1]
    if g not in pd_genes: pd_genes.append(g)
put("sfig1D_genes", ",".join(pd_genes), evs)
# ---- worked examples (Fig. 6) -----------------------------------------------------------------------------------------
we=pd.read_csv(f"{FIG}/worked_examples.tsv", sep="\t").iloc[0]; evw6=rel(f"{FIG}/worked_examples.tsv"); meta=pd.read_csv(f"{RES}/gene_meta.tsv", sep="\t", index_col=0)
for role in ("robust","extreme"):
    g=we[role]; put(f"ex_{role}_gene", we[f"{role}_name"], evw6); put(f"ex_{role}_id", g, evw6); put(f"ex_{role}_mechanism", we[f"{role}_mechanism"], evw6)
    put(f"ex_{role}_chrom", meta.loc[g,"chrom"], "analysis_v2/results/gene_meta.tsv"); put(f"ex_{role}_n_snps", int(meta.loc[g,"n_snps"]), "gene_meta.tsv"); put(f"ex_{role}_length_kb", float(meta.loc[g,"length_bp"])/1000, "gene_meta.tsv"); put(f"ex_{role}_snp_density", float(meta.loc[g,"snp_density_per_kb"]), "gene_meta.tsv")
    put(f"ex_{role}_n_tissues", int(cat.loc[g,"n_tissues"]), evc); put(f"ex_{role}_mean_bias", float(cat.loc[g,"mean_bias"]), evc); put(f"ex_{role}_residual", float(cat.loc[g,"residual_fraction"]), evc); put(f"ex_{role}_biotype", cat.loc[g,"biotype"], evc)
    for t in tissues:
        r=long[(long.gene_id==g)&(long.tissue==t)]
        if len(r):
            r=r.iloc[0]; put(f"ex_{role}_{t}_bias", float(r.bias), ev); put(f"ex_{role}_{t}_logFC_naive", float(r.logFC_naive), ev); put(f"ex_{role}_{t}_logFC_wasp", float(r.logFC_wasp), ev); put(f"ex_{role}_{t}_logFC_fvb", float(r.logFC_fvb), ev); put(f"ex_{role}_{t}_bias_multi", float(r.bias_multi), ev); put(f"ex_{role}_{t}_residual", float(r.residual_fraction), ev); put(f"ex_{role}_{t}_bias_wasp", float(r.bias_wasp), ev)
put("ex_extreme_tissue", we.extreme_tissue, evw6); put("ex_extreme_bias", float(we.extreme_bias), evw6); put("ex_extreme_fold", float(2**we.extreme_bias), evw6)
rx=long[(long.gene_id==we["extreme"])&(long.tissue==we.extreme_tissue)].iloc[0]; put("ex_extreme_naive_fold", float(2**abs(rx.logFC_naive)), rel(f"{RES}/bias/all_tissues.primary.tsv")); put("ex_extreme_fvb_fold", float(2**abs(rx.logFC_fvb)), rel(f"{RES}/bias/all_tissues.primary.tsv"))
# per-sample CPM for the extreme gene in its tissue (Fig. 6B)
f6=open(f"{FIG}/Fig6_worked_examples_data.tsv").read().split("## panelB")[1].strip().split("\n"); pb=pd.DataFrame([l.split("\t") for l in f6[1:]],columns=f6[0].split("\t")); pb["v"]=pb.log2_cpm_plus1.astype(float)
for arm in ("naive","wasp","fvb"):
    for s in ("FVB","B6"): put(f"ex_extreme_cpm_log2_{arm}_{s}_mean", float(pb[(pb.arm==arm)&(pb.strain==s)].v.mean()), rel(f"{FIG}/Fig6_worked_examples_data.tsv"))
# ---- within-FVB (reviewer check) ---------------------------------------------------------------------------------------
wf=pd.read_csv(f"{RES}/within_fvb/summary.tsv", sep="\t"); ev=rel(f"{RES}/within_fvb/summary.tsv")
for _,r in wf.iterrows():
    for c in ("DE_FDR05_naive","DE_FDR05_wasp","DE_FDR05_fvb","DE_concordant_FDR05"): put(f"wf_{r.tissue}_{c}", int(r[c]), ev)
    put(f"wf_{r.tissue}_pearson_r", float(r.pearson_r_logFC_naive_vs_fvb), ev); put(f"wf_{r.tissue}_pct_absdelta_ge_0.5", float(r["frac_abs_delta_ge_0.5"]*100), ev); put(f"wf_{r.tissue}_catalog_median_absdelta", float(r.catalog_genes_median_abs_delta), ev)
# ---- negative control -----------------------------------------------------------------------------------------------
rg=pd.read_csv(f"{RES}/negative_control/reference_swap_by_group.tsv", sep="\t"); ev=rel(f"{RES}/negative_control/reference_swap_by_group.tsv")
for _,r in rg.iterrows():
    k=r.group.replace(" ","_").replace("/","").replace("(","").replace(")","").replace("-","_")
    put(f"nc_{k}_n", int(r.n_animals), ev); put(f"nc_{k}_catalog_median_mean", float(r.catalog_genes_median_swap_mean), ev); put(f"nc_{k}_catalog_median_sd", float(r.catalog_genes_median_swap_sd), ev); put(f"nc_{k}_catalog_median_min", float(r.catalog_median_swap_min), ev); put(f"nc_{k}_catalog_median_max", float(r.catalog_median_swap_max), ev)
    put(f"nc_{k}_gw_median_mean", float(r["median_mean"]), ev); put(f"nc_{k}_gw_pct_ge_0.5_mean", float(r["frac_ge_0.5_mean"]*100), ev); put(f"nc_{k}_gw_pct_le_-0.5_mean", float(r["frac_le_-0.5_mean"]*100), ev)
rs=pd.read_csv(f"{RES}/negative_control/reference_swap_summary.tsv", sep="\t"); evr=rel(f"{RES}/negative_control/reference_swap_summary.tsv"); py=rs[rs.tissue=="mammary_tumour_PyVT"]
put("nc_pyvt_n_genes_min", int(py.n_genes.min()), evr); put("nc_pyvt_catalog_genes_n", int(py.catalog_genes_n.iloc[0]), evr); put("nc_pyvt_gw_median_absmax", float(py["median"].abs().max()), evr); put("nc_pyvt_gw_pct_ge_0.5_max", float(py["frac_ge_0.5"].max()*100), evr)
ho=pd.read_csv(f"{RES}/negative_control/heldout_catalog_swap.by_group.tsv", sep="\t"); ev=rel(f"{RES}/negative_control/heldout_catalog_swap.by_group.tsv")
for _,r in ho.iterrows():
    k=r.group.replace(" ","_").replace("/","").replace("(","").replace(")","").replace("-","_"); put(f"ho_{r.held_out_tissue}_{k}_median_mean", float(r.median_swap_mean), ev); put(f"ho_{r.held_out_tissue}_{k}_sd", float(r.median_swap_sd), ev); put(f"ho_{r.held_out_tissue}_n_catalog", int(r.n_catalog), ev)
b6h=ho[ho.group.str.startswith("C57BL/6J")]; pyh=ho[ho.group=="PyVT (SE)"]; fvh=ho[ho.group.str.startswith("FVB/NJ")]
put("ho_B6_median_min", float(b6h.median_swap_mean.min()), ev); put("ho_B6_median_max", float(b6h.median_swap_mean.max()), ev); put("ho_PyVT_median_min", float(pyh.median_swap_mean.min()), ev); put("ho_PyVT_median_max", float(pyh.median_swap_mean.max()), ev); put("ho_FVB_median_min", float(fvh.median_swap_mean.min()), ev); put("ho_FVB_median_max", float(fvh.median_swap_mean.max()), ev)
fs=pd.read_csv(f"{RES}/negative_control/pyvt_fullsite_chr3_10_16.by_group.tsv", sep="\t"); ev=rel(f"{RES}/negative_control/pyvt_fullsite_chr3_10_16.by_group.tsv")
for _,r in fs.iterrows():
    k=r.group.replace(" ","_").replace("/","").replace("(","").replace(")","").replace("-","_"); put(f"fs_{k}_chr{r.chrom}_mean_af", float(r.mean_af), ev); put(f"fs_{k}_chr{r.chrom}_max_af", float(r.max_af), ev); put(f"fs_{k}_chr{r.chrom}_mean_n_sites", float(r.mean_n_sites), ev); put(f"fs_{k}_chr{r.chrom}_max_run", int(r.max_longest_run), ev)
ps=pd.read_csv(f"{RES}/negative_control/pyvt_fullsite_chr3_10_16.per_site.tsv", sep="\t"); ev=rel(f"{RES}/negative_control/pyvt_fullsite_chr3_10_16.per_site.tsv"); c16=ps[ps.contig==16]
put("fs_chr16_block_start", int(c16.position.min()), ev); put("fs_chr16_block_end", int(c16.position.max()), ev); put("fs_chr16_block_n_sites", int(c16.position.nunique()), ev); put("fs_chr16_block_animals", int(c16.sample_id.nunique()), ev)
put("fs_chr16_block_site_mean_af_min", float(c16.groupby("position").af.mean().min()), ev); put("fs_chr16_block_site_mean_af_max", float(c16.groupby("position").af.mean().max()), ev); put("fs_chr16_block_shared_with_B6", int(c16["also_af_ge_0.1_in_B6_comparator"].sum()), ev)
put("fs_pyvt_sites_af_ge_0.1_distinct", int(ps[["contig","position"]].drop_duplicates().shape[0]), ev); put("fs_pyvt_sites_af_ge_0.1_shared_with_B6", int(ps.drop_duplicates(["contig","position"])["also_af_ge_0.1_in_B6_comparator"].sum()), ev)
pa=pd.read_csv(f"{RES}/negative_control/pyvt_fullsite_chr3_10_16.per_animal.tsv", sep="\t"); pv=pa[(pa.group=="PyVT (SE)")&(pa.chrom==16)]; put("fs_pyvt_chr16_mean_af_min", float(pv.mean_af.min()), ev); put("fs_pyvt_chr16_mean_af_max", float(pv.mean_af.max()), ev); put("fs_pyvt_chr16_frac_mid_min", float(pv["frac_af_0.1_0.9"].min()), ev); put("fs_pyvt_chr16_frac_mid_max", float(pv["frac_af_0.1_0.9"].max()), ev)
st=pd.read_csv(f"{RES}/negative_control/catalog_swap_by_chrom_stratum.tsv", sep="\t"); ev=rel(f"{RES}/negative_control/catalog_swap_by_chrom_stratum.tsv")
gm=meta.chrom.reindex(cat.index); blk=cat.index[(gm=="16")&(meta.start.reindex(cat.index)>=37_000_000)&(meta.end.reindex(cat.index)<=50_000_000)]
put("catalog_genes_in_chr16_block", ",".join(cat.loc[blk,"gene_name"]), evc); put("catalog_genes_in_chr16_block_n", len(blk), evc)
pyg=pd.read_csv(f"{RES}/negative_control/pyvt_reference_swap_per_gene.tsv", sep="\t", index_col=0)
for g in blk:
    if g in pyg.index: put(f"nc_pyvt_median_swap_{cat.loc[g,'gene_name']}", float(pyg.loc[g].median()), rel(f"{RES}/negative_control/pyvt_reference_swap_per_gene.tsv"))
    else: put(f"nc_pyvt_median_swap_{cat.loc[g,'gene_name']}", "not expressed in PyVT", rel(f"{RES}/negative_control/pyvt_reference_swap_per_gene.tsv"))
excl=pd.read_csv(f"{RES}/negative_control/reference_swap_summary.tsv", sep="\t")
put("nc_pyvt_catalog_median_excl_chr16block_mean", float(np.mean([pyg.loc[[g for g in cat.index if g in pyg.index and g not in set(blk)],s].median() for s in pyg.columns])), rel(f"{RES}/negative_control/pyvt_reference_swap_per_gene.tsv"))
# ---- matched-reference sensitivity -----------------------------------------------------------------------------------
sm=pd.read_csv(f"{RES}/sensitivity_matched/summary.tsv", sep="\t"); ev=rel(f"{RES}/sensitivity_matched/summary.tsv")
for _,r in sm.iterrows(): put("matched_"+r["item"].replace(" ","_").replace("(","").replace(")","").replace(":","").replace("/","_").replace("=",""), r["value"], ev)
# ---- gene universe / references -------------------------------------------------------------------------------------
put("gtf_genes_total", len(meta), "analysis_v2/results/gene_meta.tsv"); put("gtf_genes_in_fvb", int((meta.in_fvb_gtf=="Y").sum()), "gene_meta.tsv"); put("gtf_genes_not_in_fvb", int((meta.in_fvb_gtf=="N").sum()), "gene_meta.tsv"); put("gtf_genes_not_in_fvb_pct", float((meta.in_fvb_gtf=="N").mean()*100), "gene_meta.tsv")
put("gtf_genes_not_in_fvb_zero_snp_pct", float((meta.n_snps[meta.in_fvb_gtf=="N"]==0).mean()*100), "gene_meta.tsv")
qcs=open(f"{PROJ}/analysis_v2/qc/qc_sites_subsample.md5").read(); put("qc_sites_full_n", 203510, "analysis_v2/qc/fvbnj_homalt_exonic.hetGT.vcf.gz (bcftools index -n)"); put("qc_sites_subset_n", 20351, "analysis_v2/qc/fvbnj_homalt_exonic.hetGT.sub10pct.seed20261016.vcf.gz")
put("junctions_fvb_index", 286280, "revision_v2/audit/junction_sets_fvb_index_vs_b6_gtf.txt"); put("junctions_b6_gtf", 290071, "same"); put("junctions_b6_only", 3791, "same"); put("junctions_b6_only_chrY", 3688, "same"); put("junctions_total_fvb_arm", 576397, "same")
# FVB/NJ variant set sizes (bcftools index -n; revision_v2/audit/fvb_vcf_counts.txt, 2026-10-20)
for k,v in (("vcf_fvbnj_snps_pass_minac1",5960448),("vcf_fvbnj_biallelic_snps",5778575),("vcf_fvbnj_indels",1214530),("vcf_fvbnj_homalt_autosomal_snps",5030214),("vcf_fvbnj_homalt_exonic_snps",203510)): put(k, v, "revision_v2/audit/fvb_vcf_counts.txt")
# Phase J: ranges of WASP attrition by strain (for the v1 sentences kept in the reconstruction) and a few catalog genes named in the text
_b6=[float(v["value"]) for v in N if v["key"].startswith("wasp_removed_") and v["key"].endswith("_B6_mean_pct") and "PyVT" not in v["key"]]
_fv=[float(v["value"]) for v in N if v["key"].startswith("wasp_removed_") and v["key"].endswith("_FVB_mean_pct")]
_pp=[float(v["value"]) for v in N if v["key"].startswith("wasp_removed_") and v["key"].endswith("_FVB_minus_B6_pp")]
ev=rel(f"{RES}/qc/mapping_stats_per_sample.tsv")
put("wasp_removed_B6_min_pct", min(_b6), ev); put("wasp_removed_B6_max_pct", max(_b6), ev); put("wasp_removed_FVB_min_pct", min(_fv), ev); put("wasp_removed_FVB_max_pct", max(_fv), ev)
put("wasp_removed_FVB_minus_B6_pp_min", min(_pp), ev); put("wasp_removed_FVB_minus_B6_pp_max", max(_pp), ev)
_cat=pd.read_csv(f"{RES}/catalog/catalog.primary.tsv", sep="\t").set_index("gene_name"); ev=rel(f"{RES}/catalog/catalog.primary.tsv")
for g in ("Gm25939","Ndufs5","Rps2","Rpl29","Klra1","Mup17","3000002C10Rik"):
    if g in _cat.index:
        r=_cat.loc[g]; put(f"catalog_gene_{g}_residual", float(r["residual_fraction"]), ev); put(f"catalog_gene_{g}_n_tissues", int(r["n_tissues"]), ev); put(f"catalog_gene_{g}_mean_bias", float(r["mean_bias"]), ev)
        for t in tissues:
            if pd.notna(r[t]): put(f"catalog_gene_{g}_{t}_bias", float(r[t]), ev)
# vW-tagged (SNP-overlapping) read fractions by group (revision_v2/audit/vw_tagged_fraction_by_group.txt; Phase J)
ev="revision_v2/audit/vw_tagged_fraction_by_group.txt"; _vw={}
for l in open(f"{PROJ}/{ev}"):
    m=re.match(r"^(.+?)\s{2,}(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$", l.rstrip())
    if m: _vw[m.group(1).strip()]=float(m.group(3))
if _vw:
    put("vw_overlap_pct_pyvt", _vw.get("PyVT (SE 101)", float("nan")), ev); _o=[v for k,v in _vw.items() if not k.startswith("PyVT")]; put("vw_overlap_pct_other_min", min(_o), ev); put("vw_overlap_pct_other_max", max(_o), ev)
put("junctions_inserted_at_mapping", int(float(dict((r["key"],r["value"]) for r in N)["junctions_total_fvb_arm"])-float(dict((r["key"],r["value"]) for r in N)["junctions_fvb_index"])), "revision_v2/audit/junction_sets_fvb_index_vs_b6_gtf.txt")
rm=pd.read_csv(f"{RES}/worked_examples/rps2_multimap_cpm.tsv", sep="\t"); ev=rel(f"{RES}/worked_examples/rps2_multimap_cpm.tsv")   # Phase K-1(a)
for _,r in rm.iterrows():
    k=f"rps2mm_{r.gene.replace('-','_')}_{r.arm}_{'unique' if r.counting=='unique-only' else 'multi'}"; put(f"{k}_FVB_cpm", float(r.FVBNJ_mean_cpm), ev); put(f"{k}_B6_cpm", float(r.C57BL6J_mean_cpm), ev)
for arm in ("naive","fvb"):
    u=float(rm[(rm.gene=="Rps2")&(rm.arm==arm)&(rm.counting=="unique-only")].FVBNJ_mean_cpm.iloc[0]); m=float(rm[(rm.gene=="Rps2")&(rm.arm==arm)&(rm.counting!="unique-only")].FVBNJ_mean_cpm.iloc[0])
    put(f"rps2_multi_gain_cpm_{arm}", m-u, ev); put(f"rps2_multi_fold_{arm}", m/u, ev)
put("rps2_multi_gain_cpm", float(rm[(rm.gene=="Rps2")&(rm.arm=="naive")&(rm.counting!="unique-only")].FVBNJ_mean_cpm.iloc[0])-float(rm[(rm.gene=="Rps2")&(rm.arm=="naive")&(rm.counting=="unique-only")].FVBNJ_mean_cpm.iloc[0]), ev)
# Phase K-1'(a): every Rps2 / pseudogene CPM in the article uses ONE definition — shared-denominator mean CPM (each arm and counting mode
# divided by the naive-arm unique-only library size of the same animal; mean over 8 FVB/NJ or 8 C57BL/6J skin animals; the normalisation of
# Fig. 6B) of Rps2 and the three processed pseudogenes Rps2-ps10, -ps13, -ps6 (analysis_v2/results/worked_examples/rps2_multimap_cpm.tsv).
# The earlier keys summed all 12 annotated Rps2 pseudogenes (revision_v2/audit/rps2_pseudogene_counts_skin.tsv) and are no longer used.
rm3=pd.read_csv(f"{RES}/worked_examples/rps2_multimap_cpm.tsv", sep="\t"); ev=rel(f"{RES}/worked_examples/rps2_multimap_cpm.tsv")
for arm in ("naive","fvb"):
    for mode,lab in (("unique","unique-only"),("multi","-M --fraction")):
        a=rm3[(rm3.arm==arm)&(rm3.counting==lab)]; r=a[a.gene=="Rps2"].iloc[0]; p=a[a.gene!="Rps2"]; k=f"{arm}_{mode}"
        put(f"rps2_skin_{k}_FVB_cpm", float(r.FVBNJ_mean_cpm), ev); put(f"rps2_skin_{k}_B6_cpm", float(r.C57BL6J_mean_cpm), ev)
        put(f"rps2_pseudogenes_skin_{k}_FVB_cpm", float(p.FVBNJ_mean_cpm.sum()), ev); put(f"rps2_pseudogenes_skin_{k}_B6_cpm", float(p.C57BL6J_mean_cpm.sum()), ev)
        put(f"rps2_family_skin_{k}_FVB_cpm", float(a.FVBNJ_mean_cpm.sum()), ev)
    u=rm3[(rm3.arm==arm)&(rm3.gene=="Rps2")]; uu=u[u.counting=="unique-only"].iloc[0]; mm=u[u.counting!="unique-only"].iloc[0]
    put(f"rps2_skin_{arm}_log2ratio_unique", float(np.log2(uu.FVBNJ_mean_cpm/uu.C57BL6J_mean_cpm)), ev); put(f"rps2_skin_{arm}_log2ratio_multi", float(np.log2(mm.FVBNJ_mean_cpm/mm.C57BL6J_mean_cpm)), ev)
    put(f"rps2_skin_{arm}_log2ratio_shift", float(np.log2(mm.FVBNJ_mean_cpm/mm.C57BL6J_mean_cpm)-np.log2(uu.FVBNJ_mean_cpm/uu.C57BL6J_mean_cpm)), ev)
fm=rm3[(rm3.arm=="fvb")&(rm3.counting!="unique-only")]; put("rps2_skin_fvb_multi_share_of_family", float(fm[fm.gene=="Rps2"].FVBNJ_mean_cpm.iloc[0]/fm.FVBNJ_mean_cpm.sum()), ev)
ev="numbers.tsv: ex_extreme_*_residual"; _res=[float(v["value"]) for v in N if v["key"] in ("ex_extreme_adipose_residual","ex_extreme_hypothalamus_residual","ex_extreme_liver_residual","ex_extreme_skin_residual")]
if _res: put("ex_extreme_residual_min", min(_res), ev); put("ex_extreme_residual_max", max(_res), ev)
al=pd.read_csv(f"{RES}/worked_examples/rps2_pseudogene_alleles.tsv", sep="\t"); ev=rel(f"{RES}/worked_examples/rps2_pseudogene_alleles.tsv")   # K-1'(c)
sn=al[al.variant.str.startswith("SNP")]
put("rps2_exonic_snps_in_pseudogene_sequence", int((sn[sn.pseudogene=="Rps2-ps10"].allele_carried!="not in pseudogene").sum()), ev)
for g in ("Rps2-ps10","Rps2-ps13","Rps2-ps6"): put(f"rps2_{g.replace('-','_')}_fvb_allele_snps", int((sn[sn.pseudogene==g].allele_carried=="FVB/NJ").sum()), ev)
put("rps2_exonic_indels_in_pseudogene_sequence", int((al[(al.variant.str.startswith("indel"))&(al.pseudogene=="Rps2-ps10")].allele_carried!="not in pseudogene").sum()), ev)
put("rps2_exonic_snps_outside_pseudogene_sequence", int((sn[sn.pseudogene=="Rps2-ps10"].allele_carried=="not in pseudogene").sum()), ev)   # Phase L: exonic SNPs in exon segments the pseudogenes do not contain
pd.DataFrame(N).to_csv(OUT, sep="\t", index=False); print(f"{len(N)} numbers -> {OUT}")
