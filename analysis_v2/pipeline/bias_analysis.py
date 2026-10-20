#!/usr/bin/env python3
"""analysis_v2 downstream analysis (Methods fixed in analysis_v2/params.yaml; rule 0.4-8).
Inputs : analysis_v2/results/counts/<tissue>.{naive,wasp,fvb}.{unique,multi}.tsv + <tissue>.samples.tsv (merge_counts.py),
         analysis_v2/results/gene_meta.tsv (gene_meta.sh), analysis_v2/pipeline/run_edger.R
Outputs: analysis_v2/results/
  edger/<tissue>.<design>.<arm>.<mode>.tsv      edgeR tables (logFC = FVB vs B6 strain coefficient)
  bias/<tissue>.<design>.tsv                    per-gene logFC per arm, bias (fvb - naive), bias_wasp (wasp - naive), bias_multi, delta
  bias/all_tissues.<design>.tsv                 long table over tissues
  catalog/catalog.<design>.tsv                  bias >= 0.5 in >= 2 tissues AND in FVB GTF (+ multimap mechanism class)
  catalog/catalog_sensitivity.<design>.tsv      thresholds 0.25/0.5/1.0 x min tissues 1..N
  predictors/gene_predictors.<design>.tsv, predictors/summary.<design>.tsv
  negative_control/pyvt_reference_swap.tsv, pyvt_per_chromosome_af.tsv
  within_fvb/<tissue>.diet.<arm>.tsv, within_fvb/summary.tsv  (reviewer 2 comment 1)
Design keys: primary = "~ diet + strain" (tissues with two diets) or "~ strain"; water = control-diet animals only, "~ strain".
"""
import os, sys, subprocess, csv, itertools, collections
import numpy as np, pandas as pd
from scipy import stats
PROJ="/usr/local/jupyter/FVB_bias"; FINAL=f"{PROJ}/analysis_v2/results"
RES=os.environ.get("RES_DIR", FINAL)   # outputs (and counts from merge_counts.py); gene_meta/f1_ase inputs always from the final results dir
CNT=f"{RES}/counts"; R=f"{PROJ}/analysis_v2/pipeline/run_edger.R"
for d in ("edger","bias","catalog","predictors","negative_control","within_fvb"): os.makedirs(f"{RES}/{d}", exist_ok=True)
BIAS_THR=0.5; MIN_TISSUES=2; DELTA_THR=0.3
def _two_strains(t):
    """partial dry runs only: a tissue enters the strain contrast once both strains have >= 2 animals (final run: always true)"""
    n=pd.read_csv(f"{CNT}/{t}.samples.tsv", sep="\t")["strain"].value_counts()
    if len(n)<2 or n.min()<2: print(f"[skip] {t}: strains present {dict(n)} (needs both strains, >= 2 each)"); return False
    return True
STRAIN_TISSUES=[t for t in ("adipose","hypothalamus","liver","skin") if os.path.exists(f"{CNT}/{t}.samples.tsv") and _two_strains(t)]
meta=pd.read_csv(f"{FINAL}/gene_meta.tsv", sep="\t", index_col=0)

def edger(tissue, design, arm, mode, formula, coef, filt=""):
    out=f"{RES}/edger/{tissue}.{design}.{arm}.{mode}.tsv"
    if not os.path.exists(out):
        r=subprocess.run(["Rscript",R,f"{CNT}/{tissue}.{arm}.{mode}.tsv",f"{CNT}/{tissue}.samples.tsv",formula,coef,out,filt],capture_output=True,text=True)
        if r.returncode!=0: sys.exit(f"edgeR failed {tissue} {design} {arm} {mode}: {r.stderr[-2000:]}")
        print(r.stdout.strip())
    return pd.read_csv(out, sep="\t", index_col=0)

CONTROL_LEVELS=("water","control")
def control_level(ss): return next((c for c in CONTROL_LEVELS if c in set(ss["diet"])), None)
def designs_for(tissue):
    """primary: ~ [diet +] [sex +] strain on all animals; water: control-diet/treatment animals only (~ [sex +] strain)."""
    ss=pd.read_csv(f"{CNT}/{tissue}.samples.tsv", sep="\t")
    sexterm=" sex +" if ss["sex"].fillna("").nunique()>1 else ""
    d={}
    if ss["diet"].nunique()>1:
        d["primary"]=(f"~ diet +{sexterm} strain","strainFVB","")
        d["water"]=(f"~{sexterm} strain","strainFVB",f'diet == "{control_level(ss)}"')
    else:
        d["primary"]=(f"~{sexterm} strain","strainFVB","")
    return d

def bias_table(tissue, design):
    formula,coef,filt=designs_for(tissue)[design]
    nv=edger(tissue,design,"naive","unique",formula,coef,filt); wp=edger(tissue,design,"wasp","unique",formula,coef,filt); fv=edger(tissue,design,"fvb","unique",formula,coef,filt)
    nvm=edger(tissue,design,"naive","multi",formula,coef,filt); fvm=edger(tissue,design,"fvb","multi",formula,coef,filt)
    g=nv.index.intersection(fv.index)
    t=pd.DataFrame({"tissue":tissue,"logFC_naive":nv.loc[g,"logFC"],"FDR_naive":nv.loc[g,"FDR"],"logCPM_naive":nv.loc[g,"logCPM"],
                    "logFC_wasp":wp["logFC"].reindex(g),"FDR_wasp":wp["FDR"].reindex(g),"logFC_fvb":fv.loc[g,"logFC"],"FDR_fvb":fv.loc[g,"FDR"]})
    t["bias"]=t["logFC_fvb"]-t["logFC_naive"]; t["bias_wasp"]=t["logFC_wasp"]-t["logFC_naive"]
    gm=nvm.index.intersection(fvm.index); bm=(fvm.loc[gm,"logFC"]-nvm.loc[gm,"logFC"]); t["bias_multi"]=bm.reindex(g)
    t["delta_multi_minus_unique"]=t["bias_multi"]-t["bias"]; t["residual_fraction"]=t["bias_multi"]/t["bias"]
    t.index.name="gene_id"; t.to_csv(f"{RES}/bias/{tissue}.{design}.tsv", sep="\t")
    return t

def run_design(design):
    tabs={}
    for tissue in STRAIN_TISSUES:
        if design not in designs_for(tissue): continue
        tabs[tissue]=bias_table(tissue,design)
    if not tabs: return
    long=pd.concat([v.assign(gene_id=v.index) for v in tabs.values()], ignore_index=True)
    long.to_csv(f"{RES}/bias/all_tissues.{design}.tsv", sep="\t", index=False)
    wide=long.pivot(index="gene_id", columns="tissue", values="bias"); tissues=list(wide.columns)
    # catalog
    hit=(wide>=BIAS_THR); n_hit=hit.sum(axis=1); in_fvb=meta["in_fvb_gtf"].reindex(wide.index).eq("Y")
    cat=wide[(n_hit>=MIN_TISSUES)&in_fvb].copy(); cat["n_tissues"]=n_hit[cat.index]; cat["n_tissues_tested"]=wide.loc[cat.index].notna().sum(axis=1)
    cat["mean_bias"]=wide.loc[cat.index].mean(axis=1); cat["tissues"]=hit.loc[cat.index].apply(lambda r:",".join(c for c in tissues if r[c]),axis=1)
    cat["gene_name"]=meta["gene_name"].reindex(cat.index); cat["biotype"]=meta["biotype"].reindex(cat.index); cat["snp_density_per_kb"]=meta["snp_density_per_kb"].reindex(cat.index)
    bm=long.pivot(index="gene_id",columns="tissue",values="bias_multi"); cat["mean_bias_multi"]=bm.loc[cat.index].mean(axis=1)
    cat["delta_multi_minus_unique"]=cat["mean_bias_multi"]-cat["mean_bias"]; cat["residual_fraction"]=cat["mean_bias_multi"]/cat["mean_bias"]
    cat["mechanism_class"]=np.where(cat["delta_multi_minus_unique"].abs()>DELTA_THR,"SNP-induced multi-mapping","direct SNP misalignment")
    cat=cat.sort_values("mean_bias",ascending=False); cat.index.name="gene_id"
    cat[["gene_name","biotype","n_tissues","n_tissues_tested","mean_bias","tissues"]+tissues+["snp_density_per_kb","mean_bias_multi","delta_multi_minus_unique","residual_fraction","mechanism_class"]].to_csv(f"{RES}/catalog/catalog.{design}.tsv", sep="\t")
    print(f"[{design}] tissues={tissues} genes_tested={len(wide)} catalog={len(cat)} (bias>={BIAS_THR} in >= {MIN_TISSUES} tissues, in FVB GTF)")
    # sensitivity grid
    rows=[]
    for thr in (0.25,0.5,1.0):
        for k in range(1,len(tissues)+1):
            n=((wide>=thr).sum(axis=1)>=k); rows.append({"threshold":thr,"min_tissues":k,"n_genes_all":int(n.sum()),"n_genes_in_fvb_gtf":int((n&in_fvb).sum())})
    pd.DataFrame(rows).to_csv(f"{RES}/catalog/catalog_sensitivity.{design}.tsv", sep="\t", index=False)
    # predictors: per-gene mean bias across the tissues in which the gene was tested (every gene tested in >= 1 tissue), as in the
    # old Methods ("per-gene mean bias across the datasets", all expressed genes). 2026-10-19 22:36: an earlier draft restricted
    # this to genes tested in all tissues (14,316 of 27,732) without a recorded decision; reverted to the old definition (rule 0.4-8).
    full=wide; mb=full.mean(axis=1); lc=long.pivot(index="gene_id",columns="tissue",values="logCPM_naive").loc[full.index].mean(axis=1)
    pr=pd.DataFrame({"mean_bias":mb,"mean_logCPM_naive":lc}).join(meta[["gene_name","biotype","length_bp","snp_density_per_kb"]])
    pr["log10_length"]=np.log10(pr["length_bp"]); pr["biotype_group"]=pr["biotype"].map(lambda b: b if b in ("protein_coding","lncRNA","miRNA") else ("pseudogene" if "pseudogene" in str(b) else "other"))
    pr.index.name="gene_id"; pr.to_csv(f"{RES}/predictors/gene_predictors.{design}.tsv", sep="\t")
    summ=[]
    for c in ("snp_density_per_kb","log10_length","mean_logCPM_naive"):
        rho,p=stats.spearmanr(pr[c],pr["mean_bias"],nan_policy="omit"); summ.append({"predictor":c,"n_genes":int(pr[c].notna().sum()),"spearman_rho":rho,"p_value":p})
    for b,gdf in pr.groupby("biotype_group"): summ.append({"predictor":f"biotype={b}","n_genes":len(gdf),"spearman_rho":np.nan,"p_value":np.nan,"mean_bias":gdf["mean_bias"].mean(),"frac_bias_ge_0.5":(gdf["mean_bias"]>=0.5).mean()})
    pd.DataFrame(summ).to_csv(f"{RES}/predictors/summary.{design}.tsv", sep="\t", index=False)
    return cat

def reference_swap(tissue, label):
    """per-sample log2((fvb_count+0.5)/(naive_count+0.5)) at genes expressed (naive CPM>1) — no strain contrast needed."""
    nv=pd.read_csv(f"{CNT}/{tissue}.naive.unique.tsv", sep="\t", index_col=0); fv=pd.read_csv(f"{CNT}/{tissue}.fvb.unique.tsv", sep="\t", index_col=0)
    ss=pd.read_csv(f"{CNT}/{tissue}.samples.tsv", sep="\t").set_index("sample_id")
    cpm=nv/nv.sum()*1e6; keep=(cpm>1).mean(axis=1)>=0.5
    ratio=np.log2((fv.loc[keep]+0.5)/(nv.loc[keep]+0.5)); ratio.index.name="gene_id"
    ratio.to_csv(f"{RES}/negative_control/{label}_reference_swap_per_gene.tsv", sep="\t")
    rows=[]
    for s in ratio.columns:
        v=ratio[s]; rows.append({"sample_id":s,"tissue":tissue,"strain":ss.loc[s,"strain"],"n_genes":int(v.notna().sum()),"median":v.median(),"mean":v.mean(),"q05":v.quantile(.05),"q95":v.quantile(.95),"frac_ge_0.5":(v>=0.5).mean(),"frac_le_-0.5":(v<=-0.5).mean()})
    return pd.DataFrame(rows), ratio

def negative_control(cat_genes):
    out=[]; per_gene={}
    if os.path.exists(f"{CNT}/mammary_tumour_PyVT.samples.tsv"):
        df,ratio=reference_swap("mammary_tumour_PyVT","pyvt"); out.append(df); per_gene["pyvt"]=ratio
    for t in STRAIN_TISSUES:
        df,ratio=reference_swap(t,t); out.append(df); per_gene[t]=ratio
    if not out: return
    allr=pd.concat(out, ignore_index=True)
    # catalog-gene summary per sample
    for t,ratio in per_gene.items():
        g=[x for x in cat_genes if x in ratio.index]
        for s in ratio.columns:
            allr.loc[allr["sample_id"]==s,"catalog_genes_n"]=len(g); allr.loc[allr["sample_id"]==s,"catalog_genes_median_swap"]=ratio.loc[g,s].median() if g else np.nan
    allr.to_csv(f"{RES}/negative_control/reference_swap_summary.tsv", sep="\t", index=False)
    # per-chromosome FVB allele fraction for PyVT (congenic interval check)
    rows=[]
    for s in os.listdir(f"{PROJ}/analysis_v2/data/mammary_tumour_PyVT") if os.path.isdir(f"{PROJ}/analysis_v2/data/mammary_tumour_PyVT") else []:
        p=f"{PROJ}/analysis_v2/data/mammary_tumour_PyVT/{s}/qc/{s}.ase.tsv"
        if not os.path.exists(p): continue
        a=pd.read_csv(p, sep="\t"); a=a[(a.refCount+a.altCount)>=20]; a["af"]=a.altCount/(a.refCount+a.altCount)
        for c,g in a.groupby("contig"): rows.append({"sample_id":s,"chrom":c,"n_sites":len(g),"mean_af":g.af.mean(),"frac_homalt":(g.af>0.9).mean(),"frac_het":((g.af>=0.1)&(g.af<=0.9)).mean()})
    if rows: pd.DataFrame(rows).to_csv(f"{RES}/negative_control/pyvt_per_chromosome_af.tsv", sep="\t", index=False)

def within_fvb(cat_genes):
    rows=[]
    for t in STRAIN_TISSUES:
        ss=pd.read_csv(f"{CNT}/{t}.samples.tsv", sep="\t")
        if ss["diet"].nunique()<2: continue
        res={}; ctrl=control_level(ss); trt=[x for x in ss["diet"].unique() if x!=ctrl][0]
        sexterm=" sex +" if ss["sex"].fillna("").nunique()>1 else ""
        for arm in ("naive","wasp","fvb"):
            out=f"{RES}/within_fvb/{t}.diet.{arm}.tsv"
            if not os.path.exists(out):
                r=subprocess.run(["Rscript",R,f"{CNT}/{t}.{arm}.unique.tsv",f"{CNT}/{t}.samples.tsv",f"~{sexterm} diet",f"diet{trt}",out,'strain == "FVB"'],capture_output=True,text=True)
                if r.returncode!=0: sys.exit(f"within-FVB edgeR failed {t} {arm}: {r.stderr[-1500:]}")
                print(r.stdout.strip())
            res[arm]=pd.read_csv(out, sep="\t", index_col=0)
        g=res["naive"].index.intersection(res["fvb"].index); d=res["fvb"].loc[g,"logFC"]-res["naive"].loc[g,"logFC"]
        gc=[x for x in cat_genes if x in g]
        rows.append({"tissue":t,"n_genes":len(g),"pearson_r_logFC_naive_vs_fvb":np.corrcoef(res["naive"].loc[g,"logFC"],res["fvb"].loc[g,"logFC"])[0,1],
                     "DE_FDR05_naive":int((res["naive"].FDR<0.05).sum()),"DE_FDR05_wasp":int((res["wasp"].FDR<0.05).sum()),"DE_FDR05_fvb":int((res["fvb"].FDR<0.05).sum()),
                     "median_abs_delta_logFC":d.abs().median(),"frac_abs_delta_ge_0.5":(d.abs()>=0.5).mean(),"catalog_genes_n":len(gc),"catalog_genes_median_abs_delta":d.loc[gc].abs().median() if gc else np.nan,
                     "DE_concordant_FDR05":int(((res["naive"].loc[g].FDR<0.05)&(res["fvb"].loc[g].FDR<0.05)&(np.sign(res["naive"].loc[g].logFC)==np.sign(res["fvb"].loc[g].logFC))).sum())})
    if rows: pd.DataFrame(rows).to_csv(f"{RES}/within_fvb/summary.tsv", sep="\t", index=False)

if __name__=="__main__":
    cat=run_design("primary"); run_design("water")
    cat_genes=list(cat.index) if cat is not None else []
    negative_control(cat_genes); within_fvb(cat_genes)
    print("bias_analysis DONE")
