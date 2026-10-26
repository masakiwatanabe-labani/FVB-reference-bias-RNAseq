#!/usr/bin/env python3
"""Phase M-2 — independent recomputation of the numbers printed in the deliverables.

Independence: this script imports no module of the pipeline (analysis_v2/pipeline) or of the text/number generators
(revision_v2/scripts), and does not read manuscript/revised_v2/generated/numbers.tsv. It reads only primary result files:
per-tissue x arm edgeR tables, the matched-reference edgeR tables, per-animal featureCounts matrices and .summary files,
STAR splice-junction lists and Log.out, F1 per-site ASEReadCounter tables and the WASP allele-removal table, identity per-site
tables, within-strain DE tables, the MGP-derived FVB_NJ VCFs and the Ensembl 110 GTF (through revision_v2/freeze/cache/, written by
cache/vcf_counts.sh and cache/gtf_vcf_cache.sh with bcftools/bedtools only), rps2_pseudogene_alleles.tsv, the vW tag counts, the
PyVT full-site table, the survey table and (for the response letter only) the submitted-version catalog under legacy/.
Every quantity is computed here from those files with the definitions stated in the Methods, not taken from a derived table.

For every claim: claim id, group, location (file; the needle string must occur in that deliverable), text value, exact value,
value at text precision, PASS/FAIL (REWRITE: printed value superseded by a Phase W correction, see REWRITE). Bounds ("within ±x", "fewer than y%", "at most z", "≤", "≥") PASS only if the exact value
satisfies the bound as printed. Output: revision_v2/freeze/recompute_report.tsv. Usage: recompute_independent.py"""
import re, os, glob, gzip, math, csv, collections
import numpy as np, pandas as pd
from scipy import stats

P = "/usr/local/jupyter/FVB_bias"; RES = f"{P}/analysis_v2/results"; DATA = f"{P}/analysis_v2/data"; C = f"{P}/revision_v2/freeze/cache"
SUB = f"{P}/revision_v2/submission"
TISSUES = ["adipose", "hypothalamus", "liver", "skin"]

# ---------------------------------------------------------------- deliverable texts (needle check) ---------------------------
def norm(s): return re.sub(r"\s+", " ", s.replace("−", "-").replace(" ", " ").replace("\xa0", " "))
TXT = {"main": norm(open(f"{SUB}/LAR_main_text_v2.plain.txt", encoding="utf-8").read()),
       "af1": norm(open(f"{SUB}/LAR_additional_file_1.plain.txt", encoding="utf-8").read()),
       "resp": norm(open(f"{SUB}/response_to_reviewers_v2.plain.txt", encoding="utf-8").read()),
       "cover": norm(open(f"{SUB}/cover_letter_v2.plain.txt", encoding="utf-8").read())}
FIG = {}
for f in glob.glob(f"{SUB}/figures/*.pdf"):
    FIG[os.path.basename(f)[:-4]] = norm(os.popen(f"pdftotext -layout '{f}' -").read())
ROWS = []
def textval(s):
    s = s.replace("−", "-").replace(",", "").replace("+", "").replace("%", "").strip()
    return float(s)
def decimals(s):
    s = s.replace("%", "").strip(); return len(s.split(".")[1]) if "." in s else 0
# Phase W (2026-10-26): text claims whose printed value is superseded by a Phase W correction. The deliverable texts are the authors'
# (2026-10-24 dumps, not rebuilt here); the exact value is recomputed with the corrected definition and the claim is reported as
# REWRITE (printed value pending the authors' rewrite) instead of FAIL. A REWRITE claim whose text already matches becomes PASS.
REWRITE = {"K5": "W-1 survey re-screen (RNA-Seq library strategy, FVB background per sample)", "K6": "W-1 survey re-screen", "K7": "W-1 survey re-screen",
           "K8r": "W-1 survey re-screen", "D11": "W-2: |bias_multi|/|bias_unique| <= 0.30 (the signed ratio also admits six sign-reversing genes)",
           "W3a": "W-3: samples, not animals (fructose-study tissues share animals)", "W3b": "W-3: samples, not animals", "W3c": "W-3: 83 samples from the animals counted here"}
def claim(cid, group, where, needle, text, exact, kind="num", note="", scale=1.0):
    """kind: num (rounded equality at printed precision) | int | le (|exact| <= text) | lt (exact < text) | ge (exact >= text) |
    gt | sci (mantissa 1 dp + exponent) | plt (p < text) | nr (not recomputable) | eq (string)"""
    src = TXT.get(where) if where in TXT else FIG.get(where, "")
    found = (norm(needle) in src) if needle else True
    ex = exact
    if kind == "nr":
        status = "NOT_RECOMPUTABLE"; atp = ""
    elif kind == "eq":
        atp = str(ex); status = "PASS" if str(ex) == str(text) else "FAIL"
    elif kind == "sci":
        atp = f"{ex:.1e}"; m, e = text.split("e"); status = "PASS" if atp == f"{float(m):.1f}e{int(e):+03d}" else "FAIL"
    else:
        x = float(ex) * scale; t = textval(text); d = decimals(text)
        atp = f"{x:.{d}f}"
        if kind == "num": status = "PASS" if abs(x - t) <= 0.5 * 10 ** -d + 1e-9 else "FAIL"
        elif kind == "int": status = "PASS" if int(round(x)) == int(t) and abs(x - round(x)) < 1e-9 else "FAIL"
        elif kind == "le": status = "PASS" if abs(x) <= t + 1e-12 else "FAIL"
        elif kind == "lt": status = "PASS" if x < t else "FAIL"
        elif kind == "ge": status = "PASS" if x >= t - 1e-12 else "FAIL"
        elif kind == "gt": status = "PASS" if x > t else "FAIL"
        elif kind == "plt": status = "PASS" if x < t else "FAIL"
        else: raise ValueError(kind)
    if not found: status = "FAIL" if status != "NOT_RECOMPUTABLE" else status; note = (note + "; " if note else "") + "NEEDLE NOT FOUND in " + where
    if cid in REWRITE and found and status == "FAIL": status = "REWRITE"; note = (note + "; " if note else "") + REWRITE[cid]
    ROWS.append(dict(claim_id=cid, group=group, location=where, needle=needle[:120], text_value=text, exact=(f"{ex:.6g}" if isinstance(ex, (float, np.floating)) else str(ex)),
                     at_text_precision=atp, kind=kind, status=status, note=note))

# ---------------------------------------------------------------- annotation / variants --------------------------------------
G = pd.read_csv(f"{C}/genes_b6.snps_all.bed", sep="\t", header=None, names=["chrom", "start0", "end", "gene_id", "gene_name", "biotype", "strand", "snps_all"])
Gh = pd.read_csv(f"{C}/genes_b6.snps_homalt.bed", sep="\t", header=None, usecols=[3, 7], names=["gene_id", "snps_homalt"])
G = G.merge(Gh, on="gene_id"); G["length"] = G.end - G.start0; G["density"] = G.snps_all / (G.length / 1000)
G = G.set_index("gene_id")
FVBGTF = set(open(f"{C}/genes_fvbgtf.ids").read().split())
name2id = {}
for gid, nm in G.gene_name.items(): name2id.setdefault(nm, gid)
def gid(name): return name2id[name]

g = "A variant set and references"
n_snp = int(open(f"{C}/snps_homalt_auto.n").read()); n_ind_all = int(open(f"{C}/indels_all.n").read()); n_ind_ha = int(open(f"{C}/indels_homalt_auto.n").read())
n_exonic = int(open(f"{C}/snps_homalt_auto_exonic.n").read())
claim("A1", g, "main", "about 5.0 million homozygous autosomal SNPs", "5.0", n_snp / 1e6, note=f"{n_snp:,} PASS biallelic hom-alt autosomal SNVs (FVB_NJ.snps.GRCm39.vcf.gz)")
claim("A2", g, "main", "and 1.1 million homozygous autosomal indels", "1.1", n_ind_ha / 1e6, note=f"{n_ind_ha:,} homozygous-alternative autosomal indels (FVB_NJ.indels.GRCm39.vcf.gz); all records {n_ind_all:,}")
claim("A2r", g, "resp", "1.1 million homozygous autosomal indels", "1.1", n_ind_ha / 1e6)
claim("A3", g, "main", "203,510 of the SNPs within annotated exons", "203510", n_exonic, "int", note="hom-alt autosomal PASS biallelic SNVs overlapping an Ensembl 110 exon (bedtools)")
sub = set();
for l in gzip.open(f"{P}/analysis_v2/qc/fvbnj_homalt_exonic.hetGT.sub10pct.seed20261016.vcf.gz", "rt"):
    if l[0] != "#": c, p = l.split("\t")[:2]; sub.add((c, int(p)))
full = set()
for l in gzip.open(f"{P}/analysis_v2/qc/fvbnj_homalt_exonic.hetGT.vcf.gz", "rt"):
    if l[0] != "#": c, p = l.split("\t")[:2]; full.add((c, int(p)))
claim("A4", g, "main", "(20,351 sites; the site list and its validation are in the repository)", "20351", len(sub), "int", note=f"subset ⊂ full set: {sub <= full}; 10% of {len(full):,} = {len(full)*0.1:.1f}")
claim("A4b", g, "main", "fixed random 10% subset", "10", 100 * len(sub) / len(full), "num")
# subset validation: every sample evaluated at both site sets (full-set ase.tsv present), metrics over sites with depth >= 20
def idmet(d):
    d = d[d.totalCount >= 20]; af = d.altCount / d.totalCount
    return np.array([af.mean(), (af > 0.9).mean(), ((af >= 0.1) & (af <= 0.9)).mean(), (af < 0.1).mean()])
dif_all, dif_animals, verd = [], [], []
def vrule(m):
    af, hi, mid, lo = m
    if af > 0.90 and hi > 0.85: return "FVB/NJ pure"
    if af < 0.10 and lo > 0.85: return "B6 pure"
    if mid > 0.5 and 0.25 <= af <= 0.55: return "F1"
    if hi > 0.15 and mid > 0.15 and lo > 0.15: return "segregating"
    if hi > 0.25 and lo > 0.25 and mid < 0.15: return "other inbred"
    return "needs review"
for f in sorted(glob.glob(f"{DATA}/*/*/qc/*.ase.tsv")):
    s = os.path.basename(f)[:-8]; d = pd.read_csv(f, sep="\t", dtype={"contig": str})
    d2 = d[[(c, p) in sub for c, p in zip(d.contig, d.position)]]
    if len(d2) == len(d[d.totalCount >= 0]): continue   # sample evaluated at the subset only
    m1, m2 = idmet(d), idmet(d2); dd = np.abs(m1 - m2).max(); dif_all.append((s, dd)); verd.append(vrule(m1) == vrule(m2))
    if not s.startswith("heart"): dif_animals.append(dd)
all_same = all(verd); mx = max(x for _, x in dif_all); who = max(dif_all, key=lambda t: t[1])[0]
# Phase T (2026-10-24): the subset-validation sentence was deleted from Methods (the validation is in identity_check/ of the
# repository), so A5n/A5v/A5 are no longer text claims; the values above are still computed (dif_all, all_same, mx) but not reported.
# claim("A5n", g, "main", "identical verdicts for the 14 animals evaluated at both site sets", "14", len(dif_all), "int")   # retired in Phase T
# claim("A5v", g, "main", "identical verdicts for the 14 animals", "identical", "identical" if all_same else "DIFFER", "eq", note="verdict rules of Methods applied to both site sets")   # retired in Phase T
# claim("A5", g, "main", "(metrics within 0.014 and identical verdicts", "0.014", mx, "le", note=f"max |Δ| over {len(dif_all)} samples evaluated at both site sets = {mx:.4f} ({who}, frac_hi); over the 9 analysis animals only = {max(dif_animals):.4f}")   # retired in Phase T
b6ids = set(G.index); miss = sorted(b6ids - FVBGTF)
claim("A6", g, "main", "1,727 of 56,941 annotated genes were removed", "1727", len(miss), "int")
claim("A7", g, "main", "1,727 of 56,941 annotated genes were removed", "56941", len(b6ids), "int")
claim("A8", g, "af1", "1,727 (3.0%) are absent", "3.0", 100 * len(miss) / len(b6ids), "num")
claim("A9", g, "af1", "all of them have zero FVB/NJ SNPs in the gene body", "0", int(G.loc[miss, "snps_all"].max()), "int", note="max SNP count in gene body among the 1,727 (all FVB_NJ SNP records)")
mb = G.loc[miss, "biotype"].value_counts()
claim("A9b", g, "main", "most excluded genes were pseudogenes or lncRNAs", "0.5", (G.loc[miss, "biotype"].str.contains("pseudogene") | (G.loc[miss, "biotype"] == "lncRNA")).mean(), "gt", note="fraction pseudogene or lncRNA among the 1,727 (> 0.5 = 'most'); top biotypes " + ", ".join(f"{k} {v}" for k, v in mb.head(4).items()))
# junctions
def sj(path): return set(tuple(l.split("\t")[:4]) for l in open(path))
sjb, sjf = sj(f"{P}/ref/star_b6/sjdbList.out.tab"), sj(f"{P}/ref/star_fvb/sjdbList.out.tab")
claim("A10", g, "af1", "(286,280 annotated junctions)", "286280", len(sjf), "int")
claim("A11", g, "af1", "Ensembl 110 B6 GTF (290,071 junctions)", "290071", len(sjb), "int")
logout = open(glob.glob(f"{DATA}/hypothalamus/hypothalamus_B6_CWHyp4/bam/*.fvb.Log.out")[0]).read()
claim("A12", g, "af1", "290,117 collapsed junctions", "290117", int(re.search(r"(\d+) collapsed junctions", logout).group(1)), "int", note="FVB-arm STAR Log.out (hypothalamus_B6_CWHyp4)")
claim("A13", g, "af1", "576,397 junctions in total", "576397", int(re.search(r"Total junctions: (\d+)", logout).group(1)), "int")
bn = sjb - sjf
claim("A14", g, "af1", "3,791 junctions of the B6 GTF are absent", "3791", len(bn), "int", note=f"FVB-not-in-B6 = {len(sjf - sjb)}")
claim("A15", g, "af1", "3,688 of them on the Y chromosome", "3688", sum(1 for j in bn if j[0] == "Y"), "int")
claim("A15b", g, "main", "3,791 junctions absent from the converted GTF, 3,688 of them", "3791", len(bn), "int")

# ---------------------------------------------------------------- edgeR: tested genes, bias ---------------------------------
def eg(t, design, arm, cnt): return pd.read_csv(f"{RES}/edger/{t}.{design}.{arm}.{cnt}.tsv", sep="\t").set_index("gene_id")
B = {}
for t in TISSUES:
    n, f, w = eg(t, "primary", "naive", "unique"), eg(t, "primary", "fvb", "unique"), eg(t, "primary", "wasp", "unique")
    nm, fm = eg(t, "primary", "naive", "multi"), eg(t, "primary", "fvb", "multi")
    mt = pd.read_csv(f"{RES}/sensitivity_matched/edger/{t}.primary.matched.unique.tsv", sep="\t").set_index("gene_id")
    d = pd.DataFrame({"lfc_naive": n.logFC, "lfc_fvb": f.logFC}).dropna()
    d["logCPM_naive"] = n.logCPM.reindex(d.index); d["lfc_wasp"] = w.logFC.reindex(d.index)
    d["bias"] = d.lfc_fvb - d.lfc_naive; d["bias_wasp"] = d.lfc_wasp - d.lfc_naive
    d["bias_multi"] = (fm.logFC - nm.logFC).reindex(d.index)
    d["lfc_matched"] = mt.logFC.reindex(d.index); d["bias_matched"] = d.lfc_matched - d.lfc_naive; d["b6_loss"] = d.lfc_fvb - d.lfc_matched
    d["tissue"] = t; B[t] = d
ALL = pd.concat(B.values()); ALL.index.name = "gene_id"; ALL = ALL.reset_index()
genes_any = sorted(set(ALL.gene_id)); per_gene = ALL.groupby("gene_id")
ntest = per_gene.size()
g = "B genome-wide"
claim("B1", g, "main", "27,732 genes were expressed and quantifiable", "27732", len(genes_any), "int")
claim("B2", g, "main", "stronger among the 14,316 genes tested in all four tissues", "14316", int((ntest == 4).sum()), "int")
for t in TISSUES: claim(f"B3_{t}", g, "Fig3_strain_comparison" if t == "liver" else "main", "19,237" if t == "liver" else "", str(len(B[t])), len(B[t]), "int", note=f"tested genes {t} (naive ∩ FVB edgeR tables)")
med = {t: B[t].bias.median() for t in TISSUES}; mx = max(abs(v) for v in med.values())
claim("B4", g, "main", "median bias was within ±0.007 log₂ units of zero in every tissue", "0.007", mx, "le", note="medians " + ", ".join(f"{t} {v:+.5f}" for t, v in med.items()))
claim("B4b", g, "main", "The genome-wide median was within ±0.007 log₂ units of zero in every tissue", "0.007", mx, "le", note="Fig. 3B legend")
mean = {t: B[t].bias.mean() for t in TISSUES}
claim("B5a", g, "main", "values ranging from -0.028 to -0.019 log₂ units", "-0.028", min(mean.values()), note="FVB-arm shift = mean bias = mean(lfc_fvb) − mean(lfc_naive); " + ", ".join(f"{t} {v:+.4f}" for t, v in mean.items()))
claim("B5b", g, "main", "values ranging from -0.028 to -0.019 log₂ units", "-0.019", max(mean.values()))
fr = {t: 100 * (B[t].bias >= 0.5).mean() for t in TISSUES}
claim("B6a", g, "main", "from 0.56% in hypothalamus", "0.56", fr["hypothalamus"], note=", ".join(f"{t} {v:.3f}%" for t, v in fr.items()))
claim("B6b", g, "main", "to 0.92% in adipose tissue", "0.92", fr["adipose"])
claim("B6c", g, "main", "hypothalamus to 0.92% in adipose", "min", "min" if min(fr, key=fr.get) == "hypothalamus" and max(fr, key=fr.get) == "adipose" else "no", "eq", note="hypothalamus lowest, adipose highest")
u5 = set(ALL.loc[ALL.bias >= 0.5, "gene_id"])
claim("B7", g, "main", "326 unique genes showed bias ≥ 0.5", "326", len(u5), "int")
claim("B8", g, "main", "This represented 1.18% of the 27,732", "1.18", 100 * len(u5) / len(genes_any))
dw = {t: B[t].lfc_wasp.mean() - B[t].lfc_naive.mean() for t in TISSUES}
for t, txt in zip(TISSUES, ["+0.014", "+0.012", "+0.005", "+0.009"]):
    claim(f"B9_{t}", g, "main", f"{txt} in", txt, dw[t], note="mean(lfc_wasp) − mean(lfc_naive) over tested genes with a WASP estimate (WASP-arm table)")
claim("B9max", g, "main", "differed from naive alignment by less than 0.015 log₂ units", "0.015", max(abs(v) for v in dw.values()), "lt", note="Results and Fig. 3C legend")
claim("B9abs", g, "main", "WASP changed genome-wide mean log₂ fold changes by less than 0.015 units", "0.015", max(abs(v) for v in dw.values()), "lt", note="Abstract")
claim("B9d", g, "main", "changed genome-wide mean log₂ fold changes by less than 0.015 log₂ units", "0.015", max(abs(v) for v in dw.values()), "lt", note="Discussion")
# alternative: mean over genes with both arms only
dw2 = {t: (B[t].lfc_wasp - B[t].lfc_naive).mean() for t in TISSUES}
ROWS[-1]["note"] = "per-gene paired version: " + ", ".join(f"{t} {v:+.4f}" for t, v in dw2.items())
# SNP density
dens = G.density.reindex(genes_any)
claim("B10", g, "main", "The median genome-wide FVB/NJ SNP density was 0.26 SNPs/kb", "0.26", dens.median(), note="all FVB_NJ SNP records in the gene body / gene length (kb), tested genes")
claim("B10r", g, "resp", "median gene-body SNP density 0.26 SNPs/kb", "0.26", dens.median())

# ---------------------------------------------------------------- catalog -----------------------------------------------------
g = "C catalog and predictors"
q = ALL[ALL.bias >= 0.5].groupby("gene_id").size()
cat = sorted(gg for gg, k in q.items() if k >= 2 and gg in FVBGTF)
CAT = set(cat); mb_gene = per_gene.bias.mean()
claim("C1", g, "main", "This criterion identified 138 genes", "138", len(cat), "int", note=f"genes with bias ≥ 0.5 in ≥ 2 tissues: {int((q>=2).sum())}, of which in the converted GTF: {len(cat)}")
claim("C2", g, "main", "representing 0.50% of the 27,732", "0.50", 100 * len(cat) / len(genes_any))
nq = q.reindex(cat)
for k, txt, nd in ((4, "46", "Forty-six genes met the bias threshold in all four"), (3, "35", "Thirty-five genes met the threshold in three"), (2, "57", "Fifty-seven genes met the threshold in two")):
    claim(f"C3_{k}", g, "main", nd, txt, int((nq == k).sum()), "int")
claim("C4", g, "main", "16 catalog genes have a mean bias", "16", int((mb_gene.reindex(cat) < 0.5).sum()), "int")
for t, txt in zip(TISSUES, ["124", "92", "103", "84"]):
    claim(f"C5_{t}", g, "main", {"adipose": "Adipose tissue contained 124", "hypothalamus": "Hypothalamus contained 92", "liver": "Liver contained 103", "skin": "Skin contained 84"}[t], txt, len(set(B[t].index[B[t].bias >= 0.5]) & CAT), "int")
bt = G.biotype.reindex(cat)
grp = np.where(bt == "protein_coding", "pc", np.where(bt == "lncRNA", "lnc", np.where(bt.str.contains("pseudogene"), "ps", "other")))
for k, txt, nd in (("pc", "102", "accounted for 102 of the 138"), ("lnc", "15", "15 lncRNAs"), ("ps", "18", "18 pseudogenes"), ("other", "3", "3 genes of other biotypes")):
    claim(f"C6_{k}", g, "main", nd, txt, int((grp == k).sum()), "int")
# control-only catalog
W = []
for t in TISSUES:
    n, f = eg(t, "water", "naive", "unique"), eg(t, "water", "fvb", "unique"); d = (f.logFC - n.logFC).dropna(); W.append(d[d >= 0.5].index)
qw = pd.Series(np.concatenate([list(x) for x in W])).value_counts(); catw = {x for x, k in qw.items() if k >= 2 and x in FVBGTF}
claim("C7", g, "main", "alone (water-fed or control-cream animals; 132 genes", "132", len(catw), "int")
claim("C8", g, "main", "132 genes, 124 shared with the primary catalog", "124", len(catw & CAT), "int")
# predictors
pg = pd.DataFrame({"mean_bias": mb_gene, "density": dens, "length": G.length.reindex(genes_any), "biotype": G.biotype.reindex(genes_any),
                   "expr": per_gene.logCPM_naive.mean()})
r = stats.spearmanr(pg.density, pg.mean_bias)
claim("C9", g, "main", "(Spearman ρ = 0.337, p < 10⁻¹⁵)", "0.337", r.statistic); claim("C9p", g, "main", "(Spearman ρ = 0.337, p < 10⁻¹⁵)", "1e-15", r.pvalue, "plt")
claim("C9f", g, "Fig4_predictors", "0.337", "0.337", r.statistic)
a4 = pg[ntest.reindex(pg.index) == 4]; r4 = stats.spearmanr(a4.density, a4.mean_bias)
claim("C10", g, "main", "tissues, for which the mean bias is least noisy (ρ = 0.506)", "0.506", r4.statistic)
q90, q10, q75 = pg.density.quantile(0.9), pg.density.quantile(0.1), pg.density.quantile(0.75)
top, bot = pg[pg.density >= q90], pg[pg.density <= q10]
mw = stats.mannwhitneyu(top.mean_bias, bot.mean_bias)
claim("C11", g, "main", "Genes in the top 10% of SNP density had higher mean bias than genes in the bottom 10% (Mann–Whitney U test, p < 10⁻¹⁵)", "1e-15", mw.pvalue, "plt", note=f"top mean {top.mean_bias.mean():+.4f} vs bottom {bot.mean_bias.mean():+.4f}")
cm = pg.loc[cat, "density"].median(); ncm = pg.loc[[x for x in pg.index if x not in CAT], "density"].median()
claim("C12", g, "main", "Their median SNP density was 7.19 SNPs/kb", "7.19", cm); claim("C13", g, "main", "compared with 0.25 SNPs/kb for non-catalog", "0.25", ncm)
mw2 = stats.mannwhitneyu(pg.loc[cat, "density"], pg.loc[[x for x in pg.index if x not in CAT], "density"])
claim("C14", g, "main", "(Mann–Whitney U test, p < 10⁻¹⁵; Fig. 4A)", "1e-15", mw2.pvalue, "plt")
claim("C15", g, "Fig4_predictors", "2.56", "2.56", q75, note="75th percentile of SNP density over tested genes (inset threshold)")
claim("C16", g, "main", "The inset shows that 117 of 138 catalog genes fall in the top quartile", "117", int((pg.loc[cat, "density"] >= q75).sum()), "int")
tab = np.array([[int((pg.loc[cat, "density"] >= q90).sum()), int((pg.loc[cat, "density"] < q90).sum())],
                [int(((pg.density >= q90) & ~pg.index.isin(cat)).sum()), int(((pg.density < q90) & ~pg.index.isin(cat)).sum())]])
orr, pf = stats.fisher_exact(tab); se = math.sqrt((1 / tab).sum()); lo, hi = math.exp(math.log(orr) - 1.96 * se), math.exp(math.log(orr) + 1.96 * se)
claim("C17", g, "main", "(odds ratio = 10.0, 95% confidence interval = 7.2–14.1, p < 10⁻¹⁵)", "10.0", orr, note=f"2x2 {tab.tolist()}; 90th percentile {q90:.3f} SNPs/kb; catalog genes ≥ it: {tab[0,0]}")
claim("C17lo", g, "main", "95% confidence interval = 7.2–14.1", "7.2", lo, note="Woolf (log) interval")
claim("C17hi", g, "main", "95% confidence interval = 7.2–14.1", "14.1", hi)
claim("C17p", g, "main", "(odds ratio = 10.0, 95% confidence interval = 7.2–14.1, p < 10⁻¹⁵)", "1e-15", pf, "plt")
claim("C17q", g, "main", "", "6.79", q90, note="90th percentile threshold (Fisher test, Note/legend)")
claim("C17n", g, "main", "", "72", tab[0, 0], "int", note="catalog genes ≥ 90th percentile")
rl = stats.spearmanr(pg.length, pg.mean_bias); claim("C18", g, "main", "(Spearman ρ = 0.051, p < 10⁻¹⁵; Fig. 4B)", "0.051", rl.statistic); claim("C18p", g, "main", "(Spearman ρ = 0.051, p < 10⁻¹⁵; Fig. 4B)", "1e-15", rl.pvalue, "plt")
re_ = stats.spearmanr(pg.expr, pg.mean_bias); claim("C19", g, "main", "(Spearman ρ = -0.063, p < 10⁻¹⁵; Fig. 4C)", "-0.063", re_.statistic); claim("C19p", g, "main", "(Spearman ρ = -0.063, p < 10⁻¹⁵; Fig. 4C)", "1e-15", re_.pvalue, "plt")
bins = pd.qcut(pg.expr.rank(method="first"), 20, labels=False); bm = pg.groupby(bins).mean_bias.mean()
claim("C20", g, "main", "the binned mean bias stayed within ±0.07 log₂ units across the expression range", "0.07", bm.abs().max(), "le", note="20 equal-percentile bins of mean naive logCPM")
bg = np.where(pg.biotype == "protein_coding", "protein_coding", np.where(pg.biotype == "lncRNA", "lncRNA", np.where(pg.biotype.fillna("").str.contains("pseudogene"), "pseudogene", "other")))
pg["bgroup"] = bg; bmean = pg.groupby("bgroup").mean_bias.mean(); bmed = pg.groupby("bgroup").mean_bias.median()
claim("C21", g, "main", "Protein-coding genes and lncRNAs had mean bias values close to zero", "0.01", max(abs(bmean["protein_coding"]), abs(bmean["lncRNA"])), "le", note="means " + ", ".join(f"{k} {v:+.4f}" for k, v in bmean.items()))
claim("C22", g, "main", "Pseudogenes showed a more negative mean bias", "neg", "neg" if bmean["pseudogene"] < min(bmean["protein_coding"], bmean["lncRNA"]) else "no", "eq")
claim("C23", g, "main", "the median bias was within ±0.003 of zero across all biotype categories", "0.003", bmed.abs().max(), "le", note="medians " + ", ".join(f"{k} {v:+.5f}" for k, v in bmed.items()))
claim("C23b", g, "main", "Median bias was within ±0.003 of zero in all biotype categories", "0.003", bmed.abs().max(), "le", note="Fig. 4D legend")
kw = stats.kruskal(*[pg.mean_bias[pg.bgroup == k] for k in ("protein_coding", "lncRNA", "pseudogene", "other")])
claim("C24", g, "main", "(Kruskal–Wallis p < 10⁻¹⁵)", "1e-15", kw.pvalue, "plt")
nd = G.loc["ENSMUSG00000028648"] if "ENSMUSG00000028648" in G.index else G.loc[gid("Ndufs5")]
claim("C25", g, "main", "(chromosome 4; 7 FVB/NJ SNPs", "7", nd.snps_all, "int")   # Phase T: needle moved to the worked example; C25b (5.5-kb gene body) retired with the deleted Ndufs5 paragraph
claim("C25c", g, "main", "(chromosome 4; 7 FVB/NJ SNPs, 1.27 SNPs/kb)", "1.27", nd.density); claim("C25d", g, "main", "(green star; 1.3 SNPs/kb)", "1.3", nd.density)
claim("C25e", g, "main", "Ndufs5 (chromosome 4;", "4", nd.chrom, "eq")   # Phase T: needle moved to the worked example
def tb(name, t): return B[t].bias.get(gid(name), np.nan)
for t, txt in zip(TISSUES, ["+3.85", "+3.32", "+5.14", "+10.24"]): claim(f"C26_{t}", g, "main", f"{t} {txt}", txt, tb("Ndufs5", t))   # Phase T: needle moved to the worked example ("adipose +3.85, …")
# claim("C26m", g, "main", "with a four-tissue mean bias of +5.63", "+5.63", mb_gene[gid("Ndufs5")])   # retired in Phase T (Ndufs5 paragraph deleted)
for t, txt in zip(["adipose", "hypothalamus", "liver"], ["+5.03", "+5.60", "+4.28"]): claim(f"C27_{t}", g, "main", "The values were +5.03, +5.60, and +4.28", txt, tb("Gm25939", t))
claim("C27s", g, "main", "It was not detected in skin", "absent", "absent" if gid("Gm25939") not in B["skin"].index else "tested", "eq")

# ---------------------------------------------------------------- multi-mapping and mechanism classes ------------------------
g = "D multi-mapping"
mm = ALL.dropna(subset=["bias_multi"]).copy(); mm["delta"] = mm.bias_multi - mm.bias
claim("D1", g, "main", "Across 82,264 expressed gene-tissue pairs", "82264", len(mm), "int")
claim("D2", g, "main", "The median change in bias was +0.008 log₂ units", "+0.008", mm.delta.median()); claim("D3", g, "main", "The mean change was +0.021", "+0.021", mm.delta.mean())
claim("D4", g, "main", "Only 1.1% of gene-tissue pairs showed an absolute change greater than 0.5", "1.1", 100 * (mm.delta.abs() > 0.5).mean())
claim("D5", g, "af1", "and 0.7% |Δ| > 1", "0.7", 100 * (mm.delta.abs() > 1).mean())
claim("D6", g, "SFig2_multimapping", "474", "474", int(mm.gene_id.isin(CAT).sum()), "int", note="catalog gene × tissue pairs (S2 figure only)")
cg = pd.DataFrame({"mb": per_gene.bias.mean(), "mbm": mm.groupby("gene_id").bias_multi.mean()}).reindex(cat)
cg["delta"] = cg.mbm - cg.mb; cg["res"] = cg.mbm / cg.mb
mul, dirc = cg[cg.delta.abs() > 0.3], cg[cg.delta.abs() <= 0.3]
claim("D7", g, "main", "Eighty-two of the 138 catalog genes showed |Δ| > 0.3", "82", len(mul), "int"); claim("D8", g, "main", "The remaining 56 catalog genes", "56", len(dirc), "int")
claim("D9", g, "main", "The mean |Δ| in this class was 1.01", "1.01", mul.delta.abs().mean()); claim("D10", g, "main", "The mean |Δ| in this class was 0.17", "0.17", dirc.delta.abs().mean())
claim("D11", g, "main", "Thirty-seven of these 82 genes had residual fractions at or below 0.30", "37", int(((mul.mbm.abs() / mul.mb.abs()) <= 0.30).sum()), "int",
      note=f"Phase W-2: absolute ratio; the signed ratio gives {int((mul.res <= 0.30).sum())}")
claim("D12a", g, "af1", "residual fractions ranged from -9.51 to 2.26", "-9.51", mul.res.min()); claim("D12b", g, "af1", "residual fractions ranged from -9.51 to 2.26", "2.26", mul.res.max())
claim("D13a", g, "af1", "residual fractions 0.36–1.14", "0.36", dirc.res.min()); claim("D13b", g, "af1", "residual fractions 0.36–1.14", "1.14", dirc.res.max())
claim("D14", g, "af1", "Nineteen of these had residual fractions ≥ 0.80", "19", int((dirc.res >= 0.80).sum()), "int")
for nm, txt, nd_ in (("Rps2", "0.97", "Rps2 had a residual fraction of 0.97"), ("Rpl29", "0.89", "Rpl29 had a residual fraction of 0.89"), ("Klra1", "0.77", "Klra1 had a residual fraction of 0.77"), ("Ndufs5", "0.24", "(mean residual fraction 0.24)")):
    claim(f"D15_{nm}", g, "main", nd_, txt, cg.res[gid(nm)])
for nm in ("Rps2", "Rpl29", "Klra1"): claim(f"D16_{nm}", g, "main", "These included Rps2, Rpl29 and Klra1", "direct", "direct" if gid(nm) in dirc.index else "multi", "eq")
ds = B["skin"].loc[gid("Rps2")]
claim("D17", g, "main", "(residual fraction 0.98 in skin)", "0.98", ds.bias_multi / ds.bias)
rt = [B[t].loc[gid("Rps2")].bias_multi / B[t].loc[gid("Rps2")].bias for t in TISSUES]
# claim("D18a", …, "0.91", min(rt)); claim("D18b", …, "1.04", max(rt))   # retired in Phase T (per-tissue range removed from the Rps2 paragraph)
claim("D19", g, "af1", "(bias_multi = +5.23 log₂; residual fraction 0.98)", "+5.23", ds.bias_multi)
dclass = [(gg, t, B[t].loc[gg].bias) for gg in dirc.index for t in TISSUES if gg in B[t].index]
best = max(dclass, key=lambda x: x[2])
claim("D20", g, "main", "the largest single-tissue bias in this class (+5.36 log₂, skin)", "Rps2/skin", f"{G.gene_name[best[0]]}/{best[1]}", "eq")
claim("D21", g, "main", "the largest single-tissue bias in this class (+5.36 log₂, skin)", "+5.36", best[2])
s1 = mul.delta.abs().idxmax(); s2 = (dirc.res - 1).abs().idxmin()
claim("D22", g, "af1", "Ndufs5 and Mup17 represent the SNP-induced", "Ndufs5,Mup17", ",".join(G.gene_name[mul.delta.abs().sort_values(ascending=False).index[:2]]), "eq", note="rule: largest |Δ| in the multi-mapping class")
claim("D23", g, "af1", "3000002C10Rik and Rps2 represent the direct", "3000002C10Rik,Rps2", ",".join(G.gene_name[(dirc.res - 1).abs().sort_values().index[:2]]), "eq", note="rule: residual fraction closest to 1 in the direct class")

# ---------------------------------------------------------------- decomposition ------------------------------------------------
g = "E decomposition"
DMr = []
for t in TISSUES:   # bias_matched is defined wherever the naive and matched models both tested the gene; b6_loss also needs the FVB arm
    n_ = eg(t, "primary", "naive", "unique").logFC; m_ = pd.read_csv(f"{RES}/sensitivity_matched/edger/{t}.primary.matched.unique.tsv", sep="\t").set_index("gene_id").logFC
    f_ = eg(t, "primary", "fvb", "unique").logFC; bm_ = (m_ - n_).dropna()
    DMr.append(pd.DataFrame({"gene_id": bm_.index, "tissue": t, "bias_matched": bm_.values, "b6_loss": (f_.reindex(bm_.index) - m_.reindex(bm_.index)).values}))
dm = pd.concat(DMr, ignore_index=True)
claim("E1", g, "main", "median bias_matched −0.0001", "-0.0001", dm.bias_matched.median(), note=f"{len(dm):,} gene × tissue pairs with a matched estimate"); claim("E2", g, "main", "and b6_loss −0.0009", "-0.0009", ALL.b6_loss.median(), note="over tested pairs (naive, FVB and matched estimates)")
qm = dm[dm.bias_matched >= 0.5].groupby("gene_id").size(); catm = {x for x, k in qm.items() if k >= 2 and x in FVBGTF}
claim("E3", g, "main", "Applying the fixed catalog criterion to bias_matched alone yields 84 genes", "84", len(catm), "int")
claim("E4", g, "main", "84 genes, 57 of which are in the primary catalog", "57", len(catm & CAT), "int")
claim("E5", g, "main", "The other 27 of the 84 do not reach the primary criterion", "27", len(catm - CAT), "int")
both = sorted(catm & CAT); oth = sorted(CAT - catm)
claim("E6", g, "main", "In the remaining 81 catalog genes", "81", len(oth), "int")
pgd = pd.DataFrame({"bm": dm.groupby("gene_id").bias_matched.mean(), "bl": dm.dropna(subset=["b6_loss"]).groupby("gene_id").b6_loss.mean()})
mbias_all = per_gene.bias.mean()
pgd["mean_bias"] = mbias_all.reindex(pgd.index); pgd["share"] = pgd.bl / pgd.mean_bias
claim("E7", g, "main", "b6_loss exceeds bias_matched in 54 of the 81", "54", int((pgd.loc[oth, "bl"] > pgd.loc[oth, "bm"]).sum()), "int")
claim("E8", g, "main", "median share of the bias due to C57BL/6J-read loss 0.61", "0.61", pgd.loc[oth, "share"].median(), note="share = b6_loss_mean / mean_bias")
claim("E9", g, "main", "against 0.42 in the 57 genes", "0.42", pgd.loc[both, "share"].median()); claim("E10", g, "main", "and 0.50 over the whole catalog", "0.50", pgd.loc[cat, "share"].median())
nd5 = [B[t].loc[gid("Ndufs5")] for t in TISSUES]
claim("E11a", g, "main", "of this change, 1.47–4.61 log₂ units are FVB/NJ reads recovered", "1.47", min(x.bias_matched for x in nd5)); claim("E11b", g, "main", "of this change, 1.47–4.61 log₂ units are FVB/NJ reads recovered", "4.61", max(x.bias_matched for x in nd5))
claim("E11c", g, "main", "b6_loss +1.85 to +5.63", "+1.85", min(x.b6_loss for x in nd5)); claim("E11d", g, "main", "b6_loss +1.85 to +5.63", "+5.63", max(x.b6_loss for x in nd5))
claim("E11e", g, "main", "it meets the catalog criterion under both definitions of the bias", "4", sum(x.bias_matched >= 0.5 for x in nd5), "int")
rp = {t: B[t].loc[gid("Rps2")] for t in TISSUES}
for t, txt in zip(["adipose", "hypothalamus", "liver"], ["+0.73", "+1.06", "+0.86"]): claim(f"E12_{t}", g, "main", "bias_matched +0.73, +1.06 and +0.86 log₂", txt, rp[t].bias_matched)
claim("E12s", g, "main", "in skin bias_matched is −0.02 log₂", "-0.02", rp["skin"].bias_matched)
claim("E13", g, "main", "(b6_loss +5.38)", "+5.38", rp["skin"].b6_loss)
claim("E13f", g, "main", "(bias_matched +0.73 to +1.06 log₂)", "+0.73", min(rp[t].bias_matched for t in TISSUES[:3]))
# Gm26740: bias_matched_mean + b6_loss_mean vs mean_bias
gm = gid("Gm26740"); tp = [t for t in TISSUES if gm in B[t].index]; tm = sorted(dm[dm.gene_id == gm].tissue)
diffs = (pgd.bm + pgd.bl - pgd.mean_bias).reindex(cat).abs()
claim("E14", g, "af1", "for one gene, Gm26740, tested in three tissues under the primary model and four under the matched model", "Gm26740 only", ",".join(G.gene_name[diffs[diffs > 1e-6].index]) + " only" if (diffs > 1e-6).sum() == 1 else ",".join(G.gene_name[diffs[diffs > 1e-6].index]), "eq",
      note=f"bias_matched_mean + b6_loss_mean − mean_bias = {(pgd.bm + pgd.bl - pgd.mean_bias)[gm]:+.4f}; tissues tested (primary) {tp}; tissues with a matched estimate {tm} — the means are over different tissue sets, not a computation error")
claim("E15", g, "af1", "they differ from mean_bias by 0.02", "0.02", abs((pgd.bm + pgd.bl - pgd.mean_bias)[gm]), "num", note="S2 caption (M-3c)")
claim("E16", g, "af1", "tested in three tissues under the primary model and four under the matched model", "3/4", f"{len(tp)}/{len(tm)}", "eq")

# ---------------------------------------------------------------- worked examples ------------------------------------------------
g = "F worked examples"
nl = [B[t].loc[gid("Ndufs5")] for t in TISSUES]
claim("F1a", g, "main", "(log₂FC −0.92 to −5.06)", "-0.92", max(x.lfc_naive for x in nl)); claim("F1b", g, "main", "(log₂FC −0.92 to −5.06)", "-5.06", min(x.lfc_naive for x in nl))
claim("F2a", g, "main", "the contrast became +1.58 to +5.18", "+1.58", min(x.lfc_fvb for x in nl)); claim("F2b", g, "main", "the contrast became +1.58 to +5.18", "+5.18", max(x.lfc_fvb for x in nl))
claim("F3a", g, "main", "(bias_wasp +1.55 to +5.52)", "+1.55", min(x.bias_wasp for x in nl)); claim("F3b", g, "main", "(bias_wasp +1.55 to +5.52)", "+5.52", max(x.bias_wasp for x in nl))
claim("F3c", g, "main", "(bias +3.32 to +10.24 log₂; mean +5.63)", "+3.32", min(x.bias for x in nl))
claim("F4", g, "main", "(log₂FC −5.89) and the FVB arm reported log₂FC −0.53", "-5.89", rp["skin"].lfc_naive); claim("F5", g, "main", "the FVB arm reported log₂FC −0.53", "-0.53", rp["skin"].lfc_fvb)
claim("F6", g, "main", "at about 1/59 of the C57BL/6J level", "59", 2 ** -rp["skin"].lfc_naive, "num")
claim("F6b", g, "main", "(bias +5.36 log₂)", "+5.36", rp["skin"].bias)
for t, txt in zip(TISSUES, ["+1.12", "+1.36", "+1.14"]):
    claim(f"F7_{t}", g, "main", "(bias +1.12, +1.36 and +1.14 log₂", txt, rp[t].bias)
# per-animal CPM, Rps2 skin; normalisation: every arm divided by the naive-arm unique-only library size of the same animal
def cmat(t, arm, cnt): return pd.read_csv(f"{RES}/counts/{t}.{arm}.{cnt}.tsv", sep="\t").set_index("gene_id")
smp = pd.read_csv(f"{RES}/counts/skin.samples.tsv", sep="\t").set_index("sample_id")
fvbs, b6s = [s for s in smp.index if smp.strain[s] == "FVB/NJ"], [s for s in smp.index if smp.strain[s] == "C57BL/6J"]
lib = cmat("skin", "naive", "unique").sum()
r2 = gid("Rps2"); ids = {"Rps2": r2, "Rps2-ps10": gid("Rps2-ps10"), "Rps2-ps13": gid("Rps2-ps13"), "Rps2-ps6": gid("Rps2-ps6")}
lc = {}
for arm in ("naive", "wasp", "fvb"):
    m = cmat("skin", arm, "unique"); cpm = m.loc[r2] / lib * 1e6; lc[arm] = (np.log2(cpm[fvbs] + 1).mean(), np.log2(cpm[b6s] + 1).mean())
for arm, txt in zip(("naive", "wasp", "fvb"), ("2.2", "2.1", "2.2")): claim(f"F8_{arm}", g, "main", "(mean log₂(CPM + 1) 2.2 naive, 2.1 WASP, 2.2 FVB reference)", txt, lc[arm][0])
for arm, txt in zip(("naive", "wasp", "fvb"), ("7.7", "2.5", "2.6")): claim(f"F9_{arm}", g, "main", "C57BL/6J counts fell from 7.7 (naive) to 2.5 (WASP) and 2.6 (FVB reference)", txt, lc[arm][1])
# Rps2 / pseudogene CPM table (Additional file 1, Note 2), columns: naive uniq FVB, naive frac FVB, naive uniq B6, naive frac B6, fvb uniq FVB, fvb frac FVB, fvb uniq B6, fvb frac B6
tabtxt = {"Rps2": "3.71 13.00 209.76 232.82 3.67 41.84 5.06 19.90", "Rps2-ps10": "0.05 49.42 0.01 9.68 0.05 37.83 0.01 14.69",
          "Rps2-ps13": "0.13 44.41 0.08 7.24 0.13 32.58 0.11 9.48", "Rps2-ps6": "0.10 7.64 0.11 2.49 0.09 6.23 0.12 6.43"}
CPM = {}
for arm in ("naive", "fvb"):
    for cnt in ("unique", "multi"):
        m = cmat("skin", arm, cnt)
        for gn, gi in ids.items():
            v = m.loc[gi] / lib * 1e6; CPM[(gn, arm, cnt, "FVB")] = v[fvbs].mean(); CPM[(gn, arm, cnt, "B6")] = v[b6s].mean()
af1tab = re.sub(r"\s+", " ", TXT["af1"])
for gn, row in tabtxt.items():
    vals = row.split(); keys = [(gn, a, c, s) for a in ("naive", "fvb") for s in ("FVB", "B6") for c in ("unique", "multi")]
    for k, tv in zip(keys, vals):
        claim(f"F10_{gn}_{k[1]}_{k[2]}_{k[3]}", g, "af1", f"{gn} {row}", tv, CPM[k])
# claim("F11a", g, "main", "(fractional counts 49.4 and 44.4 counts per million in FVB/NJ skin against 13.0 at Rps2)", "49.4", CPM[("Rps2-ps10", "naive", "multi", "FVB")])   # retired in Phase T (sentence deleted)
# claim("F11b", g, "main", "(fractional counts 49.4 and 44.4", "44.4", CPM[("Rps2-ps13", "naive", "multi", "FVB")]); claim("F11c", g, "main", "against 13.0 at Rps2", "13.0", CPM[("Rps2", "naive", "multi", "FVB")])   # retired in Phase T (sentence deleted)
# for gn, txt in (("Rps2", "41.8"), ("Rps2-ps10", "37.8"), ("Rps2-ps13", "32.6")): claim(f"F12_{gn}", g, "main", "(41.8, 37.8 and 32.6)", txt, CPM[(gn, "fvb", "multi", "FVB")])   # retired in Phase T (sentence deleted)
def shift(arm): return math.log2(CPM[("Rps2", arm, "multi", "FVB")] / CPM[("Rps2", arm, "multi", "B6")]) - math.log2(CPM[("Rps2", arm, "unique", "FVB")] / CPM[("Rps2", arm, "unique", "B6")])
# Phase T: F13a/b, F22, F23, F24/F24n — the sentences were removed from the Results Rps2 paragraph and remain in Supplementary Note 2
claim("F13a", g, "af1", "(by 1.66 and 1.54 log₂ units)", "1.66", shift("naive")); claim("F13b", g, "af1", "(by 1.66 and 1.54 log₂ units)", "1.54", shift("fvb"), note="instruction listed 1.53; the text prints 1.54")
claim("F13c", g, "af1", "(by 1.66 and 1.54 log₂ units)", "1.54", shift("fvb"))
# Rps2 variants from the VCF / GTF
rg = G.loc[r2]
snp_bed = pd.read_csv(f"{C}/snps_all.bed", sep="\t", header=None, names=["c", "s0", "e"], dtype={"c": str}) if False else None
os.system(f"awk -v c={rg.chrom} -v s={rg.start0} -v e={rg.end} '$1==c && $3>s && $2<e' {C}/snps_all.bed > {C}/rps2_snps.bed")
os.system(f"awk -v c={rg.chrom} -v s={rg.start0} -v e={rg.end} '$1==c && $3>s && $2<e' {C}/indels_homalt.bed > {C}/rps2_indels.bed")
os.system(f"awk -v g={r2} '$4==g' {C}/exons_by_gene.bed | sort -k1,1 -k2,2n | bedtools merge -i - > {C}/rps2_exons.bed")
rs = pd.read_csv(f"{C}/rps2_snps.bed", sep="\t", header=None); ri = pd.read_csv(f"{C}/rps2_indels.bed", sep="\t", header=None, names=["c", "s0", "e", "ref", "alt"])
ex = pd.read_csv(f"{C}/rps2_exons.bed", sep="\t", header=None, names=["c", "s0", "e"])
def inex(p0, p1): return any((p1 > a) and (p0 < b) for a, b in zip(ex.s0, ex.e))
claim("F14", g, "main", "(ribosomal protein S2; chromosome 17; 53 FVB/NJ SNPs in a 3.8-kb gene body, 13.9 SNPs/kb)", "53", len(rs), "int")
claim("F15", g, "main", "53 FVB/NJ SNPs in a 3.8-kb gene body", "3.8", rg.length / 1000); claim("F16", g, "main", "13.9 SNPs/kb", "13.9", rg.density)
claim("F17", g, "af1", "(16 SNPs and 5 indels within exons; exonic: within the union of the exons of all Ensembl 110 Rps2 transcripts)", "16", sum(inex(a, b) for a, b in zip(rs[1], rs[2])), "int", note="exons = union of all Rps2 transcripts")
claim("F18", g, "af1", "53 FVB/NJ SNPs and 12 homozygous indels in 3.8 kb", "12", len(ri), "int")
exi = ri[[inex(a, b) for a, b in zip(ri.s0, ri.e)]]
claim("F19", g, "af1", "(16 SNPs and 5 indels within exons", "5", len(exi), "int", note="exonic hom-alt indels: " + "; ".join(f"{a+1} {r_}>{b_}" for a, r_, b_ in zip(exi.s0, exi.ref, exi.alt)))
# F20 (12-bp insertion): the sentence was removed from the deliverables in Phase L; the largest exonic insertion is checked in the report note of F19
al = pd.read_csv(f"{RES}/worked_examples/rps2_pseudogene_alleles.tsv", sep="\t"); sn = al[al.variant.str.startswith("SNP")]
ux = pd.read_csv(f"{RES}/worked_examples/rps2_exonic_variants_union.tsv", sep="\t")
my_snp = set(int(b) for a, b in zip(rs[1], rs[2]) if inex(a, b)); my_ind = set(int(a) + 1 for a, b in zip(exi.s0, exi.e))
claim("F17x", g, "af1", "", "same sets", "same sets" if (set(ux[ux.type == "SNP"].pos) == my_snp and set(ux[ux.type != "SNP"].pos) == my_ind) else "DIFFER", "eq", note="exonic SNP/indel positions recomputed here (VCF + GTF union) = positions in rps2_exonic_variants_union.tsv")
claim("F21", g, "main", "carry the FVB/NJ allele at all five exonic SNP positions", "5", int(((ux.type == "SNP") & (ux.in_copied_part_any == "Y")).sum()), "int", note="union SNPs inside the copied part; each of them is in rps2_pseudogene_alleles.tsv with the same call: " + str(all(len(sn[(sn.pseudogene == "Rps2-ps10") & (sn.variant == f"SNP 17:{p} {r_}>{a_}") & (sn.allele_carried != "not in pseudogene")]) == 1 for p, r_, a_ in zip(ux.pos[(ux.type=="SNP")&(ux.in_copied_part_any=="Y")], ux.ref[(ux.type=="SNP")&(ux.in_copied_part_any=="Y")], ux.alt[(ux.type=="SNP")&(ux.in_copied_part_any=="Y")]))))
for p_, txt in (("Rps2-ps10", "5"), ("Rps2-ps13", "5"), ("Rps2-ps6", "3")): claim(f"F22_{p_}", g, "af1", "(Rps2-ps6 at three;", txt, int((sn[sn.pseudogene == p_].allele_carried == "FVB/NJ").sum()), "int")
claim("F23", g, "af1", "the other eleven exonic SNPs", "11", len(my_snp) - int(((ux.type == "SNP") & (ux.in_copied_part_any == "Y")).sum()), "int")
claim("F24", g, "af1", "The five exonic indels of the FVB/NJ allele all lie outside that part", "0", int(((ux.type != "SNP") & (ux.in_copied_part_any == "Y")).sum()), "int")
claim("F24n", g, "af1", "The five exonic indels of the FVB/NJ allele", "5", len(my_ind), "int")

# ---------------------------------------------------------------- F1 ASE ------------------------------------------------------------
g = "G F1 allele-specific expression"
def f1af(path):
    d = pd.read_csv(path, sep="\t", dtype={"contig": str}); d = d[(d.totalCount >= 20) & (d.refCount >= 3) & (d.altCount >= 3)]
    d = d.assign(af=d.altCount / d.totalCount); return d
F1S = sorted({os.path.basename(f).split(".")[0] for f in glob.glob(f"{RES}/f1_ase/v2/*.naive.ase.tsv")})
ds_ = {s: ("GSE200632" if s.startswith("keratinocytes") else "GSE293440") for s in F1S}
nv = {s: f1af(f"{RES}/f1_ase/v2/{s}.naive.ase.tsv") for s in F1S}; wv = {s: f1af(f"{RES}/f1_ase/v2/{s}.wasp.ase.tsv") for s in F1S}
an = np.array([nv[s].af.mean() for s in F1S]); aw = np.array([wv[s].af.mean() for s in F1S])
def tt(x): t_, p_ = stats.ttest_1samp(x, 0.5); return x.mean(), x.std(ddof=1) / math.sqrt(len(x)), t_, p_
mn, se_, t_, p_ = tt(an)
claim("G1", g, "main", "across all 12 F1 samples was 0.480 ± 0.0003", "0.480", mn); claim("G2", g, "main", "0.480 ± 0.0003", "0.0003", se_, note="standard error over 12 animals")
claim("G3", g, "main", "(t = -61.7, p = 2.5 × 10⁻¹⁵; Fig. 2A)", "-61.7", t_); claim("G4", g, "main", "p = 2.5 × 10⁻¹⁵", "2.5e-15", p_, "sci")
claim("G4r", g, "resp", "(p = 2.5 × 10⁻¹⁵)", "2.5e-15", p_, "sci")
claim("G5", g, "main", "This corresponds to a 2.0-percentage-point deficit", "2.0", 100 * (0.5 - mn))
for gs, txt in (("GSE200632", "0.481"), ("GSE293440", "0.480")):
    claim(f"G6_{gs}", g, "main", "The mean FVB allele fraction was 0.481 in GSE200632 and 0.480 in GSE293440", txt, an[[ds_[s] == gs for s in F1S]].mean())
ns = np.array([len(nv[s]) for s in F1S])
claim("G7a", g, "main", "Each sample contributed 12,816–37,003 heterozygous sites", "12816", ns.min(), "int"); claim("G7b", g, "main", "12,816–37,003", "37003", ns.max(), "int")
claim("G8", g, "main", "The mean number of sites per sample was 26,966", "26966", ns.mean(), "num")
mw_, sew, tw, pw = tt(aw)
claim("G9", g, "main", "increased to 0.512 ± 0.0005", "0.512", mw_); claim("G10", g, "main", "0.512 ± 0.0005", "0.0005", sew)
claim("G11", g, "main", "(t = 25.0, p = 4.8 × 10⁻¹¹; Fig. 2A)", "25.0", tw); claim("G12", g, "main", "p = 4.8 × 10⁻¹¹", "4.8e-11", pw, "sci")
for gs, txt in (("GSE293440", "0.511"), ("GSE200632", "0.514")):
    claim(f"G12f_{gs}", g, "Fig2_f1_ase", txt, txt, aw[[ds_[s] == gs for s in F1S]].mean(), note=f"Fig. 2A per-dataset WASP mean {gs} (no text counterpart)")
# Fig. 2B: the GSE293440 animal with the median number of sites (lower median), sites retained in both arms; Fig. 2C: the eight animals pooled
cc = sorted([s for s in F1S if ds_[s] == "GSE293440"], key=lambda s: len(nv[s])); sB = cc[(len(cc) - 1) // 2]
bo = {s: nv[s].merge(wv[s], on=["contig", "position"], suffixes=("_n", "_w")) for s in cc}
claim("G29", g, "Fig2_f1_ase", "n = 32,067 sites (both arms)", "32067", len(bo[sB]), "int", note=f"animal {sB}")
claim("G30", g, "Fig2_f1_ase", "naive mean 0.482", "0.482", bo[sB].af_n.mean()); claim("G31", g, "Fig2_f1_ase", "WASP mean 0.511", "0.511", bo[sB].af_w.mean())
pool = pd.concat(bo.values())
claim("G32", g, "Fig2_f1_ase", "257,514", "257514", len(pool), "int", note="site × animal pairs retained in both arms, 8 GSE293440 animals")
claim("G33", g, "Fig2_f1_ase", "0.483", "0.483", pool.af_n.mean(), note="Fig. 2C pooled naive mean"); claim("G34", g, "Fig2_f1_ase", "0.511", "0.511", pool.af_w.mean(), note="Fig. 2C pooled WASP mean")
wr = pd.read_csv(f"{RES}/f1_ase/v2/f1_wasp_allele_removal.tsv", sep="\t")
claim("G13a", g, "main", "WASP removed 13–18% of the reads carrying the B6", "13", 100 * wr.frac_ref_reads_removed_by_wasp.min()); claim("G13b", g, "main", "WASP removed 13–18% of the reads", "18", 100 * wr.frac_ref_reads_removed_by_wasp.max())
claim("G14a", g, "main", "but only 1–2% of the reads carrying the FVB allele", "1", 100 * wr.frac_alt_reads_removed_by_wasp.min()); claim("G14b", g, "main", "but only 1–2% of the reads", "2", 100 * wr.frac_alt_reads_removed_by_wasp.max())
claim("G15a", g, "main", "so the FVB fraction rose at 69–78% of sites", "69", 100 * wr.frac_sites_af_increased.min()); claim("G15b", g, "main", "rose at 69–78% of sites", "78", 100 * wr.frac_sites_af_increased.max())
claim("G16", g, "main", "and fell at fewer than 1.1%", "1.1", 100 * wr.frac_sites_af_decreased.max(), "lt", note=f"maximum over animals {100*wr.frac_sites_af_decreased.max():.3f}%")
# catalog-gene sites: sites inside the gene body of a catalog gene
cb = G.loc[cat, ["chrom", "start0", "end"]]
iv = collections.defaultdict(list)
for c_, a_, b_ in cb.itertuples(index=False): iv[str(c_)].append((a_, b_))
def incat(d): return np.array([any(a_ < p <= b_ for a_, b_ in iv.get(c_, ())) for c_, p in zip(d.contig, d.position)])
cn, on, cw, ow, ncs = [], [], [], [], []
for s in F1S:
    m = incat(nv[s]); cn.append(nv[s].af[m].mean()); on.append(nv[s].af[~m].mean()); ncs.append(m.sum())
    m2 = incat(wv[s]); cw.append(wv[s].af[m2].mean()); ow.append(wv[s].af[~m2].mean())
cn, on, cw, ow = map(np.array, (cn, on, cw, ow))
claim("G17a", g, "main", "(230–450 per animal)", "230", min(ncs), "int"); claim("G17b", g, "main", "(230–450 per animal)", "450", max(ncs), "int")
claim("G18", g, "main", "the naive FVB allele fraction was 0.357, against 0.482", "0.357", cn.mean()); claim("G19", g, "main", "against 0.482 at the other sites", "0.482", on.mean())
claim("G20", g, "main", "p = 1.1 × 10⁻¹⁰; Fig. 2D", "1.1e-10", stats.ttest_rel(cn, on).pvalue, "sci")
claim("G21", g, "main", "a deficit of 14 percentage points", "14", 100 * (0.5 - cn.mean()))
claim("G22", g, "main", "Under WASP the two values were 0.497 and 0.512", "0.497", cw.mean()); claim("G23", g, "main", "Under WASP the two values were 0.497 and 0.512", "0.512", ow.mean())
claim("G24", g, "main", "(p = 0.017)", "0.017", stats.ttest_rel(cw, ow).pvalue)
claim("G25", g, "main", "is seven times the genome-wide deficit", "7", (0.5 - cn.mean()) / (0.5 - an.mean()), "num", note="ratio of the catalog-site deficit to the all-site deficit")
# submitted (old) site set, response letter
old_n = sum(1 for l in open(f"{P}/revision_v2/audit/f1_old/old_sites.bed") if l.strip() and not l.startswith(("#", "track")))
claim("G26", g, "resp", "at 94,351 positions", "94351", old_n, "int")
on_ = np.array([f1af(f"{RES}/f1_ase/v2/{s}.naive.oldsites.ase.tsv").af.mean() for s in F1S]); ow_ = np.array([f1af(f"{RES}/f1_ase/v2/{s}.wasp.oldsites.ase.tsv").af.mean() for s in F1S])
claim("G27", g, "resp", "reproduces the submitted values (0.329 and 0.382)", "0.329", on_.mean()); claim("G28", g, "resp", "(0.329 and 0.382)", "0.382", ow_.mean())

# ---------------------------------------------------------------- negative control (reference swap), from count matrices ----------
g = "H negative control"
def swap(t):
    n, f = cmat(t, "naive", "unique"), cmat(t, "fvb", "unique"); f = f.reindex(n.index).fillna(0)[n.columns]
    cpm = n / n.sum() * 1e6; keep = (cpm > 1).sum(axis=1) >= n.shape[1] / 2
    return np.log2((f[keep] + 0.5) / (n[keep] + 0.5))
SW = {t: swap(t) for t in TISSUES + ["mammary_tumour_PyVT"]}
STR = {}
for t in TISSUES: STR.update(pd.read_csv(f"{RES}/counts/{t}.samples.tsv", sep="\t").set_index("sample_id").strain.to_dict())
gmed = pd.concat([SW[t].median() for t in SW]); claim("H1", g, "main", "The genome-wide median of this ratio was 0.000 in every PyVT tumour", "0.000", gmed.abs().max(), "num")
fge = {}
for t in SW:
    for s in SW[t].columns: fge[s] = (SW[t][s] >= 0.5).mean()
def grp_of(s, t): return "PyVT" if t == "mammary_tumour_PyVT" else (STR[s] + (" skin" if t == "skin" else " PE"))
gg_ = collections.defaultdict(list)
for t in SW:
    for s in SW[t].columns: gg_[grp_of(s, t)].append(fge[s])
gm_ = {k: 100 * np.mean(v) for k, v in gg_.items()}
claim("H2a", g, "main", "with 0.3–0.5% of genes at or above +0.5 in every group", "0.3", min(gm_.values()), note="group means: " + ", ".join(f"{k} {v:.2f}%" for k, v in gm_.items()))
claim("H2b", g, "main", "with 0.3–0.5% of genes at or above +0.5 in every group", "0.5", max(gm_.values()))
cmed = {}; ncat_expr = {}
for t in SW:
    cg_ = [x for x in cat if x in SW[t].index]; ncat_expr[t] = len(cg_)
    for s in SW[t].columns: cmed[(t, s)] = SW[t].loc[cg_, s].median()
fv = [v for (t, s), v in cmed.items() if t in TISSUES and STR[s] == "FVB/NJ"]; b6 = [v for (t, s), v in cmed.items() if t in TISSUES and STR[s] == "C57BL/6J"]; py = [v for (t, s), v in cmed.items() if t == "mammary_tumour_PyVT"]
claim("H3a", g, "main", "(median reference-swap ratio per animal +0.24 to +0.49 over the catalog genes expressed", "+0.24", min(fv)); claim("H3b", g, "main", "+0.24 to +0.49", "+0.49", max(fv))
claim("H3n1", g, "main", "over the catalog genes expressed in the tissue (76–86 of 138)", "76", min(ncat_expr[t] for t in TISSUES), "int", note="catalog genes with naive CPM > 1 in ≥ half the samples: " + ", ".join(f"{t} {ncat_expr[t]}" for t in TISSUES))
claim("H3n2", g, "main", "over the catalog genes expressed in the tissue (76–86 of 138)", "86", max(ncat_expr[t] for t in TISSUES), "int")
claim("H4a", g, "main", "C57BL/6J animals lost reads (−0.39 to −0.69)", "-0.39", max(b6)); claim("H4b", g, "main", "C57BL/6J animals lost reads (−0.39 to −0.69)", "-0.69", min(b6))
claim("H4c", g, "main", "is smaller than in the C57BL/6J animals of the four tissues (−0.39 to −0.69)", "-0.69", min(b6))
claim("H5a", g, "main", "the PyVT tumours also lost reads (−0.08 to −0.16; 77 catalog genes expressed)", "-0.08", max(py)); claim("H5b", g, "main", "(−0.08 to −0.16; 77", "-0.16", min(py))
claim("H5n", g, "main", "77 catalog genes expressed", "77", ncat_expr["mammary_tumour_PyVT"], "int")
claim("H5d", g, "main", "(median reference-swap ratio −0.08 to −0.16 over the catalog genes)", "-0.16", min(py))
# held-out
HO = {}
for ho in TISSUES:
    oth3 = [t for t in TISSUES if t != ho]; qq = ALL[(ALL.tissue.isin(oth3)) & (ALL.bias >= 0.5)].groupby("gene_id").size()
    hc = [x for x, k in qq.items() if k >= 2 and x in FVBGTF]
    def med_(t, s): cg_ = [x for x in hc if x in SW[t].index]; return SW[t].loc[cg_, s].median()
    b6m = [med_(ho, s) for s in SW[ho].columns if STR[s] == "C57BL/6J"]; fvm = [med_(ho, s) for s in SW[ho].columns if STR[s] == "FVB/NJ"]
    pym = [med_("mammary_tumour_PyVT", s) for s in SW["mammary_tumour_PyVT"].columns]
    HO[ho] = dict(b6_mean=np.mean(b6m), b6_med=np.median(b6m), fv_mean=np.mean(fvm), fv_med=np.median(fvm), py_mean=np.mean(pym), py_med=np.median(pym), n=len(hc))
ho = pd.DataFrame(HO).T
note_ho = "group mean of per-animal medians: B6 " + ", ".join(f"{k} {v:+.3f}" for k, v in ho.b6_mean.items()) + "; median of per-animal medians: B6 " + ", ".join(f"{k} {v:+.3f}" for k, v in ho.b6_med.items())
claim("H6a", g, "main", "lost reads at those genes (group mean of the per-animal median −0.43 to −0.66 across the four held-out tissues)", "-0.43", ho.b6_mean.max(), note=note_ho)
claim("H6b", g, "main", "(group mean of the per-animal median −0.43 to −0.66 across the four held-out tissues)", "-0.66", ho.b6_mean.min(), note=f"as medians of per-animal medians: {ho.b6_med.max():+.3f} to {ho.b6_med.min():+.3f}")
claim("H7a", g, "main", "FVB/NJ animals gained (+0.29 to +0.44)", "+0.29", ho.fv_mean.min(), note=f"medians of per-animal medians: {ho.fv_med.min():+.3f} to {ho.fv_med.max():+.3f}")
claim("H7b", g, "main", "FVB/NJ animals gained (+0.29 to +0.44)", "+0.44", ho.fv_mean.max())
claim("H8a", g, "main", "and the PyVT tumours lost (−0.09 to −0.18)", "-0.09", ho.py_mean.max(), note=f"medians of per-animal medians: {ho.py_med.max():+.3f} to {ho.py_med.min():+.3f}")
claim("H8b", g, "main", "and the PyVT tumours lost (−0.09 to −0.18)", "-0.18", ho.py_mean.min())
# chr16 block
ps = pd.read_csv(f"{RES}/negative_control/pyvt_fullsite_chr3_10_16.per_site.tsv", sep="\t", dtype={"contig": str})
c16 = ps[(ps.contig == "16") & (ps.af >= 0.1) & ps.group.str.startswith("PyVT")]
per_t = c16.groupby("sample_id").position.agg(["min", "max", "size"])
claim("H9a", g, "main", "on chromosome 16 (37.7–49.7 Mb)", "37.7", c16.position.min() / 1e6, note=f"PyVT sites with FVB-allele fraction ≥ 0.1 on chr16; tumours with such sites: {len(per_t)}")
claim("H9b", g, "main", "(37.7–49.7 Mb)", "49.7", c16.position.max() / 1e6)
claim("H9c", g, "main", "in all 13 tumours", "13", len(per_t), "int")
b6c = ps[(ps.contig == "16") & (ps.af >= 0.1) & ~ps.group.str.startswith("PyVT") & (ps.position.between(c16.position.min(), c16.position.max()))]
claim("H9d", g, "main", "and in none of 16 C57BL/6J comparison animals", "0", b6c.sample_id.nunique(), "int")
inblk = [x for x in cat if str(G.chrom[x]) == "16" and G.end[x] > c16.position.min() and G.start0[x] < c16.position.max()]
claim("H10", g, "main", "two catalog genes lie in this block", "Cd200r4,Gm6030", ",".join(sorted(G.gene_name[inblk])), "eq")
cd = gid("Cd200r4"); claim("H11", g, "main", "(Cd200r4, whose PyVT reference-swap ratio was intermediate at −0.26", "-0.26", SW["mammary_tumour_PyVT"].loc[cd].median(), note="median over the 13 tumours")
claim("H12", g, "main", "Gm6030, not expressed in the tumours", "absent", "absent" if gid("Gm6030") not in SW["mammary_tumour_PyVT"].index else "expressed", "eq")
cg2 = [x for x in cat if x in SW["mammary_tumour_PyVT"].index]; cg3 = [x for x in cg2 if x not in inblk]
m_all = np.median([SW["mammary_tumour_PyVT"].loc[cg2, s].median() for s in SW["mammary_tumour_PyVT"].columns]); m_ex = np.median([SW["mammary_tumour_PyVT"].loc[cg3, s].median() for s in SW["mammary_tumour_PyVT"].columns])
claim("H13", g, "main", "left the PyVT catalog-gene median unchanged (−0.12)", "-0.12", m_ex, note=f"median over tumours of the per-tumour catalog median; with the two genes {m_all:+.4f}, without {m_ex:+.4f}")
# SNP-overlapping fraction from vW tags (naive arm)
vw = collections.defaultdict(list)
for f in glob.glob(f"{DATA}/*/*/qc/vW_tag_counts.tsv"):
    s = f.split("/")[-3]; t = f.split("/")[-4]
    if t not in TISSUES + ["mammary_tumour_PyVT"]: continue
    d = dict(l.split("\t") for l in open(f).read().split("\n") if l.strip()); tot = sum(int(v) for v in d.values()); tag = tot - int(d.get("none", 0))
    k = "PyVT" if t == "mammary_tumour_PyVT" else f"{STR[s]} {t}"; vw[k].append(100 * tag / tot)
vm = {k: np.mean(v) for k, v in vw.items()}; oth_ = [v for k, v in vm.items() if k != "PyVT"]
claim("H14", g, "main", "overlapping an FVB/NJ SNP (7.7%)", "7.7", vm["PyVT"]); claim("H15a", g, "main", "the other libraries (4.8–10.8%)", "4.8", min(oth_)); claim("H15b", g, "main", "(4.8–10.8%)", "10.8", max(oth_), note="group means by tissue × strain")

# ---------------------------------------------------------------- WASP attrition from featureCounts summaries ----------------------
g = "I WASP attrition"
att = collections.defaultdict(list)
for t in TISSUES:
    for s in pd.read_csv(f"{RES}/counts/{t}.samples.tsv", sep="\t").sample_id:
        def asg(arm): return int([l.split("\t")[1] for l in open(f"{DATA}/{t}/{s}/counts/{s}.{arm}.unique.txt.summary") if l.startswith("Assigned\t")][0])
        att[(t, STR[s])].append(100 * (1 - asg("wasp") / asg("naive")))
am = {k: np.mean(v) for k, v in att.items()}
claim("I1a", g, "main", "the removed fraction ranged from 0.23% to 2.55%", "0.23", min(am.values()), note="group means (tissue × strain) of per-animal 100 × (1 − WASP assigned / naive assigned): " + ", ".join(f"{k[0]}/{k[1]} {v:.2f}" for k, v in am.items()))
claim("I1b", g, "main", "the removed fraction ranged from 0.23% to 2.55%", "2.55", max(am.values()))
claim("I2a", g, "main", "(C57BL/6J 1.19–2.55% versus FVB/NJ 0.23–0.37%)", "1.19", min(v for k, v in am.items() if k[1] == "C57BL/6J")); claim("I2b", g, "main", "(C57BL/6J 1.19–2.55%", "2.55", max(v for k, v in am.items() if k[1] == "C57BL/6J"))
claim("I2c", g, "main", "versus FVB/NJ 0.23–0.37%)", "0.23", min(v for k, v in am.items() if k[1] == "FVB/NJ")); claim("I2d", g, "main", "versus FVB/NJ 0.23–0.37%)", "0.37", max(v for k, v in am.items() if k[1] == "FVB/NJ"))
dd_ = [am[(t, "FVB/NJ")] - am[(t, "C57BL/6J")] for t in TISSUES]
claim("I3a", g, "main", "ranged from -2.23 to -0.89 percentage points", "-2.23", min(dd_)); claim("I3b", g, "main", "ranged from -2.23 to -0.89 percentage points", "-0.89", max(dd_))
claim("I4", g, "main", "WASP removed only 0.23–2.55% of assigned reads", "2.55", max(am.values()))

# ---------------------------------------------------------------- within-strain ---------------------------------------------------
g = "J within-strain"
WS = {(t, a): pd.read_csv(f"{RES}/within_fvb/{t}.diet.{a}.tsv", sep="\t").set_index("gene_id") for t in TISSUES for a in ("naive", "wasp", "fvb")}
de = {k: int((v.FDR < 0.05).sum()) for k, v in WS.items()}
claim("J1", g, "main", "no gene was differentially expressed between fructose- and water-fed FVB/NJ mice at FDR < 0.05 under the naive and WASP arms", "0", sum(de[(t, a)] for t in TISSUES[:3] for a in ("naive", "wasp")), "int")
claim("J2", g, "main", "three genes in liver under the FVB arm", "3", de[("liver", "fvb")], "int", note=f"adipose fvb {de[('adipose','fvb')]}, hypothalamus fvb {de[('hypothalamus','fvb')]}")
for a, txt in (("naive", "598"), ("wasp", "593"), ("fvb", "539")): claim(f"J3_{a}", g, "main", "(598 genes under naive, 593 under WASP, 539 under the FVB reference", txt, de[("skin", a)], "int")
sn_, sf_ = WS[("skin", "naive")], WS[("skin", "fvb")]; sh = set(sn_.index[sn_.FDR < 0.05]) & set(sf_.index[sf_.FDR < 0.05])
claim("J4", g, "main", "519 shared and concordant in direction", "519", sum(np.sign(sn_.logFC[x]) == np.sign(sf_.logFC[x]) for x in sh), "int", note=f"shared {len(sh)}")
rr, f05, cmd = [], [], []
for t in TISSUES:
    a, b = WS[(t, "naive")], WS[(t, "fvb")]; com = a.index.intersection(b.index); d_ = b.logFC[com] - a.logFC[com]
    rr.append(stats.pearsonr(a.logFC[com], b.logFC[com]).statistic); f05.append(100 * (d_.abs() >= 0.5).mean()); cmd.append(d_[[x for x in com if x in CAT]].abs().median())
claim("J5a", g, "main", "(Pearson r = 0.991–0.997;", "0.991", min(rr)); claim("J5b", g, "main", "(Pearson r = 0.991–0.997;", "0.997", max(rr))
# claim("J5c", g, "main", "nearly identical gene by gene (Pearson r ≥ 0.99)", "0.99", min(rr), "ge")   # retired in Phase T (sentence deleted)
claim("J6a", g, "main", "0.09–0.38% of genes differed by ≥ 0.5 log₂", "0.09", min(f05)); claim("J6b", g, "main", "0.09–0.38% of genes differed", "0.38", max(f05))
claim("J7a", g, "main", "(median |Δ log₂FC| 0.02–0.08)", "0.02", min(cmd)); claim("J7b", g, "main", "(median |Δ log₂FC| 0.02–0.08)", "0.08", max(cmd))

# ---------------------------------------------------------------- letters / survey ---------------------------------------------------
g = "K letters and survey"
old = pd.read_csv(f"{P}/legacy/results/stage6/catalog_final.tsv", sep="\t")
claim("K1", g, "resp", "Of the 32 genes of the submitted catalog", "32", len(old), "int")
claim("K2", g, "resp", "24 are in the new catalog", "24", int(old.gene_id.isin(CAT).sum()), "int")
claim("K3", g, "resp", "26 reach the threshold in at least one tissue", "26", int(old.gene_id.isin(u5).sum()), "int")
negall = [nm for nm in ("Fau", "Rsph3a") if all(B[t].bias.get(gid(nm), np.nan) < 0 for t in TISSUES)]
claim("K4", g, "resp", "(one of them now shows negative bias in all four tissues)", "1", len(negall), "int", note=f"submitted worked examples Fau, Rsph3a; negative in all four: {negall}; in new catalog: {[nm for nm in ('Fau','Rsph3a') if gid(nm) in CAT]}")
sv = pd.read_csv(f"{RES}/survey/survey_of_practice_rescreen.tsv", sep="\t"); vc = sv.category.value_counts()   # Phase W-1: the re-screened survey (was survey_of_practice_v2.tsv, 154 series)
claim("K5", g, "main", "We retrieved data-processing metadata from 154 GEO series", "154", len(sv), "int")
claim("K6", g, "main", "Among these 154 series, 140 reported", "140", vc.get("stated_B6", 0), "int"); claim("K7", g, "main", "The remaining 14 studies did not report", "14", vc.get("not_reported", 0), "int")
claim("K8", g, "main", "and none reported a strain-aware alignment strategy", "0", vc.get("strain_aware", 0), "int")
claim("K8r", g, "resp", "all 140 FVB RNA-seq studies that reported a reference used a B6 assembly", "140", vc.get("stated_B6", 0), "int")
for k in ("stated_B6", "not_reported", "strain_aware"):   # Phase X: the printed percentage is read from the regenerated Fig. 1A at the position of the count
    n_ = int(vc.get(k, 0)); m_ = re.search(rf"\b{n_} \((\d+\.\d)%\)", FIG["Fig1_survey"]); txt = m_.group(1) if m_ else "absent"
    claim(f"K9_{k}", g, "Fig1_survey", f"{n_} ({txt}%)", txt, 100 * n_ / len(sv))
claim("K9_total", g, "Fig1_survey", f"{len(sv)} FVB RNA-seq series", str(len(sv)), len(sv), "int", note="Fig. 1A title")
# Phase W-3: animals of the strain comparison, from the sample identifiers of the manifest (fructose study: strain letter, diet
# letter, tissue code, animal number numbered per strain; skin: one sample per animal); not from the manifest animal_id column
man = pd.read_csv(f"{P}/revision_v2/manifest/manifest.tsv", sep="\t", dtype=str, keep_default_na=False); man = man[(man.include == "Y") & (man.role == "strain_comparison")]
def _animal(sid):
    t = sid.split("_")[-1]; m = re.match(r"^([CF])[FW](?:M|Hyp|L)(\d+)$", t)
    return f"{m.group(1)}{m.group(2)}" if m else t
an_ = {st: man[man.strain == st].sample_id.map(_animal).nunique() for st in ("FVB/NJ", "C57BL/6J")}
claim("W3a", g, "main", "(36 FVB/NJ and 34 C57BL/6J animals)", "36", an_["FVB/NJ"], "int", note=f"FVB/NJ samples {int((man.strain=='FVB/NJ').sum())}")
claim("W3b", g, "main", "(36 FVB/NJ and 34 C57BL/6J animals)", "34", an_["C57BL/6J"], "int", note=f"C57BL/6J samples {int((man.strain=='C57BL/6J').sum())}")
claim("W3c", g, "main", "All 83 animals of the three-arm analysis", "83", an_["FVB/NJ"] + an_["C57BL/6J"] + 13, "int", note="13 PyVT mice, one tumour each")
claim("K10", g, "main", "This search returned 1,047 GEO series", "1047", "", "nr", note="search-hit count: no search log is recorded in the project (the value is carried verbatim from the submitted text); cannot be recomputed")
claim("K11", g, "main", "taken from a 290-sample series", "290", "", "nr", note="GEO series size of GSE293440: 290 !Series_sample_id entries in the GEO SOFT record (checked online 2026-10-23); not a result file")
claim("K12", g, "Fig1_survey", "naive B6 mapping 0.480; WASP 0.512 (Fig. 2)", "0.480", an.mean(), "num", note="Fig. 1B regenerated in Phase M-3 (D1)")
claim("K13", g, "Fig1_survey", "naive B6 mapping 0.480; WASP 0.512 (Fig. 2)", "0.512", aw.mean(), "num")
claim("K13b", g, "Fig1_survey", "", "absent", "absent" if ("0.333" not in FIG["Fig1_survey"] and "0.386" not in FIG["Fig1_survey"]) else "present", "eq", note="no submitted-version F1 value left in Fig. 1")
v1 = pd.read_csv(f"{RES}/f1_ase/f1_ase_per_sample.tsv", sep="\t"); v1n, v1w = v1[v1.arm == "naive"].mean_fvb_fraction.mean(), v1[v1.arm == "wasp"].mean_fvb_fraction.mean()
claim("K16", g, "resp", "0.480 (not 0.333)", "0.333", v1n, note="submitted-version value quoted; submitted per-animal table f1_ase/f1_ase_per_sample.tsv (old site set)")
claim("K17", g, "resp", "0.512 (not 0.386)", "0.386", v1w); claim("K18", g, "cover", "0.333", "0.333", v1n, note="submitted-version value quoted in the cover letter")
claim("K14", g, "resp", "the FVB allele fraction is 0.480 (not 0.333)", "0.480", an.mean()); claim("K15", g, "resp", "it is 0.512 (not 0.386)", "0.512", aw.mean())

g = "L Supplementary Fig. S1"
T1 = f"{RES}/identity/full_site_tables"
fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(f"{P}/ref/b6_fasta/GRCm39.primary.fa.fai")}
def binfrac(path):
    a = pd.read_csv(path, sep="\t", dtype={"contig": str}); a = a[((a.refCount + a.altCount) >= 20) & a.contig.isin([str(i) for i in range(1, 20)])]
    a = a.assign(af=a.altCount / (a.refCount + a.altCount), b=a.position // 20_000_000)
    q = a.groupby(["contig", "b"]).af.agg(n="size", hi=lambda x: (x > 0.9).mean()); return q[q.n >= 10].hi
hs = [binfrac(f"{T1}/heart_B6_C57_M_12wk_mRNA_{i}.ase.tsv") for i in range(1, 5)]
hm = pd.concat(hs, axis=1).dropna(); shared = sorted({c for (c, b), r in hm.iterrows() if (r >= 0.5).all()}, key=int)
claim("L1", g, "af1", "on chromosomes 2, 3, 8, 10 and 11 at least half of the sites in a bin", "2, 3, 8, 10 and 11", ", ".join(shared[:-1]) + " and " + shared[-1], "eq")
claim("L1b", g, "af1", "(for example chromosome 12, 0–20 Mb, in one animal)", "1", int(sum(h_.get(("12", 0), 0) >= 0.5 for h_ in hs)), "int")
cnd = []
for f_ in glob.glob(f"{DATA}/*/*/qc/*.ase.tsv"):
    s_ = os.path.basename(f_)[:-8]
    if s_.startswith("heart") or s_ not in STR: continue
    d_ = pd.read_csv(f_, sep="\t", dtype={"contig": str})
    if all((c, p) in sub for c, p in zip(d_.contig, d_.position)): continue
    cnd.append((STR[s_], -int((d_.totalCount >= 20).sum()), s_))
pick = {st: sorted(x for x in cnd if x[0] == st)[0][2] for st in ("C57BL/6J", "FVB/NJ")}
claim("L2", g, "af1", "(both mesenteric adipose tissue, GSE123875; chosen as the animal of each strain with the most sites of depth ≥ 20)", "adipose_B6_CWM8,adipose_FVB_FWM9", f"{pick['C57BL/6J']},{pick['FVB/NJ']}", "eq", note=f"{len(cnd)} candidates")
b6r, fvr = binfrac(f"{T1}/{pick['C57BL/6J']}.ase.tsv"), binfrac(f"{T1}/{pick['FVB/NJ']}.ase.tsv")
claim("L3", g, "af1", "The C57BL/6J reference animal has at most 0.033 of sites", "0.033", b6r.max()); claim("L3b", g, "af1", "has at most 0.033", "0.033", b6r.max(), "le")
claim("L4", g, "af1", "the FVB/NJ reference animal at least 0.69", "0.69", fvr.min()); claim("L4b", g, "af1", "at least 0.69", "0.69", fvr.min(), "ge")
g = "N figure text (Phase N)"
claim("N1", g, "Fig3_strain_comparison", "catalog genes tested in liver: 117 of 138", "117", len(set(B["liver"].index) & CAT), "int", note="Fig. 3A key (N-2b)")
for k_, n_ in ((2, "40"), (3, "31"), (4, "46")):
    claim(f"N2_{k_}", g, "Fig3_strain_comparison", f"(n = {n_} tested in liver)", n_, int(sum(1 for x in set(B["liver"].index) & CAT if int(nq[x]) == k_)), "int")
claim("N3", g, "Fig4_predictors", "Spearman ρ = 0.337, p < 1e-15", "1e-15", r.pvalue, "plt", note="Fig. 4A (N-1b): no underflowed p = 0")
claim("N3b", g, "Fig4_predictors", "", "absent", "absent" if not re.search(r"p = 0\b|0\.0e\+00|0e\+00", FIG["Fig4_predictors"]) and all(not re.search(r"p = 0\b|0\.0e\+00|0e\+00", t_) for t_ in FIG.values()) else "present", "eq", note="no printed p of 0 in any figure")
claim("N4", g, "Fig6_worked_examples", "(skin +5.36 log₂)", "+5.36", rp["skin"].bias, note="Fig. 6D title (N-1a); title lines: 'largest single-tissue' / 'bias of the direct-misalignment class' / '(skin +5.36 log₂)'")
claim("N4b", g, "Fig6_worked_examples", "bias of the direct-misalignment class", "qualified", "qualified" if "largest single-tissue bias\n" not in FIG["Fig6_worked_examples"] else "unqualified", "eq")
_bt = G.biotype.reindex(genes_any)
for k_, (lab_, sel_) in enumerate((("17,486", _bt == "protein_coding"), ("5,403", _bt == "lncRNA"), ("2,169", _bt.fillna("").str.contains("pseudogene")), ("141", _bt == "miRNA"))):
    claim(f"N6_{k_}", g, "Fig4_predictors", "", lab_, int(sel_.sum()), "int", note="Fig. 4D biotype n (rotated tick label; pdftotext may split the digits)")
claim("N6_4", g, "Fig4_predictors", "", "2,533", int(len(_bt) - (_bt == "protein_coding").sum() - (_bt == "lncRNA").sum() - _bt.fillna("").str.contains("pseudogene").sum() - (_bt == "miRNA").sum()), "int")
claim("N5", g, "af1", "lie outside the part that the pseudogenes copy", "11", len(my_snp) - int(((ux.type == "SNP") & (ux.in_copied_part_any == "Y")).sum()), "int", note="N-2a: 10 in exon segments absent from the copied transcript ENSMUST00000146867, 1 (17:24,940,701) in that transcript just past the copied span")
g = "P Phase P (genome-level differences, clusters)"
# exonic SNP count per tested gene: union of all exons of the gene (GTF, cache/exons_by_gene.bed) x every FVB_NJ SNP record (cache/snps_all.bed)
if not os.path.exists(f"{C}/gene_exon_union_snps.tsv"):
    os.system(f"sort -k4,4 -k1,1 -k2,2n -S 2G {C}/exons_by_gene.bed | awk -F'\\t' 'BEGIN{{OFS=\"\\t\"}} {{if($4==g && $1==c && $2<=e){{if($3>e)e=$3}} else {{if(g!=\"\")print c,s,e,g; c=$1;s=$2;e=$3;g=$4}}}} END{{print c,s,e,g}}' | sort -k1,1 -k2,2n -S 2G | bedtools intersect -c -sorted -a - -b {C}/snps_all.bed | awk -F'\\t' '{{n[$4]+=$5}} END{{for(k in n) print k\"\\t\"n[k]}}' > {C}/gene_exon_union_snps.tsv")
ge = pd.read_csv(f"{C}/gene_exon_union_snps.tsv", sep="\t", header=None, names=["gene_id", "n"]).set_index("gene_id").n.reindex(genes_any)
claim("P1", g, "main", "48.1% of the genes tested here carry no exonic FVB/NJ SNP", "48.1", 100 * (ge == 0).mean(), note=f"{int((ge == 0).sum()):,} of {len(ge):,} tested genes; exons = union of all Ensembl 110 exons of the gene; SNPs = every FVB_NJ SNP record")
wx = open(f"{P}/revision_v2/audit/wong2012_excerpts.txt").read()
claim("P2", g, "main", "(30,048 of at least 50 bp in the first FVB/NJ genome study [4])", "30048", "30048" if "catalogue of 30,048 structural variants (≥50 bp)" in wx else "not in source", "eq", note="literature value; verified against the verbatim source sentence (revision_v2/audit/wong2012_excerpts.txt, Wong et al. 2012)")
claim("P2b", g, "main", "the Pde6b nonsense mutation that causes early-onset retinal degeneration and the 2-bp Hc deletion that causes C5 deficiency [4]", "2", "2" if "2-bp deletion in the Hc gene" in wx and "nonsense mutation in the Pde6b gene" in wx else "not in source", "eq", note="literature statement checked against the source sentences")
claim("P2c", g, "resp", "about 50-fold Illumina coverage aligned to the C57BL/6J reference", "50", "50" if "approximately 50× coverage" in wx and "NCBIM37/mm9" in wx else "not in source", "eq", note="response letter; Wong et al. 2012")
# clusters: catalog genes (independent CAT), Ensembl gene start = start0 + 1, single linkage with gaps < 1 Mb, clusters of >= 3 genes
cg_ = G.loc[cat, ["chrom", "start0", "gene_name"]].assign(start=lambda d: d.start0 + 1).sort_values(["chrom", "start"])
CL = []
for ch, sdf in cg_.groupby("chrom"):
    cur = []
    for gi_, st in zip(sdf.index, sdf.start):
        if cur and st - cur[-1][1] >= 1_000_000:
            if len(cur) >= 3: CL.append((ch, [x for x, _ in cur]))
            cur = []
        cur.append((gi_, st))
    if len(cur) >= 3: CL.append((ch, [x for x, _ in cur]))
ing = [x for _, l in CL for x in l]
def fam(pred): sel = [l for ch, l in CL if any(pred(G.gene_name[x], str(ch)) for x in l)]; return sum(len(l) for l in sel), len(sel)
claim("P3", g, "main", "39 of the 138 lie in 10 clusters of at least three catalog genes", "10", len(CL), "int", note="; ".join(f"chr{ch}: {', '.join(G.gene_name[x] for x in l)}" for ch, l in CL))
claim("P3b", g, "main", "39 of the 138 lie in 10 clusters", "39", len(ing), "int")
claim("P3c", g, "main", "26 of these 39 genes are in the SNP-induced multi-mapping class", "26", sum(1 for x in ing if x in mul.index), "int")
for k_, txt_, pr_ in (("mhc", "9", lambda n, ch: ch == "17" and n.startswith("H2-")), ("hist", "11", lambda n, ch: ch == "13" and bool(re.match(r"H(2a|2b|3|4)c\d", n))), ("mup", "3", lambda n, ch: n.startswith("Mup")), ("ifi", "3", lambda n, ch: n.startswith("Ifi"))):
    n_, c_ = fam(pr_); claim(f"P4_{k_}", g, "main", {"mhc": "the MHC on chromosome 17 (9 genes in two clusters)", "hist": "(11 genes in two clusters)", "mup": "the major urinary protein cluster on chromosome 4 (3 genes)", "ifi": "the Ifi200 cluster on chromosome 1 (3 genes)"}[k_], txt_, n_, "int", note=f"{c_} cluster(s)")
s2 = pd.read_excel(f"{SUB}/supplementary_tables_v2.xlsx", sheet_name="S2", dtype=str).fillna("")
claim("P5", g, "af1", "", "39", int((s2.cluster != "").sum()), "int", note="Supplementary Table S2, column cluster: genes with a cluster label")
claim("P5b", g, "af1", "", "same genes", "same genes" if set(s2.gene_id[s2.cluster != ""]) == set(ing) else "DIFFER", "eq", note="S2 cluster genes = clusters recomputed here")
fs_ = pd.read_csv(f"{P}/revision_v2/freeze/figure_sizes.tsv", sep="\t")
claim("N7", "N figure text (Phase N)", "resp", "(170 mm, full page width; Supplementary Fig. S3, 85 mm)", "170", sorted(set(fs_[fs_.figure != "SFig3_bias_decomposition"].width_mm.round(0).astype(int)))[0] if fs_[fs_.figure != "SFig3_bias_decomposition"].width_mm.round(0).nunique() == 1 else -1, "int", note="width of every full-width figure as saved (freeze/figure_sizes.tsv)")
claim("N7b", "N figure text (Phase N)", "resp", "(170 mm, full page width; Supplementary Fig. S3, 85 mm)", "85", int(round(float(fs_[fs_.figure == "SFig3_bias_decomposition"].width_mm.iloc[0]))), "int")
claim("N7c", "N figure text (Phase N)", "resp", "with all lettering at least 7 pt at that size", "7", float(fs_.min_text_pt.min()), "ge", note="smallest text over all figures (freeze/figure_sizes.tsv)")
g = "Q Phase Q (decomposition wording, within-strain wording)"
claim("Q1a", g, "main", "of this change, 1.47–4.61 log₂ units are FVB/NJ reads recovered on their own reference (bias_matched)", "1.47", min(x.bias_matched for x in nd5))
claim("Q1b", g, "main", "of this change, 1.47–4.61 log₂ units", "4.61", max(x.bias_matched for x in nd5))
claim("Q1c", g, "main", "the rest is loss of C57BL/6J reads (b6_loss +1.85 to +5.63)", "+1.85", min(x.b6_loss for x in nd5)); claim("Q1d", g, "main", "the rest is loss of C57BL/6J reads (b6_loss +1.85 to +5.63)", "+5.63", max(x.b6_loss for x in nd5))
claim("Q2a", g, "main", "at 57 of the 138 catalog genes FVB/NJ expression is demonstrably underestimated", "57", len(catm & CAT), "int")
claim("Q2b", g, "main", "at 57 of these genes FVB/NJ reads were demonstrably recovered when only the FVB/NJ animals were moved to their own reference", "57", len(catm & CAT), "int", note="Abstract")
claim("Q2c", g, "main", "it demonstrably recovers FVB/NJ reads at 57 of the 138 catalog genes", "57", len(catm & CAT), "int")
claim("Q2d", g, "main", "at 57 of them a demonstrable underestimation of FVB/NJ expression", "57", len(catm & CAT), "int")
claim("Q2e", g, "cover", "demonstrably recovers FVB/NJ reads at 57 of the 138 catalog genes", "57", len(catm & CAT), "int")
claim("Q2f", g, "cover", "mean FVB allele fraction 0.480 rather than 0.333", "0.480", an.mean()); claim("Q2g", g, "cover", "(0.512 rather than 0.386)", "0.512", aw.mean())
claim("Q2h", g, "cover", "(0.512 rather than 0.386)", "0.386", v1w, note="submitted-version value quoted")
# claim("Q3a", g, "main", "(598 under the naive arm and 539 under the FVB/NJ reference in skin)", "598", de[("skin", "naive")], "int"); claim("Q3b", g, "main", "(598 under the naive arm and 539 under the FVB/NJ reference in skin)", "539", de[("skin", "fvb")], "int")   # retired in Phase T (sentence deleted)
# claim("Q4", g, "main", "all 140 FVB RNA-seq studies that reported a reference used a B6 assembly", "140", vc.get("stated_B6", 0), "int")   # retired in Phase T (sentence deleted)
out = pd.DataFrame(ROWS); out.to_csv(f"{P}/revision_v2/freeze/recompute_report.tsv", sep="\t", index=False)
print(out.status.value_counts().to_string()); print(out.groupby(["group", "status"]).size().unstack(fill_value=0).to_string())
print(out[out.status != "PASS"][["claim_id", "location", "text_value", "exact", "at_text_precision", "status", "note"]].to_string(max_colwidth=140))
