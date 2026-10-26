#!/usr/bin/env python3
"""Phase W-1 re-screen, step 2: decide for every candidate series whether it is an FVB RNA-seq series, from the cached GEO/SRA records
(raw/, rescreen_fetch.py), and re-run the unchanged reference-genome classifier of stage3_survey_v2.py on the RNA-Seq samples only.

Rule (authors' decision, Phase W-1):
  * an FVB RNA-seq series contains >= 1 Mus musculus sample with library_strategy == "RNA-Seq" (SRA run table; the GEO field where the
    SRX is absent from SRA) whose strain is an FVB background: FVB/N, FVB/NJ, FVB/NCrl, FVB/NTac, FVB/NJcl, unspecified "FVB", an
    FVB-congenic designation (FVB.xxx), or a transgenic/knockout line stated to be on an FVB background. A B6-congenic line (B6.FVB-...),
    a mixed background in which FVB is one contributor and an F1/F2/F3 hybrid are not FVB background. miRNA-Seq, ncRNA-Seq, ATAC-seq,
    ChIP-Seq, Bisulfite-Seq, WGS, AMPLICON and OTHER do not count.
  * strain per sample from the characteristics fields naming strain/background/line/genotype; where no sample carries strain
    information, from the series summary/overall design (strain_evidence = series_level).
  * SuperSeries: dropped when one of its SubSeries is in the candidate list (each study counted once, as its SubSeries).
  * judgement cases (strain only in a title or protocol, FVB-derived cell lines or grafts, ...): rescreen_decisions.tsv (gse_acc,
    decision Y/N, proposed_decision, borderline_type, evidence, and for an inclusion that the automatic rule does not make: fvb_gsm =
    the FVB RNA-Seq samples (semicolon list of GSM, or "all_rnaseq"), strain_evidence, strain_text) overrides the automatic decision;
    every row of that file is reported as a borderline series.
  * reference-genome category: classify() of stage3_survey_v2.py (keyword lists unchanged) on the de-duplicated
    !Sample_data_processing text of the series' RNA-Seq mouse samples (of all its RNA-Seq samples if none is FVB, for the record).
Outputs (OUT_DIR): rescreen_series.tsv (all candidates, decision and evidence), survey_of_practice_rescreen.tsv (survivors; Supplementary
Table S5 and Fig. 1A source), excluded_on_rescreen.tsv (Supplementary Table S5b source), survey_classification.tsv (survivors, short
form), summary_stats_rescreen.txt. Titles are the GEO series titles from the cached records (verbatim).
Usage: rescreen_classify.py CANDIDATES_TSV RAW_DIR DECISIONS_TSV OUT_DIR"""
import sys, os, re, csv, collections
from rescreen_parse import series, samples
from stage3_survey_v2 import classify

CAND, RAW, DEC, OUT = sys.argv[1:5]
STRAIN_KEY = re.compile(r"^\s*(strain|background|genetic background|background strain|strain background|mouse strain|mouse line|mouseline|"
                        r"genotype|genotype/variation|line|transgene|donor mouse|mouse)\b[^:]*:", re.I)
OTHER = re.compile(r"c57|\bb6|bl6|bl/6|\b129|balb|\bcd-?1\b|\bsjl|\bdba|swiss|\bc3h|nmri|\bicr\b|\bcast\b|\bmolf|\bnzb|\bsv\b|sv/?129|\bg4\b|\bc67")
HYBRID = re.compile(r"\bf[1-9]\b|hybrid|mixed|\bmix\b|mix of|backcross|\bcross|\bx\b|×|\*")
SC = re.compile(r"single[- ]?(cell|nucleus|nuclei)[- ]?(rna|transcriptom|sequenc|librar|3|5|gene|resolution|atlas)|\bscrna|\bsnrna|\bsc-?rna-?seq|"
                r"10x genomics|chromium|drop-?seq|sci-rna", re.I)   # not "single-cell suspension" (bulk protocols) or "single nucleotide"

def strain_status(text):
    """'fvb' | 'not_fvb' (FVB with another strain, a hybrid or a B6-congenic) | 'other' (another strain only) | '' (no strain named)"""
    t = text.lower()
    if "fvb" not in t: return "other" if OTHER.search(t) else ""
    t = re.sub(r"fvb(/n\w*)?\.[\w/()+-]+", "fvb", t)            # FVB-congenic designations (FVB.B6-Tg(...), FVB.129P2-...) are FVB background
    return "not_fvb" if (OTHER.search(t) or HYBRID.search(t)) else "fvb"

def sample_strain(s):
    vals = [c.split(":", 1)[1].strip() if ":" in c else c for c in s["characteristics"] if STRAIN_KEY.match(c)]
    st = [strain_status(v) for v in vals]
    if "not_fvb" in st: return "not_fvb", vals
    if "fvb" in st: return "fvb", vals
    if "other" in st: return "other", vals
    return "", vals

def series_level(S):
    txt = S["summary"] + " " + S["overall_design"]
    sents = [x for x in re.split(r"(?<=[.;])\s+", txt) if "fvb" in x.lower()]
    st = {strain_status(x) for x in sents}
    return ("not_fvb" if "not_fvb" in st else "fvb" if "fvb" in st else ""), " | ".join(sents)

def main():
    cand = [r for r in csv.DictReader(open(CAND), delimiter="\t")]
    dec = {r["gse_acc"]: r for r in csv.DictReader(open(DEC), delimiter="\t")} if os.path.exists(DEC) else {}
    gset = {r["gse_acc"] for r in cand}; rows = []
    for r in cand:
        g = r["gse_acc"]; S = series(RAW, g); sm = samples(RAW, g)
        mouse = [s for s in sm if s["organism"] and all(o == "Mus musculus" for o in s["organism"])]
        for s in mouse: s["strain_status"], s["strain_vals"] = sample_strain(s)
        rna = [s for s in mouse if s["strategy"] == "RNA-Seq"]
        strategies = sorted({s["strategy"] or ("(none: " + (s["type"] or "?") + ")") for s in sm})
        any_sample_strain = any(s["strain_status"] for s in mouse)
        if any_sample_strain:
            evidence = "sample_level"; fvb_rna = [s for s in rna if s["strain_status"] == "fvb"]
            strain_text = "; ".join(sorted({v for s in (fvb_rna or rna) for v in s["strain_vals"] if v}))[:400]
        else:
            stl, sent = series_level(S); evidence = "series_level"
            fvb_rna = rna if stl == "fvb" else []; strain_text = sent[:400]
        d = dec.get(g)
        if d and d["decision"] == "Y" and d.get("fvb_gsm"):
            fvb_rna = rna if d["fvb_gsm"] == "all_rnaseq" else [x for x in rna if x["gsm"] in d["fvb_gsm"].split(";")]
            evidence, strain_text = d["strain_evidence"], d["strain_text"]
        n_fvb_rna = len(fvb_rna)
        sc_flags = [("SINGLE CELL" in (s["source_lib"] or "").upper()) or bool(SC.search(" ".join([s["title"], s["source"]] + s["description"] + s["extract"])))
                    for s in fvb_rna]
        if not fvb_rna: data_type = ""
        elif all(sc_flags): data_type = "single-cell"
        elif any(sc_flags): data_type = "mixed"
        else: data_type = "single-cell" if SC.search(S["title"] + " " + S["overall_design"]) and not re.search("bulk", S["title"] + S["overall_design"], re.I) else "bulk"
        non_rna_strat = sorted({s["strategy"] for s in sm if s["strategy"] and s["strategy"] != "RNA-Seq"})
        subs_in = [x for x in S["superseries_of"] if x in gset]
        # automatic decision
        if subs_in: auto, reason = "N", "duplicate of SubSeries (" + ", ".join(subs_in) + ")"
        elif not rna: auto, reason = "N", "non-RNA-seq strategy (" + ", ".join(strategies) + ")"
        elif n_fvb_rna == 0:
            auto = "N"
            if any(s["strain_status"] == "not_fvb" for s in rna): reason = "not FVB background (FVB only as part of a mixed, hybrid or B6-congenic background)"
            elif evidence == "series_level" and series_level(S)[0] == "not_fvb": reason = "not FVB background (series description: mixed, hybrid or B6-congenic)"
            elif any(s["strain_status"] == "fvb" for s in mouse): reason = "non-RNA-seq strategy for the FVB samples (" + ", ".join(sorted({s['strategy'] or '-' for s in mouse if s['strain_status']=='fvb'})) + ")"
            else: reason = "not FVB background (no FVB strain in the sample or series records)"
        else: auto, reason = "Y", "FVB RNA-seq"
        final, btype, bev, prop = auto, "", "", ""
        if d:
            final, btype, bev, prop = d["decision"], d["borderline_type"], d["evidence"], d["proposed_decision"]
            if final == "N" and (auto == "Y" or d.get("reason")): reason = d.get("reason") or "excluded on judgement: " + btype
            if final == "Y" and auto == "N": reason = "included on judgement: " + btype
        # reference-genome classification on the RNA-Seq samples (FVB ones when present)
        basis = fvb_rna or rna
        dp = " | ".join(dict.fromkeys(x for s in basis for x in s["data_processing"]))
        cat, kw = classify(dp) if dp else ("not_reported", "")
        t = dp.lower()
        asm = ";".join(a for a, pat in (("GRCm38/mm10", r"grcm38|mm10|ncbi38"), ("GRCm39/mm39", r"grcm39|mm39"), ("mm9", r"\bmm9|ncbi37")) if re.search(pat, t))
        rows.append(dict(gse_acc=g, title=S["title"] or r.get("title", ""), include=final, reason=reason, auto_decision=auto, borderline_type=btype,
            borderline_evidence=bev, proposed_decision=prop, series_type="; ".join(S["types"]), library_strategies=";".join(strategies),
            mixed_strategy="Y" if (rna and non_rna_strat) else "N", data_type=data_type, n_samples=len(sm), n_mouse_rnaseq_samples=len(rna),
            n_fvb_rnaseq_samples=n_fvb_rna, strain_evidence=evidence, strain_text=strain_text, superseries_of=";".join(S["superseries_of"]),
            subseries_of=";".join(S["subseries_of"]), previous_category=r.get("category", ""), category=cat, matched_keyword=kw, assemblies_named=asm,
            evidence_source="GEO_sample_dp" if dp else "none", evidence_text=dp[:600], validation_url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}"))
    os.makedirs(OUT, exist_ok=True); cols = list(rows[0])
    def wr(name, rr, cc):
        with open(f"{OUT}/{name}", "w", newline="") as o:
            w = csv.DictWriter(o, cc, delimiter="\t", extrasaction="ignore", lineterminator="\n"); w.writeheader(); w.writerows(rr)
    wr("rescreen_series.tsv", rows, cols)
    keep = [x for x in rows if x["include"] == "Y"]; drop = [x for x in rows if x["include"] != "Y"]
    wr("survey_of_practice_rescreen.tsv", keep, ["gse_acc", "title", "category", "evidence_source", "matched_keyword", "assemblies_named", "evidence_text",
       "validation_url", "series_type", "library_strategies", "data_type", "n_fvb_rnaseq_samples", "strain_evidence", "strain_text", "borderline_type"])
    wr("survey_classification.tsv", keep, ["gse_acc", "title", "category", "matched_keyword", "assemblies_named", "data_type", "n_fvb_rnaseq_samples", "strain_evidence"])
    wr("excluded_on_rescreen.tsv", drop, ["gse_acc", "title", "reason", "previous_category", "series_type", "library_strategies", "strain_evidence", "strain_text",
       "borderline_type", "validation_url"])
    c = collections.Counter(x["category"] for x in keep); n = len(keep)
    with open(f"{OUT}/summary_stats_rescreen.txt", "w") as o:
        o.write(f"Candidates: {len(rows)}; FVB RNA-seq series after re-screen: {n}; excluded: {len(drop)}\n")
        for k in ("stated_B6", "not_reported", "strain_aware"): o.write(f"{k}: {c.get(k,0)} ({100*c.get(k,0)/n:.1f}%)\n")
    print(open(f"{OUT}/summary_stats_rescreen.txt").read())

if __name__ == "__main__":
    main()
