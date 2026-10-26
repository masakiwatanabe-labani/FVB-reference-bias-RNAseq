#!/usr/bin/env python3
"""Survey screen and classification (Phase W-1 rules; Phase X: broad candidate list, judgement patterns) from the cached GEO/SRA
records (raw/, rescreen_fetch.py); the reference-genome category of the RNA-Seq samples is assigned by classify_reference() of
survey/reference_classifier.py (sentence-level patterns; manual decisions in survey/rescreen_decisions.tsv).

Rule (authors' decision, Phases W-1 and X):
  * an FVB RNA-seq series contains >= 1 Mus musculus sample with library_strategy == "RNA-Seq" (SRA run table; the GEO field where the
    SRX is absent from SRA) whose strain is an FVB background: FVB/N, FVB/NJ, FVB/NCrl, FVB/NTac, FVB/NJcl, unspecified "FVB", an
    FVB-congenic designation (FVB.xxx), or a transgenic/knockout line on an FVB background; FVB-derived cell lines, primary cells and
    FVB-derived tumours grown in FVB hosts count. A B6-congenic line (B6.FVB-...), a mixed background in which FVB is one contributor,
    an F1/F2/F3 hybrid and a multi-strain pool are not FVB background. miRNA-Seq, ncRNA-Seq, ATAC-seq, ChIP-Seq, Bisulfite-Seq, WGS,
    AMPLICON and OTHER do not count.
  * strain per sample from the characteristics fields naming strain/background/line/genotype/cell line; where no sample carries
    strain information: the series summary/overall design (strain_evidence = series_level), then the sample titles/source names
    (sample_title), then the sample growth/treatment protocols (protocol).
  * SuperSeries: dropped when one of its SubSeries is a candidate, unless it holds FVB RNA-Seq samples outside the candidate SubSeries
    (then it represents the study). Sibling SubSeries are counted separately.
Judgement patterns (column `pattern`; proposed decision applied provisionally; each series with a pattern gets a row in
rescreen_decisions.tsv, source = pattern): see PATTERNS. Rows with source = manual (Phase W) override the pattern; a pattern row whose
`decision` differs from its `proposed_decision` (an author's change) is kept and applied on the next run.
Reference-genome category (Phase Z): classify_reference() of reference_classifier.py on the full de-duplicated !Sample_data_processing
text of the FVB RNA-Seq samples (of all RNA-Seq samples of the series if none is FVB, for the record); a manual row may set
`category_override` (with the reason in `evidence`). Phase Z-2: a retained SuperSeries of a retained series is removed at the end
(pairs in OUT_DIR/superseries_pairs_removed.tsv).
Outputs (OUT_DIR): rescreen_series.tsv (all candidates), survey_of_practice_rescreen.tsv (survivors; Supplementary Table S5 and
Fig. 1A source), excluded_on_rescreen.tsv (Supplementary Table S5b source), survey_classification.tsv (survivors, short form),
summary_stats_rescreen.txt; DECISIONS_TSV rewritten (manual rows kept, pattern rows regenerated).
Usage: rescreen_classify.py CANDIDATES_TSV RAW_DIR DECISIONS_TSV OUT_DIR"""
import sys, os, re, csv, collections
from rescreen_parse import series, samples
from reference_classifier import classify_reference   # Phase Z (replaces classify() of stage3_survey_v2.py)

CAND, RAW, DEC, OUT = sys.argv[1:5]
STRAIN_KEY = re.compile(r"^\s*(strain|background|genetic background|background strain|strain background|mouse strain|mouse line|mouseline|"
                        r"genotype|genotype/variation|line|transgene|donor mouse|mouse|mice strain|host strain|cell line)\b[^:]*:", re.I)
OTHER = re.compile(r"c57|\bb6|bl6|bl/6|\b129|balb|\bcd-?1\b|\bsjl|\bdba|swiss|\bc3h|nmri|\bicr\b|(?<![a-z])cast(?![a-wyz])|\bmolf|\bnzb|\bsv\b|sv/?129|\bg4\b|\bc67|\bnod\b|\bnsg\b|\bnude\b|athymic|cast/?ei|\bpwk|\bwsb|\bnzo|\ba/j\b|\bakr|\bcba|\bspret|\bmsm")
HYBRID = re.compile(r"\bf[1-9]\b|hybrid|mixed|\bmix\b|mix of|backcross|\bcross|\bx\b|×|\*|fvb[\w/]*x[a-z(]|[a-z/)]x ?fvb|outcross|crossed (to|with|onto)")
SC = re.compile(r"single[- ]?(cell|nucleus|nuclei)[- ]?(rna|transcriptom|sequenc|librar|3|5|gene|resolution|atlas)|\bscrna|\bsnrna|\bsc-?rna-?seq|"
                r"10x genomics|chromium|drop-?seq|sci-rna", re.I)   # not "single-cell suspension" (bulk protocols) or "single nucleotide"
CELL_LINE = re.compile(r"cell line|\bcells? lines?\b", re.I)
def is_cell_line(s):   # a cell line named in title/source, or a "cell line" field whose value is not a strain name
    if CELL_LINE.search(s["title"] + " " + s["source"]): return True
    return any(c.lower().startswith("cell line") and ":" in c and not strain_status(c.split(":", 1)[1]) and c.split(":", 1)[1].strip().lower() not in ("", "n/a", "na", "none", "-")
               for c in s["characteristics"])
GRAFT = re.compile(r"allograft|syngeneic|orthotopic|transplant|xenograft|engraft|injected (with|into)|inject\w* .{0,40}cells", re.I)
PATTERNS = {   # pattern -> (proposed decision, rule)
 "series_level":        ("Y", "strain only in the series summary/overall design (no sample field names a strain): series text used"),
 "sample_title":        ("Y", "strain only in the sample titles/source names, FVB without another strain: included"),
 "protocol_fvb":        ("Y", "strain only in the sample protocols, FVB without another strain: included"),
 "protocol_ambiguous":  ("N", "strain only in the sample protocols, FVB together with another strain (cross, or 'C57BL/6 or FVB'): excluded"),
 "cell_line":           ("Y", "all FVB RNA-Seq samples are FVB-derived cell lines: included (Phase W decision)"),
 "graft":               ("Y", "FVB samples are grafts/allografts of FVB-derived cells in FVB hosts: included (Phase W decision)"),
 "two_strain_field":    ("N", "a strain field names FVB and another strain without saying F1/mixed (e.g. 'C57Bl6; FVB'): excluded as mixed"),
 "multi_strain_pool":   ("N", "a sample is a pool of several strains (multiplexed embryos/cells): excluded"),
 "conflict":            ("N", "sample fields name another strain while the series text or protocols name FVB: the sample field is followed"),
 "superseries_outside": ("Y", "SuperSeries with FVB RNA-Seq samples outside its candidate SubSeries: kept as the study's record"),
}
def strain_status(text):
    """'fvb' | 'not_fvb' (FVB with another strain, a hybrid or a B6-congenic) | 'other' (another strain only) | '' (no strain named)"""
    t = text.lower()
    if "fvb" not in t: return "other" if OTHER.search(t) else ""
    t = re.sub(r"fvb(/n\w*)?\.[\w/()+-]+", "fvb", t)            # FVB-congenic designations (FVB.B6-Tg(...), FVB.129P2-...) are FVB background
    return "not_fvb" if (OTHER.search(t) or HYBRID.search(t)) else "fvb"
def not_fvb_kind(text):
    t = text.lower()
    if re.search(r"multiplex|\bpool", t) or (len(re.findall(r",", t)) >= 2 and OTHER.search(t) and not HYBRID.search(t)): return "multi-strain pool"
    if re.search(r"\bb6\.|c57bl/?6[a-z]*\.|backcross\w* (on)?to c57|backcross\w* .{0,30}(c57|b6)|onto (a )?c57", t): return "B6-congenic or backcrossed to C57BL/6"
    if re.search(r"\bf[1-9]\b|hybrid|×|\bx\b|\bcross|\*", t): return "F1/F2/F3 hybrid or cross"
    if re.search(r"mix", t) or re.search(r"fvb[^;,]{0,6}[/:+-]|[/:+-][^;,]{0,6}fvb", t): return "mixed background"
    return "two strains named in one field"
def sample_strain(s):
    vals = [c.split(":", 1)[1].strip() if ":" in c else c for c in s["characteristics"] if STRAIN_KEY.match(c)]
    st = [strain_status(v) for v in vals]
    if "not_fvb" in st: return "not_fvb", vals
    if "fvb" in st: return "fvb", vals
    if "other" in st: return "other", vals
    return "", vals
def fvb_sentences(text):
    return [x for x in re.split(r"(?<=[.])\s+(?=[A-Z0-9(])", text) if "fvb" in x.lower()]   # sentence ends only ("FVB; 129S6-..." stays one strain statement)
def series_level(S):
    sents = fvb_sentences(S["summary"] + " " + S["overall_design"]); st = {strain_status(x) for x in sents}
    return ("not_fvb" if "not_fvb" in st else "fvb" if "fvb" in st else ""), " | ".join(sents)

def main():
    cand = [r for r in csv.DictReader(open(CAND), delimiter="\t")]
    dec_rows = list(csv.DictReader(open(DEC), delimiter="\t")) if os.path.exists(DEC) else []
    manual = {r["gse_acc"]: r for r in dec_rows if r.get("source", "manual") == "manual"}
    changed = {r["gse_acc"]: r for r in dec_rows if r.get("source", "manual") != "manual" and r["decision"] != r["proposed_decision"]}
    gset = {r["gse_acc"] for r in cand}; rows = []; pattern_rows = []
    sub_samples = {}
    for r in cand:
        g = r["gse_acc"]; S = series(RAW, g); sm = samples(RAW, g)
        mouse = [s for s in sm if s["organism"] and all(o == "Mus musculus" for o in s["organism"])]
        for s in mouse: s["strain_status"], s["strain_vals"] = sample_strain(s)
        rna = [s for s in mouse if s["strategy"] == "RNA-Seq"]
        strategies = sorted({s["strategy"] or ("(none: " + (s["type"] or "?") + ")") for s in sm})
        any_sample_strain = any(s["strain_status"] for s in mouse)
        pattern, pev, excl_kind = "", "", ""
        if any_sample_strain:
            evidence = "sample_level"; fvb_rna = [s for s in rna if s["strain_status"] == "fvb"]
            strain_text = "; ".join(sorted({v for s in (fvb_rna or rna) for v in s["strain_vals"] if v}))[:400]
            if rna and not fvb_rna:
                nf = [v for s in rna for v in s["strain_vals"] if strain_status(v) == "not_fvb"]
                if nf:
                    excl_kind = not_fvb_kind(nf[0])
                    if excl_kind == "two strains named in one field": pattern, pev = "two_strain_field", nf[0]
                    if excl_kind == "multi-strain pool": pattern, pev = "multi_strain_pool", nf[0]
                else:
                    sl = fvb_sentences(S["summary"] + " " + S["overall_design"] + " " + " ".join(dict.fromkeys(x for s in rna for x in s["protocol"])))
                    if any(strain_status(x) == "fvb" for x in sl) and any(s["strain_status"] == "other" for s in rna):
                        pattern, pev = "conflict", " | ".join(x for x in sl if strain_status(x) == "fvb")[:300]
        else:
            stl, sent = series_level(S)
            if stl:
                evidence = "series_level"; fvb_rna = rna if stl == "fvb" else []; strain_text = sent[:400]
                if stl == "fvb" and rna: pattern, pev = "series_level", sent[:300]
                if stl == "not_fvb": excl_kind = not_fvb_kind(sent)
            else:
                tt = [s for s in rna if strain_status(s["title"] + " " + s["source"]) == "fvb"]
                if tt:
                    evidence = "sample_title"; fvb_rna = tt; strain_text = "; ".join(sorted({s["title"] for s in tt}))[:400]
                    pattern, pev = "sample_title", strain_text[:300]
                else:
                    ps = {s["gsm"]: fvb_sentences(" ".join(s["protocol"])) for s in rna}; allp = [x for v in ps.values() for x in v]
                    evidence = "protocol" if allp else "none"; strain_text = " | ".join(dict.fromkeys(allp))[:400]
                    if allp and all(strain_status(x) == "fvb" for x in allp):
                        fvb_rna = [s for s in rna if ps[s["gsm"]]]; pattern, pev = "protocol_fvb", strain_text[:300]
                    else:
                        fvb_rna = []
                        if allp: pattern, pev = "protocol_ambiguous", strain_text[:300]; excl_kind = not_fvb_kind(" ".join(allp))
        if fvb_rna and not pattern:
            txt = lambda s: " ".join([s["title"], s["source"]] + s["characteristics"] + s["description"])
            if all(is_cell_line(s) for s in fvb_rna): pattern, pev = "cell_line", "; ".join(sorted({s["source"] for s in fvb_rna}))[:300]
            elif all(GRAFT.search(txt(s)) for s in fvb_rna): pattern, pev = "graft", "; ".join(sorted({s["source"] for s in fvb_rna}))[:300]
        d = manual.get(g)
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
        if subs_in and fvb_rna:   # FVB RNA-Seq samples of the SuperSeries not in any candidate SubSeries
            inside = set().union(*[set(sub_samples.setdefault(x, series(RAW, x)["sample_ids"])) for x in subs_in])
            outside = [s for s in fvb_rna if s["gsm"] not in inside]
            if outside: pattern, pev = "superseries_outside", ";".join(s["gsm"] for s in outside)[:300]; fvb_rna_out = outside
        # automatic decision
        if subs_in and pattern != "superseries_outside": auto, reason = "N", "duplicate of SubSeries (" + ", ".join(subs_in) + ")"
        elif not rna: auto, reason = "N", "non-RNA-seq strategy (" + ", ".join(strategies) + ")"
        elif n_fvb_rna == 0:
            auto = "N"
            if excl_kind: reason = f"not FVB background ({excl_kind})"
            elif any(s["strain_status"] == "fvb" for s in mouse): reason = "non-RNA-seq strategy for the FVB samples (" + ", ".join(sorted({s['strategy'] or '-' for s in mouse if s['strain_status']=='fvb'})) + ")"
            elif evidence == "none" or not any(s["strain_status"] for s in rna): reason = "not FVB background (no FVB strain in the sample or series records)"
            else: reason = "not FVB background (other strains only)"
        else: auto, reason = "Y", "FVB RNA-seq"
        if pattern in ("conflict", "two_strain_field", "protocol_ambiguous", "multi_strain_pool"): auto = "N"
        final, btype, bev, prop, source = auto, "", "", "", ""
        if pattern:
            prop = PATTERNS[pattern][0]; final = prop; btype = pattern; bev = pev; source = "pattern"
            c = changed.get(g)
            if c: final = c["decision"]
            if final == "Y" and auto == "N" and pattern != "superseries_outside": reason = f"included by pattern {pattern}"
            if final == "N" and reason == "FVB RNA-seq": reason = f"not FVB background (pattern {pattern})"
            if final == "Y" and pattern == "superseries_outside": reason = "FVB RNA-seq (SuperSeries; FVB samples outside the candidate SubSeries)"
            pattern_rows.append(dict(gse_acc=g, decision=final, proposed_decision=prop, borderline_type=PATTERNS[pattern][1], evidence=pev, reason="" if final == "Y" else reason,
                                     fvb_gsm="", strain_evidence=evidence, strain_text=strain_text[:200], pattern=pattern, source="pattern"))
        if d:
            final, btype, bev, prop, source = d["decision"], d["borderline_type"], d["evidence"], d["proposed_decision"], "manual"
            if final == "N": reason = d.get("reason") or "excluded on judgement: " + btype
            if final == "Y" and (auto == "N" or reason.startswith("not FVB")): reason = "included on judgement: " + btype
        if final == "N" and reason == "FVB RNA-seq": reason = "excluded on judgement"
        basis = fvb_rna or rna
        dpl = list(dict.fromkeys(x for s in basis for x in s["data_processing"])); dp = " | ".join(dpl)
        cr = classify_reference(dpl); cat, kw, asm = cr["category"], cr["matched_keyword"], cr["assemblies_named"]
        if d and d.get("category_override"):   # manual verification overrule (Phase Z-1.5)
            cat, kw = d["category_override"], "manual: " + d["borderline_type"]
            if cat != "strain_aware": cr["strain_aware_method"] = ""
        rows.append(dict(gse_acc=g, title=S["title"] or r.get("title", ""), include=final, reason=reason, auto_decision=auto, pattern=pattern if source == "pattern" else (d.get("pattern", "") if d else ""),
            decision_source=source, borderline_type=btype, borderline_evidence=bev, proposed_decision=prop, series_type="; ".join(S["types"]),
            library_strategies=";".join(strategies), mixed_strategy="Y" if (rna and non_rna_strat) else "N", data_type=data_type, n_samples=len(sm),
            n_mouse_rnaseq_samples=len(rna), n_fvb_rnaseq_samples=n_fvb_rna, strain_evidence=evidence, strain_text=strain_text, superseries_of=";".join(S["superseries_of"]),
            subseries_of=";".join(S["subseries_of"]), source_query=r.get("source_query", ""), previous_category=r.get("category", ""), category=cat, matched_keyword=kw,
            assemblies_named=asm, strain_aware_method=cr["strain_aware_method"], reference_named_without_assembly=cr["reference_named_without_assembly"], human_assembly_named=cr["human_assembly_named"], assembly_name_ambiguous=cr["assembly_name_ambiguous"],
            decisive_sentence=cr["sentence"][:500], release_date=r.get("release_date", ""), evidence_source="GEO_sample_dp" if dp else "none", evidence_text=dp[:600], validation_url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}"))
    # Phase Z-2: a retained series that is the SuperSeries of another retained series is a duplicate (relations read from both sides)
    kept = {x["gse_acc"] for x in rows if x["include"] == "Y"}; byg = {x["gse_acc"]: x for x in rows}; pairs = []
    for x in rows:
        if x["include"] != "Y": continue
        subs = set(filter(None, x["superseries_of"].split(";"))) | {y["gse_acc"] for y in rows if x["gse_acc"] in y["subseries_of"].split(";")}
        live = sorted(subs & kept - {x["gse_acc"]})
        if live:
            pairs.append((x["gse_acc"], live)); x["include"] = "N"; x["reason"] = "duplicate of SubSeries (" + ", ".join(live) + ", retained)"
            for p_ in pattern_rows:
                if p_["gse_acc"] == x["gse_acc"]: p_["decision"] = p_["proposed_decision"] = "N"; p_["reason"] = x["reason"]   # rule outcome, not an author change
    with open(f"{OUT}/superseries_pairs_removed.tsv", "w") as o:
        o.write("superseries\tretained_subseries\n"); [o.write(f"{a}\t{';'.join(b)}\n") for a, b in pairs]
    os.makedirs(OUT, exist_ok=True); cols = list(rows[0])
    def wr(path, rr, cc):
        with open(path, "w", newline="") as o:
            w = csv.DictWriter(o, cc, delimiter="\t", extrasaction="ignore", lineterminator="\n"); w.writeheader(); w.writerows(rr)
    wr(f"{OUT}/rescreen_series.tsv", rows, cols)
    keep = [x for x in rows if x["include"] == "Y"]; drop = [x for x in rows if x["include"] != "Y"]
    wr(f"{OUT}/survey_of_practice_rescreen.tsv", keep, ["gse_acc", "title", "category", "evidence_source", "matched_keyword", "assemblies_named", "evidence_text",
       "validation_url", "series_type", "library_strategies", "data_type", "n_fvb_rnaseq_samples", "strain_evidence", "strain_text", "borderline_type",
       "strain_aware_method", "reference_named_without_assembly", "human_assembly_named", "assembly_name_ambiguous", "decisive_sentence", "release_date"])
    wr(f"{OUT}/survey_classification.tsv", keep, ["gse_acc", "title", "category", "matched_keyword", "assemblies_named", "data_type", "n_fvb_rnaseq_samples", "strain_evidence"])
    wr(f"{OUT}/excluded_on_rescreen.tsv", drop, ["gse_acc", "title", "reason", "previous_category", "series_type", "library_strategies", "strain_evidence", "strain_text",
       "borderline_type", "validation_url"])
    dcols = ["gse_acc", "decision", "proposed_decision", "borderline_type", "evidence", "reason", "fvb_gsm", "strain_evidence", "strain_text", "pattern", "source", "category_override"]
    man_rows = [dict(r, source="manual", pattern=r.get("pattern") or "manual") for r in dec_rows if r.get("source", "manual") == "manual"]
    wr(DEC, man_rows + [p for p in pattern_rows if p["gse_acc"] not in manual], dcols)
    c = collections.Counter(x["category"] for x in keep); n = len(keep)
    with open(f"{OUT}/summary_stats_rescreen.txt", "w") as o:
        o.write(f"Candidates: {len(rows)}; FVB RNA-seq series after screening: {n}; excluded: {len(drop)}\n")
        for k in ("stated_B6", "not_reported", "strain_aware"): o.write(f"{k}: {c.get(k,0)} ({100*c.get(k,0)/n:.1f}%)\n")
    print(open(f"{OUT}/summary_stats_rescreen.txt").read())

if __name__ == "__main__":
    main()
