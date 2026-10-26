#!/usr/bin/env python3
"""Phase W (2026-10-26): text items for the authors, appended to revision_v2/change_ledger.tsv.
The manuscript, Additional file 1 and the letters are maintained by the authors outside the repository, so the ledger cannot be
regenerated from the section files (make_change_ledger.py). This script finds, in the last plain-text dumps of the assembled documents
(revision_v2/submission/LAR_main_text_v2.plain.txt and LAR_additional_file_1.plain.txt, 2026-10-24), every sentence that a Phase W
correction affects, and writes one ledger row per sentence: section = document, status = REWRITE (the authors rewrite it),
trigger = W-1/W-2/W-3, new_text = the current sentence, number_check = what changes (old -> new, numbers.tsv key). Rows of earlier
runs with the same trigger are replaced; all other rows are kept. Usage: phase_w_ledger.py W-1 [W-2 W-3]  (the items to write)"""
import csv, re, sys
P = "/usr/local/jupyter/FVB_bias"; LED = f"{P}/revision_v2/change_ledger.tsv"
DOCS = {"main text": f"{P}/revision_v2/submission/LAR_main_text_v2.plain.txt", "Additional file 1": f"{P}/revision_v2/submission/LAR_additional_file_1.plain.txt"}
N = {r["key"]: r["value"] for r in csv.DictReader(open(f"{P}/manuscript/revised_v2/generated/numbers.tsv"), delimiter="\t")}
def n(k): return N.get(k, "?")
ITEMS = {   # trigger -> [(regex on the sentence, what changes)]
 "W-1": [
  (r"returned 1,047 GEO series", "search description: no record of a 1,047-series search; the candidate list (154) is the output of scripts/create_survey_tsv.py (three FVB AND C57BL/6 queries, retmax 100 each, plus 6 series added by hand)"),
  (r"FVB was not the main experimental strain|transgenic host background or as a cell-line origin|fewer than three samples", "criteria not applied to the 154 candidates; superseded by the Phase W-1 rule (RNA-Seq library strategy, FVB background per sample)"),
  (r"After manual review, 154", f"154 -> {n('survey_total')} (survey_total); candidates 154, excluded on re-screen {154 - int(n('survey_total')) if n('survey_total') != '?' else '?'}"),
  (r"154 GEO series|154 public|154 publicly|154 series|these 154|of 154", f"154 -> {n('survey_total')} (survey_total)"),
  (r"\b140\b", f"140 -> {n('survey_stated_B6_n')} (survey_stated_B6_n); assemblies: GRCm38/mm10 {n('survey_b6_grcm38_mm10_n')}, mm9 {n('survey_b6_mm9_n')}, GRCm39/mm39 {n('survey_b6_grcm39_mm39_n')}"),
  (r"\b14 (studies|series|did)|remaining 14", f"14 -> {n('survey_not_reported_n')} (survey_not_reported_n)"),
  (r"GRCm38/mm10 or GRCm39/mm39|legacy mm9|\bmm9\b", f"assemblies among stated B6: GRCm38/mm10 {n('survey_b6_grcm38_mm10_n')}, mm9 {n('survey_b6_mm9_n')}, GRCm39/mm39 {n('survey_b6_grcm39_mm39_n')} (survey_b6_*_n)"),
  (r"Studies were retrieved from NCBI GEO", "S5 caption: search terms as in the main text; add the re-screen rule and the new columns (series_type, library_strategies, data_type, n_fvb_rnaseq_samples, strain_evidence, strain_text, assemblies_named, borderline_type) and sheet S5b"),
 ],
 "W-2": [
  (r"Thirty-seven of these 82|37 genes had residual fractions", f"37 -> {n('mech_multi_absresidual_le_0.30')} with |bias_multi|/|bias_unique| <= 0.30 (mech_multi_absresidual_le_0.30); the signed ratio also admits the {n('mech_multi_signed_not_abs_le_0.30_n')} sign-reversing genes {n('mech_multi_signed_not_abs_le_0.30_genes')}"),
  (r"more than 70% of", "the 70% reading holds only for |residual| <= 0.30 (abs_residual_fraction, Supplementary Table S2)"),
  (r"Negative values arise when|reverses sign or exceeds", f"{n('mech_multi_sign_reversed_n')} of 82 genes change sign: {n('mech_multi_sign_reversed_genes')} (mech_multi_sign_reversed_*); |residual| range {n('mech_multi_absresidual_min')}–{n('mech_multi_absresidual_max')}"),
  (r"Lower values indicate that more of the unique-only bias|value near 0 indicates bias", "true for 0 <= residual <= 1 only; negative values are sign reversals (abs_residual_fraction, sign_reversed in Supplementary Table S2)"),
  (r"residual_fraction, mechanism_class|residual_fraction, mean_bias_multi / mean_bias", "S2 column list: add abs_residual_fraction (|mean_bias_multi| / |mean_bias|) and sign_reversed (Y/N) after residual_fraction"),
 ],
 "W-3": [
  (r"36 FVB/NJ and 34 C57BL/6J animals", f"samples, not animals: {n('animals_fvb_samples')} FVB/NJ and {n('animals_b6_samples')} C57BL/6J samples from {n('animals_fvb_n')} FVB/NJ and {n('animals_b6_n')} C57BL/6J animals (animals_*)"),
  (r"n denotes biological animals", "n per tissue counts samples; adipose, hypothalamus and liver of the fructose study share animals (S1c)"),
  (r"83 animals", f"83 samples from {n('animals_threearm_n')} animals ({n('animals_strain_comparison_n')} FVB/NJ + C57BL/6J animals and 13 PyVT mice) (animals_threearm_n)"),
  (r"12 F1 animals|the 12 animals|over the 12 animals", f"12 F1 samples: 8 GSE293440 animals and 4 GSE200632 keratinocyte cultures (donor animals not recorded)"),
  (r"all four tissues and in both laboratories|all four tissues and both laboratories", f"three of the four tissues (adipose, hypothalamus, liver) come from the same fructose-study animals (animals_shared_*); tissues are not independent at the animal level"),
  (r"at least two of the four|≥ 2 of 4|in ≥ 2 tissues|in at least two", "catalog criterion counts tissues; two qualifying tissues may come from the same animals (fructose study); see animals_pairs_* in numbers.tsv"),
  (r"independent laboratory|two laboratories", "the fructose-study tissues are one laboratory and shared animals; skin is the only animal-independent tissue"),
  (r"One row per animal × arm|per animal\b|negative control per animal", "rows are samples (sample × arm); 83 samples from 59 animals"),
  (r"Sheet S1a — manifest", "S1a gains animal_id and n_tissues_from_animal; new sheet S1c — animals (animal_id, strain, diet_or_condition, dataset, tissues_sampled, n_samples, role; one row per animal with an included sample)"),
 ],
}
def sentences(t):
    t = re.sub(r"[ \t]+", " ", t)
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(“0-9])|\n\s*\n", t) if s.strip()]
want = sys.argv[1:] or list(ITEMS)
rows = [r for r in csv.DictReader(open(LED), delimiter="\t") if r["trigger"] not in want]
cols = list(csv.DictReader(open(LED), delimiter="\t").fieldnames)
new = []
for trig in want:
    for doc, path in DOCS.items():
        S = sentences(open(path, encoding="utf-8").read())
        for i, s in enumerate(S):
            for pat, what in ITEMS[trig]:
                if re.search(pat, s):
                    new.append(dict(section=doc, para_no="", sent_no=str(i + 1), status="REWRITE", trigger=trig, v1_text="", new_text=s[:1500],
                                    v1_section="", v1_subsection="", v1_para="", v1_sub_para="", number_check=what)); break
with open(LED, "w", newline="") as o:
    w = csv.DictWriter(o, cols, delimiter="\t"); w.writeheader(); w.writerows(rows + new)   # csv default line terminator (CRLF), as make_change_ledger.py writes
import collections
print(dict(collections.Counter((r["trigger"], r["section"]) for r in new)))
