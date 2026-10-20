# Phase A audit summary — all 83 SRRs of the old analysis

Sources: `revision_v2/audit/sample_audit.tsv` (per run), `old_runs_metadata.tsv`, `series_structure.tsv`, raw ENA/GEO records in `revision_v2/audit/raw/` (2026-10-13), sequence QC in `analysis_v2/qc/` (`qc_calibration.md` for the rules). All 83 BAMs were scored; 0 failures.

Overall: of the 71 pure-strain runs, **0 match the manuscript description** on all of strain, tissue, biological n and QC. The 12 F1 runs match. Breakdown of QC verdicts across 83 runs: FVB/NJ pure 12, B6 pure 31, F1 32, other inbred strain 6, needs review 2.

## GSE123875 — mesenteric adipose (fructose study, SubSeries of GSE123896)
- Manuscript: FVB/NJ n = 6 vs B6 n = 6, paired-end.
- Reality: the six "FVB" runs are **two DBA/2J animals** (GSM3514981 DWM6, GSM3514982 DWM8; water group; 3 technical runs each). The six "B6" runs are **two C57BL/6J animals** (GSM3514966 CFM15, GSM3514967 CFM9; **fructose** group; 3 runs each).
- Sequence QC: "FVB" arm = other inbred strain (mean alt fraction 0.55–0.60, bimodal); "B6" arm = B6 pure (0.002–0.003).
- Mismatch type: wrong strain (DBA/2J ≠ FVB/NJ), n counted as runs (2 vs 2 animals, not 6 vs 6), diet confounded with "strain" (water vs fructose).
- Impact: the adipose arm of Figs. 3–6, Table 1, S2–S4 and the catalog; the Discussion paragraph explaining the low FVB mapping rate in adipose (47.4 % vs 63.5 %) was explaining a DBA/2J-vs-B6 difference. The series does contain FVB/NJ adipose samples (water 3, fructose 5 animals) and C57BL/6J (water 4, fructose 2), which were not used.

## GSE123893 — described as liver; is hypothalamus
- Manuscript: liver, FVB/NJ n = 6 vs B6 n = 6.
- Reality: **hypothalamus** (series title "[Hypothalamus]", source "Hypothalamus tissue"). FVB arm = 3 FVB/NJ animals (FFHyp6 fructose ×1 run; FWHyp1 water ×3; FWHyp2 water ×2). B6 arm = 2 C57BL/6J animals (CFHyp10, CFHyp11; fructose; 3 runs each).
- Sequence QC: FVB/NJ pure ×6 (0.992–0.993); B6 pure ×6 (0.001–0.002).
- Mismatch type: wrong tissue; n = runs not animals (3 vs 2); diet unbalanced.
- Impact: everything labelled "liver" in the old manuscript refers to hypothalamus. Strain assignment itself was correct.

## GSE123894 — described as kidney; is liver
- Manuscript: kidney, FVB/NJ n = 6 vs B6 n = 6.
- Reality: **liver** (series "[Liver]"). FVB arm = 3 FVB/NJ animals (FFL10 ×1, FFL12 ×3, FFL18 ×2; all fructose). B6 arm = 2 C57BL/6J animals (CFL10, CFL11; fructose; 3 runs each). No kidney data exist in this study.
- Sequence QC: FVB/NJ pure ×6 (0.983–0.988); B6 pure ×6 (0.003–0.004).
- Mismatch type: wrong tissue; n = runs (3 vs 2). Strain correct; both arms fructose-fed (balanced).

## GSE135230 — described as cardiac muscle FVB/N vs B6; is Bin1 cKO brain microglia on one mixed background
- Manuscript: cardiac muscle, FVB/N n = 10 vs B6 n = 12, single-end.
- Reality: CD11b+ CD45med FACS-sorted **microglia** from 11 mice, all "B6.129S6-Bin1tm2Gcp/J; <FVB/N-Tg(Thy1-cre)1Vln/J>", 2 technical runs each. "FVB" arm = 5 Het-cKO mice; "B6" arm = 3 Het-cKO + 3 Homo-cKO mice. The "FVB/N" string in the manuscript came from the Thy1-cre transgenic line name.
- Sequence QC: Het-cKO (16 runs) = F1-like (alt fraction 0.44–0.46, 93 % heterozygous sites); Homo-cKO (6 runs) = backcross-like (0.24–0.26; 50 % het, 43–49 % B6-homozygous) — 4 "F1", 2 "needs review" under the fixed rules (see qc_calibration.md).
- Mismatch type: wrong tissue, wrong strains, no FVB-vs-B6 contrast exists; the two arms differ by Bin1 genotype and generation.
- Impact: the "cardiac" tissue in all figures/tables, including the Fau worked example (Fig. 6, "cardiac +3.67") and the catalog criterion "≥ 2 of 4 non-PyVT tissues". Unusable for any purpose in the revision.

## GSE175625 — described as PyVT FVB/N vs B6; is B6-congenic PyVT tumours with paraben treatments
- Manuscript: mammary tumour PyVT, FVB/N n = 7 vs B6 n = 6, single-end.
- Reality: 13 mammary tumours from **B6.FVB-Tg(MMTV-PyVT)634Mul/LellJ** mice (1 run each). "FVB" arm = tumours 29–55 (vehicle 2, methylparaben 3, propylparaben 2); "B6" arm = tumours 60–105 (vehicle 1, propylparaben 4, methylparaben 1).
- Sequence QC: **B6 pure ×13** (alt fraction 0.011–0.016).
- Mismatch type: wrong strain for the "FVB" arm; arms differ only by animal ID / treatment.
- Impact: PyVT column of Figs. 3–5, Table 1, S2–S4, the "pvt_discordant" catalog flag and the Discussion CNV paragraph. Candidate negative control in Phase B (all samples carry essentially no FVB alleles).

## GSE200632 — F1 keratinocytes (OK)
- Manuscript: four C57BL/6 × FVB/N F1, primary keratinocytes; F3 backcross excluded.
- Reality: GSM6040397–400 "KrasC57_replicate1–4", background strain "C57BL/6:FVB/N (F1:50/50)", LSL-KrasG12D, 1 run each; the four F3 samples were not used. ✔
- Sequence QC: F1 ×4 (0.474–0.477; 96 % heterozygous).

## GSE293440 — F1 cortex (OK, two wording points)
- Manuscript: eight C57BL/6J × FVB/NJ F1, frontal cortex.
- Reality: one GEO series of 290 samples (no SubSeries; BioProject PRJNA1243442; ENA SRP573960). The eight runs used are exactly the eight "B6xFVB, Control, WT" samples (4 young / 4 old; 4 M / 4 F), 1 run each. ✔
- Wording: GEO sample attribute tissue = "Cerebral cortex" (series overall design says frontal cortex); BioSample strain = "B6xFVB" — the FVB **substrain is not stated** in GEO, so "FVB/NJ" and "substrain-matched" are inferences, not metadata.
- Sequence QC: F1 ×8 (0.474–0.481; 98 % heterozygous).

## n in the old manuscript
Stated n (6/6, 6/6, 6/6, 10/12, 7/6) are run counts. Biological n: adipose 2 vs 2; hypothalamus 3 vs 2; liver 3 vs 2; microglia 5 vs 6 (same background); PyVT 7 vs 6 (same strain). Technical runs per BioSample: 1–3.

## Methods verification sentences
"confirmed from GEO series metadata" and "recorded in GEO sample attributes" are contradicted by the attributes retrieved here; no verification had been performed. These sentences are withdrawn in the revision and replaced by the Sample identity verification subsection backed by `sample_audit.tsv`.

## Samples available for Phase B (fructose study, GEO attributes; QC after download)
| Tissue (GSE) | FVB/NJ water | FVB/NJ fructose | C57BL/6J water | C57BL/6J fructose |
|---|---|---|---|---|
| Adipose (GSE123875) | 3 | 5 | 4 | 2 |
| Hypothalamus (GSE123893) | 4 | 4 | 4 | 4 |
| Liver (GSE123894) | 6 | 6 | 6 | 6 |
