# Survey of reference-genome practice in public FVB RNA-seq studies (Fig. 1A; Supplementary Table S5)

## Files
- `candidates_154.tsv` — the 154 candidate GEO series and their first classification (June 2026): `gse_acc`, `title`, `category`,
  `matched_keyword`, `matched_field`. The list is the output of the GEO search script of the project (three E-utilities queries that
  required both FVB and C57BL/6, at most 100 records each) plus six series added by hand; the titles of those six were typed in by the
  project, not taken from GEO. No other pre-review list exists.
- `raw/` — the GEO and SRA records of every candidate, fetched on 2026-10-26 (`FETCHED.txt`): `<GSE>.series.soft.gz` (series record),
  `<GSE>.samples.soft.gz` (all sample records: characteristics, library strategy and source, data processing, SRA/BioSample links) and
  `<GSE>.sra_runinfo.csv.gz` (SRA run table of the linked SRX experiments). The submitters' contact lines (`!Series_contact_*`,
  `!Sample_contact_*`) are removed; nothing else is changed.
- `rescreen_fetch.py` — writes `raw/` (usage: `python3 rescreen_fetch.py candidates_154.tsv raw`).
- `rescreen_parse.py` — parsers for the cached records.
- `rescreen_classify.py` — the re-screen and the classification (usage: `python3 rescreen_classify.py candidates_154.tsv raw
  rescreen_decisions.tsv OUT_DIR`). A series is kept when at least one *Mus musculus* sample has library strategy RNA-Seq (SRA run
  table; the GEO field where the SRX is not in SRA) and an FVB background (FVB/N, FVB/NJ, FVB/NCrl, FVB/NTac, FVB/NJcl, "FVB",
  FVB-congenic designations, transgenic or knockout lines stated to be on FVB); a B6-congenic line, a mixed background with FVB as one
  contributor and F1/F2/F3 hybrids are not FVB background. The strain is read per sample from the characteristics fields naming
  strain, background, line or genotype, and from the series summary/overall design only where no sample carries strain information.
  A SuperSeries is dropped when one of its SubSeries is a candidate. The reference-genome category is assigned by `classify()` of
  `stage3_survey_v2.py` (keyword lists unchanged) on the sample-level data-processing text of the FVB RNA-Seq samples.
- `rescreen_decisions.tsv` — the series whose decision needed judgement (strain only in a title or protocol, FVB-derived cell lines
  and grafts, strain fields naming two strains, SuperSeries with an FVB sample outside the listed SubSeries, conflicting records), with
  the evidence and the proposed decision; `rescreen_classify.py` applies it.
- `rescreen_series.tsv` — every candidate with the decision, the reason, the evidence columns (series type, library strategies,
  data type, number of FVB RNA-Seq samples, strain evidence and text) and the category.
- `survey_classification.tsv` — the series kept (Supplementary Table S5 has the full columns; the excluded ones are in sheet S5b).
- `summary_stats_rescreen.txt` — counts.
- `stage3_survey_v2.py` — the original survey script (online fetch; `classify()` and the keyword lists `STRAIN_AWARE_KW`, `B6_KW` are
  used unchanged by the re-screen).

## Result (2026-10-26)
See `summary_stats_rescreen.txt` and `rescreen_series.tsv`. The re-run classification uses the GEO records as fetched on 2026-10-26,
which the submitters may have edited since June 2026.
