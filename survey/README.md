# Survey of reference-genome practice in public FVB RNA-seq studies (Fig. 1A; Supplementary Table S5)

## Files
- `broad_search.py` — the candidate list (Phase X, 2026-10-26): E-utilities searches of GEO DataSets (`db=gds`), with the exact
  query strings, date and hit counts in `search_log.tsv`. Series level: `(FVB[All Fields]) AND "Mus musculus"[Organism] AND
  gse[Entry Type]` with DataSet Type "Expression profiling by high throughput sequencing", "Other" and "Non-coding RNA profiling by
  high throughput sequencing"; sample level: `(FVB[All Fields]) AND "Mus musculus"[Organism] AND gsm[Entry Type]`, each GSM mapped to
  its GSE(s), a GSE kept when its type includes high-throughput sequencing. Union → `candidates_broad.tsv` (gse_acc, source_query,
  series_type, title, release_date, n_samples). Usage: `python3 broad_search.py OUT_DIR [--variants]` (`--variants` tests FVB/N,
  FVB/NJ, FVB/NCrl, FVB/NTac, FVBN, FVB-N in place of FVB and reports the accessions they add).
- `candidates_154.tsv` — the earlier 154-series candidate list and its first classification (June 2026), kept for comparison: the
  output of a narrower search that required both FVB and C57BL/6, plus six series added by hand (titles typed by the project).
- `raw/` — GEO and SRA records of every candidate (`fetch.log`, `FETCHED.txt`): `<GSE>.series.soft.gz`, `<GSE>.samples.soft.gz`
  (characteristics, library strategy and source, protocols, data processing, SRA/BioSample links) and `<GSE>.sra_runinfo.csv.gz`.
  The submitters' contact lines (`!Series_contact_*`, `!Sample_contact_*`) are removed; nothing else is changed.
- `rescreen_fetch.py` — writes `raw/` (usage: `python3 rescreen_fetch.py candidates_broad.tsv raw`).
- `rescreen_parse.py` — parsers for the cached records.
- `rescreen_classify.py` — the screen and the classification (usage: `python3 rescreen_classify.py candidates_broad.tsv raw
  rescreen_decisions.tsv OUT_DIR`). A series is kept when at least one *Mus musculus* sample has library strategy RNA-Seq and an FVB
  background (FVB/N, FVB/NJ, FVB/NCrl, FVB/NTac, FVB/NJcl, "FVB", FVB-congenic designations, transgenic or knockout lines on FVB;
  FVB-derived cell lines, primary cells and tumours grown in FVB hosts); B6-congenic lines, mixed backgrounds with FVB as one
  contributor, F1/F2/F3 hybrids and multi-strain pools are not FVB background. Strain per sample from the characteristics; where no
  sample names a strain, from the series summary/design, then the sample titles, then the sample protocols. A SuperSeries is dropped
  when one of its SubSeries is a candidate, unless it holds FVB RNA-Seq samples outside them. The reference-genome category is
  assigned by `classify()` of `stage3_survey_v2.py` (keyword lists unchanged) on the data-processing text of the FVB RNA-Seq samples.
- `rescreen_decisions.tsv` — every series whose decision needs judgement: `source = manual` (the 16 Phase W decisions, applied as
  written) and `source = pattern` (one row per series matched by a judgement pattern, `pattern` column; the proposed decision of the
  pattern is applied provisionally; a row whose `decision` is changed from `proposed_decision` is kept and applied on the next run).
- `rescreen_series.tsv` — every candidate with the decision, reason, pattern and evidence columns and the category.
- `survey_classification.tsv` — the series kept (Supplementary Table S5 has the full columns; the excluded ones are in sheet S5b).
- `summary_stats_rescreen.txt` — counts.
- `stage3_survey_v2.py` — the original survey script; `classify()` and its keyword lists are used unchanged.

The classification uses the GEO records as fetched on 2026-10-26, which the submitters may have edited since earlier fetches.
