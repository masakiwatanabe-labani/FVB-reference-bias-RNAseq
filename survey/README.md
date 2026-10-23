# Survey of reference-genome practice in public FVB RNA-seq studies (Fig. 1A; Supplementary Table S5)

- `stage3_survey_v2.py` — the survey script. For each GEO series it fetches the GEO SOFT text (sample level first,
  `!Sample_data_processing`; series level, `!Series_data_processing`, only if no sample-level text exists) and classifies the
  series by keyword matching in priority order strain_aware > stated_B6 > not_reported (function `classify`). The complete
  keyword lists are `STRAIN_AWARE_KW` and `B6_KW` in the script (case-insensitive substring matching).
  Usage: `python3 stage3_survey_v2.py [INPUT_TSV] [OUTPUT_DIR]`; the input needs a column `gse_acc`
  (`survey_classification.tsv` can be used as the input list).
- `survey_classification.tsv` — the classification of the 154 series: `gse_acc`, `title`, `category`, `matched_keyword`
  (the first keyword of the lists that occurred in the text; empty for not_reported), `matched_field` (the GEO field the keyword
  was found in; for not_reported, the field that was searched). Result: 140 stated_B6, 14 not_reported, 0 strain_aware;
  all 154 series had sample-level data-processing text.

The classification was made on the full GEO data-processing text fetched in June 2026. That text was not cached, so it is not
included here; re-running the script fetches the current GEO records, which may have been edited since. The `evidence_text`
column of Supplementary Table S5 holds only the first 600 characters of that text and is not the classifier input (for 19 series
the matched keyword lies beyond the first 600 characters).
