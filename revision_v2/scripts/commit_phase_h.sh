#!/bin/bash
# Phase H-8 (rule 0.4-9: preserve history, no force-push, no push from here). Copies the Phase H changes into github_repo on branch
# v2-corrected (same file set as prepare_github_v2.sh, plus the Phase H audit records and the submission sources) and commits them
# as a NEW commit on top of the Phase F commits. The local tag "v2-corrected" created in Phase F pointed at the Phase F commit and
# was never pushed; it is removed so that the authors create the submission tag at the final commit (the instruction assigns push
# and tag to the authors). Usage: commit_phase_h.sh
set -euo pipefail
PROJ=/usr/local/jupyter/FVB_bias; R=$PROJ/github_repo; cd $R
git checkout -q v2-corrected
mkdir -p analysis_v2/pipeline revision_v2/scripts revision_v2/manifest revision_v2/audit/f1_old revision_v2/submission analysis_v2/results/qc analysis_v2/results/catalog analysis_v2/results/negative_control analysis_v2/results/sensitivity_matched/bias analysis_v2/results/predictors analysis_v2/results/within_fvb analysis_v2/results/f1_ase/v2 analysis_v2/figures manuscript/revised_v2
cp $PROJ/analysis_v2/pipeline/*.py $PROJ/analysis_v2/pipeline/*.sh $PROJ/analysis_v2/pipeline/*.R analysis_v2/pipeline/; cp $PROJ/analysis_v2/params.yaml analysis_v2/
cp $PROJ/revision_v2/scripts/*.py $PROJ/revision_v2/scripts/*.sh revision_v2/scripts/; cp $PROJ/revision/scripts/check_figure_order.py revision_v2/scripts/
cp $PROJ/revision_v2/manifest/manifest.tsv $PROJ/revision_v2/manifest/design_proposal.md revision_v2/manifest/
cp $PROJ/revision_v2/CHANGELOG.md $PROJ/revision_v2/text_changes.md $PROJ/revision_v2/summary_of_changes.md $PROJ/revision_v2/final_check.md $PROJ/revision_v2/change_ledger.tsv $PROJ/revision_v2/ledger_paragraph_triggers.tsv $PROJ/revision_v2/response_to_reviewers.md $PROJ/revision_v2/cover_letter.md revision_v2/
cp $PROJ/revision_v2/audit/sample_audit.tsv $PROJ/revision_v2/audit/audit_summary.md $PROJ/revision_v2/audit/gse55792_qc_note.md $PROJ/revision_v2/audit/qc_calibration.md $PROJ/revision_v2/audit/junction_sets_fvb_index_vs_b6_gtf.txt $PROJ/revision_v2/audit/fvb_vcf_counts.txt revision_v2/audit/ 2>/dev/null || true
cp $PROJ/revision_v2/audit/f1_old_site_genotypes.md $PROJ/revision_v2/audit/rps2_indels_pseudogenes.md $PROJ/revision_v2/audit/rps2_pseudogene_counts_skin.tsv $PROJ/revision_v2/audit/rps2_exons.GRCm39.110.bed $PROJ/revision_v2/audit/doi_check_final.txt $PROJ/revision_v2/audit/vw_tagged_fraction_by_group.txt revision_v2/audit/ 2>/dev/null || true
cp $PROJ/revision_v2/audit/f1_old/old_sites.bed $PROJ/revision_v2/audit/f1_old/sample_identity_verification.phaseA_F1rows.tsv revision_v2/audit/f1_old/ 2>/dev/null || true
cp $PROJ/revision_v2/submission/title_page.md $PROJ/revision_v2/submission/submission_checklist.md revision_v2/submission/ 2>/dev/null || true
cp $PROJ/analysis_v2/results/qc/*.tsv analysis_v2/results/qc/; cp $PROJ/analysis_v2/results/catalog/*.tsv analysis_v2/results/catalog/; cp $PROJ/analysis_v2/results/negative_control/*.tsv analysis_v2/results/negative_control/
cp $PROJ/analysis_v2/results/sensitivity_matched/*.tsv analysis_v2/results/sensitivity_matched/; cp $PROJ/analysis_v2/results/sensitivity_matched/bias/*.tsv analysis_v2/results/sensitivity_matched/bias/ 2>/dev/null || true
cp $PROJ/analysis_v2/results/predictors/*.tsv analysis_v2/results/predictors/; cp $PROJ/analysis_v2/results/within_fvb/*.tsv analysis_v2/results/within_fvb/
cp $PROJ/analysis_v2/results/f1_ase/v2/f1_ase_per_sample_v2.tsv $PROJ/analysis_v2/results/f1_ase/v2/f1_ase_summary_v2.tsv $PROJ/analysis_v2/results/f1_ase/v2/f1_ase_per_sample_v2.oldsites.tsv $PROJ/analysis_v2/results/f1_ase/v2/f1_ase_summary_v2.oldsites.tsv $PROJ/analysis_v2/results/f1_ase/v2/f1_wasp_allele_removal.tsv $PROJ/analysis_v2/results/f1_ase/v2/f1_sites_catalog_vs_other.tsv $PROJ/analysis_v2/results/f1_ase/v2/*.identity.tsv analysis_v2/results/f1_ase/v2/
cp $PROJ/analysis_v2/results/gene_meta.tsv $PROJ/analysis_v2/results/RESULTS.md5 analysis_v2/results/; cp $PROJ/analysis_v2/figures/*_data.tsv $PROJ/analysis_v2/figures/worked_examples.tsv $PROJ/analysis_v2/figures/FIGURES.md5 analysis_v2/figures/
cp -r $PROJ/manuscript/revised_v2/*.md $PROJ/manuscript/revised_v2/generated $PROJ/manuscript/revised_v2/supplementary_tables manuscript/revised_v2/
# Phase M (2026-10-23): identity per-site tables of Supplementary Fig. S1, Rps2 union recount, and the freeze/verification records
# (revision_v2/freeze/ without cache/, which holds large intermediate BEDs and is never copied)
mkdir -p analysis_v2/results/identity/full_site_tables analysis_v2/results/worked_examples revision_v2/freeze
cp $PROJ/analysis_v2/results/identity/full_site_tables/* analysis_v2/results/identity/full_site_tables/; cp $PROJ/analysis_v2/results/worked_examples/*.tsv analysis_v2/results/worked_examples/
find $PROJ/revision_v2/freeze -maxdepth 1 -type f -exec cp {} revision_v2/freeze/ \;
cp $PROJ/revision_v2/audit/fvb_vcf_counts.txt revision_v2/audit/; cp $PROJ/revision_v2/audit/wong2012_excerpts.txt revision_v2/audit/; cp $PROJ/analysis_v2/results/catalog/catalog_clusters.tsv analysis_v2/results/catalog/; cp $PROJ/analysis_v2/results/predictors/gene_exonic_snp_counts.tsv analysis_v2/results/predictors/
# Phase S2 (2026-10-23): GEO survey script (absolute paths replaced by arguments) and classification table, under survey/
mkdir -p survey; cp $PROJ/survey_release/stage3_survey_v2.py $PROJ/survey_release/survey_classification.tsv $PROJ/survey_release/README.md survey/
# Phase W-1 (2026-10-26): re-screen scripts, decisions, candidate list, per-series table and the cached GEO/SRA records
cp $PROJ/revision_v2/audit/sex_marker_counts.tsv revision_v2/audit/   # Phase W-3: sex markers of the F1 and PyVT samples
mkdir -p survey/raw; cp $PROJ/survey_release/rescreen_*.py $PROJ/survey_release/broad_search.py $PROJ/survey_release/search_log.tsv $PROJ/survey_release/candidates_broad.tsv $PROJ/survey_release/rescreen_decisions.tsv $PROJ/survey_release/rescreen_series.tsv $PROJ/survey_release/candidates_154.tsv $PROJ/survey_release/summary_stats_rescreen.txt survey/; cp $PROJ/survey_release/raw/* survey/raw/
# Phase S4 (2026-10-23): fixed 10% site subset of the identity check (identity_check_release/ -> identity_check/)
mkdir -p identity_check; cp $PROJ/identity_check_release/* identity_check/
git add analysis_v2 revision_v2 manuscript survey identity_check
# Phase S3 (2026-10-23): README.md and CITATION.cff generated by make_repo_readme.py; CITATION.cff is staged with mode 644 (the
# working tree carries an uncommitted v1-era mode change that is left as it is)
git add scripts/pipeline/01_build_references.sh   # Phase S3: LOG fix (edited in place; the file had no uncommitted changes)
cp $PROJ/repo_release/README.md README.md; cp $PROJ/repo_release/CITATION.cff CITATION.cff; git add README.md CITATION.cff; git update-index --chmod=-x CITATION.cff
# v1 scripts superseded by the revised pipeline -> legacy/v1/ on the branch. Index-only move: the committed v1 version is written to
# legacy/v1/, and the old path is removed from the index only (git rm --cached), so the working-tree files and the authors' uncommitted
# v1-era modifications in them stay on disk untouched. Idempotent. scripts/pipeline/01_build_references.sh stays (it built the
# references the revised analysis uses).
for p in scripts/figures/fig01_survey.py scripts/figures/fig02_f1_ase.py scripts/figures/fig03_strain_comp.py scripts/figures/fig04_predictors.py scripts/figures/fig05_catalog.py scripts/figures/fig06_worked_example.py scripts/figures/sfig_multimap.py scripts/figures/sfig_s1_multimap.py scripts/pipeline/02_align_wasp_featurecounts.sh scripts/pipeline/04_bias_analysis.py scripts/pipeline/05_multimapping_sensitivity.py tables/README_tables.md; do
  if git ls-files --error-unmatch "$p" > /dev/null 2>&1; then mkdir -p "legacy/v1/$(dirname $p)"; git show "HEAD:$p" > "legacy/v1/$p"; git add "legacy/v1/$p"; git rm --cached -q "$p"; fi
done
[ -f legacy/v1/README.md ] || { echo "Scripts of the originally submitted version (branch main), kept for reference; the revised analysis uses analysis_v2/ and revision_v2/scripts/ (see the top-level README)." > legacy/v1/README.md; git add legacy/v1/README.md; }
git commit -q -m "${COMMIT_MSG_HEAD:-v2-corrected: Phase H — F1 ASE re-derived at the identity-QC site set, bias decomposition, LAR format}

F1 allele-specific expression recomputed on re-aligned BAMs at the 203,510 FVB/NJ-homozygous exonic sites (the v1 values were
an artefact of a genotype-unfiltered site set; revision_v2/audit/f1_old_site_genotypes.md). Bias decomposed into FVB/NJ-read
gain and C57BL/6J-read loss (catalog unchanged, 138 genes). Manuscript assembled in LAR/BMC format with numbered references;
review-process mentions removed from the text. See revision_v2/CHANGELOG.md (Phase H) and revision_v2/text_changes.md.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01FXsgFVfisPBeFZ5TqiVozF"
git tag -d v2-corrected > /dev/null 2>&1 && echo "local tag v2-corrected (Phase F commit, never pushed) removed; the authors create the submission tag at the final commit" || true
echo "committed on v2-corrected: $(git log --oneline -1); NOT pushed"
