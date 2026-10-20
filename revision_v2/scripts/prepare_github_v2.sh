#!/bin/bash
# Rule 0.4-9: preserve history, no force-push. Creates branch v2-corrected in github_repo/ from main, copies the v2 code, manifest,
# QC results, figure scripts, supplementary tables and the revision record, inserts the README correction paragraph, commits and
# tags LOCALLY. It does NOT push: pushing and the tag are outward-facing and are done by the authors after the Phase F STOP
# (git -C github_repo push -u origin v2-corrected && git -C github_repo push origin v2-corrected-tag). Pre-existing uncommitted
# v1-era changes in github_repo (README.md, CITATION.cff, LICENSE, environment.yml, .gitignore) are left as they are and are NOT
# added, except README.md, which receives the correction paragraph (its other local edits come along; review the diff).
set -euo pipefail
PROJ=/usr/local/jupyter/FVB_bias; R=$PROJ/github_repo; cd $R
git rev-parse --verify v2-corrected >/dev/null 2>&1 && { echo "branch v2-corrected exists; checking it out"; git checkout v2-corrected; } || git checkout -b v2-corrected main
mkdir -p analysis_v2/pipeline revision_v2/scripts revision_v2/manifest revision_v2/audit analysis_v2/results/qc analysis_v2/results/catalog analysis_v2/results/negative_control analysis_v2/results/sensitivity_matched analysis_v2/results/predictors analysis_v2/results/within_fvb analysis_v2/figures manuscript/revised_v2
cp $PROJ/analysis_v2/pipeline/*.py $PROJ/analysis_v2/pipeline/*.sh $PROJ/analysis_v2/pipeline/*.R analysis_v2/pipeline/; cp $PROJ/analysis_v2/params.yaml analysis_v2/
cp $PROJ/revision_v2/scripts/*.py $PROJ/revision_v2/scripts/*.sh revision_v2/scripts/; cp $PROJ/revision/scripts/check_figure_order.py revision_v2/scripts/
cp $PROJ/revision_v2/manifest/manifest.tsv $PROJ/revision_v2/manifest/design_proposal.md revision_v2/manifest/
cp $PROJ/revision_v2/CHANGELOG.md $PROJ/revision_v2/text_changes.md $PROJ/revision_v2/summary_of_changes.md $PROJ/revision_v2/final_check.md revision_v2/
cp $PROJ/revision_v2/audit/sample_audit.tsv $PROJ/revision_v2/audit/audit_summary.md $PROJ/revision_v2/audit/gse55792_qc_note.md $PROJ/revision_v2/audit/qc_calibration.md $PROJ/revision_v2/audit/junction_sets_fvb_index_vs_b6_gtf.txt $PROJ/revision_v2/audit/fvb_vcf_counts.txt revision_v2/audit/ 2>/dev/null || true
cp $PROJ/analysis_v2/results/qc/*.tsv analysis_v2/results/qc/; cp $PROJ/analysis_v2/results/catalog/*.tsv analysis_v2/results/catalog/; cp $PROJ/analysis_v2/results/negative_control/*.tsv analysis_v2/results/negative_control/
cp $PROJ/analysis_v2/results/sensitivity_matched/*.tsv analysis_v2/results/sensitivity_matched/; cp $PROJ/analysis_v2/results/predictors/*.tsv analysis_v2/results/predictors/; cp $PROJ/analysis_v2/results/within_fvb/*.tsv analysis_v2/results/within_fvb/
cp $PROJ/analysis_v2/results/gene_meta.tsv $PROJ/analysis_v2/results/RESULTS.md5 analysis_v2/results/; cp $PROJ/analysis_v2/figures/*_data.tsv $PROJ/analysis_v2/figures/worked_examples.tsv $PROJ/analysis_v2/figures/FIGURES.md5 analysis_v2/figures/
cp -r $PROJ/manuscript/revised_v2/*.md $PROJ/manuscript/revised_v2/generated $PROJ/manuscript/revised_v2/supplementary_tables manuscript/revised_v2/
grep -q "## Correction (revision 2" README.md || { awk -v para="$(cat $PROJ/revision_v2/README_correction_paragraph.md)" 'NR==1{print; print ""; print para; next} {print}' README.md > README.tmp && mv README.tmp README.md; }
git add analysis_v2 revision_v2 manuscript README.md
git commit -m "v2-corrected: re-derived pure-strain analysis from manifest-selected, sequence-verified samples

Dataset selection and sample identification of the original submission were found to be in error during revision
(see revision_v2/CHANGELOG.md and README correction paragraph). This commit adds the corrected pipeline (analysis_v2/),
the sample manifest, identity-QC results, figure scripts, generated tables and the revised manuscript. The original code
remains unchanged in the history of main.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01FXsgFVfisPBeFZ5TqiVozF"
git tag -a v2-corrected -m "Revision 2: corrected dataset identification; pure-strain analysis re-derived" 2>/dev/null || echo "tag v2-corrected already exists"
echo "branch v2-corrected committed and tagged locally; NOT pushed. Review: git -C $R log --stat -1 | head; push with: git -C $R push -u origin v2-corrected && git -C $R push origin v2-corrected"
