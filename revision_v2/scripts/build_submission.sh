#!/bin/bash
# Phase H-8: build the deliverables in revision_v2/submission/.
#   build     : LAR docx (assemble_lar.py) + blue-marked copy LAR_main_text_v2_marked.docx (mark_changes_docx.py, H-9), supplementary_tables_v2.xlsx (S1–S5 sheets from manuscript/revised_v2/supplementary_tables/*.tsv),
#               figures/ (Fig 1–6 + SFig S1–S3, PDF + 600-dpi TIFF copied from analysis_v2/figures), response_to_reviewers_v2.docx,
#               cover_letter_v2.docx (md_to_docx.py)
#   checklist : submission_checklist.md from revision_v2/final_check.md + the github_repo state (run after final_check.sh and commit_phase_h.sh)
# Usage: build_submission.sh build|checklist
set -euo pipefail
PROJ=/usr/local/jupyter/FVB_bias; cd $PROJ; S=revision_v2/scripts; M=manuscript/revised_v2; SUB=revision_v2/submission; mkdir -p $SUB/figures
case "${1:-build}" in
build)
  python3 $S/assemble_lar.py
  python3 $S/make_change_ledger.py | tail -3 && python3 $S/mark_from_ledger.py $SUB/LAR_main_text_v2.docx $SUB/LAR_main_text_v2_marked.docx   # Phase J-6: blue-marked copy from the change ledger (clean docx untouched)
  python3 - <<'EOF'
import pandas as pd, re
P="/usr/local/jupyter/FVB_bias"; T=f"{P}/manuscript/revised_v2/supplementary_tables"
sheets=[("S1a","TableS1a_manifest"),("S1b","TableS1b_identity_verification"),("S2","TableS2_catalog"),("S3a","TableS3a_sensitivity_primary"),("S3b","TableS3b_sensitivity_control_only"),
        ("S3c","TableS3c_sensitivity_matched_reference"),("S3d","TableS3d_negative_control_per_animal"),("S3e","TableS3e_negative_control_heldout"),("S3f","TableS3f_pyvt_residual_segments"),
        ("S4","TableS4_alignment_stats_per_animal_arm"),("S4b","TableS4b_wasp_attrition_by_group"),("S5","TableS5_survey_of_practice")]
cap=open(f"{P}/manuscript/revised_v2/supplementary_table_captions.md",encoding="utf-8").read()
titles={m.group(1):m.group(2) for m in re.finditer(r"\*\*Supplementary Table (S\d)\. ([^*]+)\*\*", cap)}
with pd.ExcelWriter(f"{P}/revision_v2/submission/supplementary_tables_v2.xlsx", engine="openpyxl") as w:
    pd.DataFrame([{"sheet":s,"file":f+".tsv","table":titles.get(re.match(r"S\d",s).group(0),"")} for s,f in sheets]).to_excel(w, sheet_name="Contents", index=False)
    for s,f in sheets:
        df=pd.read_csv(f"{T}/{f}.tsv", sep="\t", dtype=str, keep_default_na=False); df.to_excel(w, sheet_name=s, index=False); print(s, f, df.shape)
EOF
  for f in Fig1_survey Fig2_f1_ase Fig3_strain_comparison Fig4_predictors Fig5_catalog Fig6_worked_examples SFig1_GSE55792_haplotype_blocks SFig2_multimapping SFig3_bias_decomposition; do cp analysis_v2/figures/$f.pdf analysis_v2/figures/$f.tiff $SUB/figures/; done
  python3 $S/md_to_docx.py revision_v2/response_to_reviewers.md $SUB/response_to_reviewers_v2.docx "Response to Reviewers"
  python3 $S/md_to_docx.py revision_v2/cover_letter.md $SUB/cover_letter_v2.docx "Cover letter"
  for f in response_to_reviewers_v2 cover_letter_v2 LAR_main_text_v2_marked; do pandoc $SUB/$f.docx -t plain --wrap=none -o $SUB/$f.plain.txt; done
  ls -la $SUB $SUB/figures ;;
checklist)
  R=github_repo
  { echo "# Submission checklist (Phase H-8, $(date '+%Y-%m-%d %H:%M %Z'))"; echo
    echo "## Deliverables in revision_v2/submission/"; for f in LAR_main_text_v2.docx LAR_additional_file_1.docx supplementary_tables_v2.xlsx response_to_reviewers_v2.docx cover_letter_v2.docx; do echo "- $f ($(stat -c %s $SUB/$f) bytes)"; done
    echo "- figures/: $(ls $SUB/figures/*.pdf | wc -l) PDF + $(ls $SUB/figures/*.tiff | wc -l) TIFF ($(ls $SUB/figures | sed 's/\..*//' | sort -u | tr '\n' ' '))"; echo
    echo "## Checks (from revision_v2/final_check.md)"; sed -n '/^## 3\./,/^## 4\./p' revision_v2/final_check.md | grep '^- \[' ; sed -n '/^## 1\./,/^## 2\./p' revision_v2/final_check.md | grep '^- \['; grep '^\*\*Result' revision_v2/final_check.md; echo
    echo "## GitHub state (rule 0.4-9)"; echo "- repository: $R, branch: $(git -C $R rev-parse --abbrev-ref HEAD)"; echo "- commits on v2-corrected not on main: $(git -C $R rev-list --count main..v2-corrected)"; git -C $R log --format='- %h %ad %s' --date=short main..v2-corrected | head -5
    echo "- tags: $(git -C $R tag | tr '\n' ' ' | sed 's/ *$//'); none pushed (NOT pushed: the authors push the branch and create the submission tag — \`git -C github_repo push -u origin v2-corrected\`, then tag)"
    echo "- uncommitted v1-era modifications left in place (review before push): $(git -C $R status --short | wc -l) paths"; echo
    echo "## Items for the authors before submission"; echo "- Read LAR_main_text_v2.docx and LAR_additional_file_1.docx in full (plain-text dumps: LAR_main_text_v2.plain.txt, LAR_additional_file_1.plain.txt)."
    echo "- Confirm the structured Abstract (≤ 350 words), the Rps2 positioning (limit of a SNP-only personalised reference) and the F1 wording (2.0-point deficit; WASP overcorrects)."
    echo "- Push branch v2-corrected and create the tag; then confirm the repository URL in the Availability statement resolves to the tagged version."
  } > $SUB/submission_checklist.md; cat $SUB/submission_checklist.md ;;
esac
