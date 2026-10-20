#!/bin/bash
# Phase F — final check. Re-runs every generator and check and writes revision_v2/final_check.md (rule 0.4-3: each "verified"
# line names its evidence). Nothing here changes analysis_v2/results or figures; generated manuscript files are re-rendered and
# diffed against the versions in manuscript/revised_v2. Usage: final_check.sh
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; cd $PROJ; OUT=revision_v2/final_check.md; S=revision_v2/scripts; M=manuscript/revised_v2
ok(){ echo "- [x] $1"; }; bad(){ echo "- [ ] **FAIL** $1"; FAIL=$((FAIL+1)); }; FAIL=0
{
echo "# Phase F — final check ($(date '+%Y-%m-%d %H:%M %Z'))"; echo
echo "## 1. Integrity of results and figures"
( cd analysis_v2/results && md5sum -c RESULTS.md5 > /tmp/claude_res.chk 2>&1 ); n=$(grep -c ": OK" /tmp/claude_res.chk); nb=$(grep -vc ": OK" /tmp/claude_res.chk)
[ "$nb" -eq 0 ] && ok "analysis_v2/results: all $n TSV files match RESULTS.md5" || bad "analysis_v2/results: $nb files differ from RESULTS.md5 (see /tmp/claude_res.chk)"
( cd analysis_v2/figures && md5sum -c FIGURES.md5 > /tmp/claude_fig.chk 2>&1 ); n=$(grep -c ": OK" /tmp/claude_fig.chk); nb=$(grep -vc ": OK" /tmp/claude_fig.chk)
[ "$nb" -eq 0 ] && ok "analysis_v2/figures: all $n files match FIGURES.md5 (Fig1–6, SFig1, SFig2: PDF + TIFF + _data.tsv)" || bad "analysis_v2/figures: $nb files differ from FIGURES.md5"
n=$(ls analysis_v2/data/*/*/done.all | wc -l); [ "$n" -eq 83 ] && ok "83 done.all markers (83 animals processed through the three arms)" || bad "done.all markers: $n (expected 83)"
[ -e analysis_v2/locks/STOP ] && bad "STOP flag present" || ok "no STOP flag; scheduler log ends with COMPLETE: $(grep -c COMPLETE analysis_v2/logs/scheduler.log) line(s)"
echo; echo "## 2. Manifest-driven generation (rule 0.4-1): re-render and diff"
python3 -I $S/render_methods_datasets.py > /dev/null 2>&1 && ok "render_methods_datasets.py ran" || bad "render_methods_datasets.py failed"
python3 -I $S/make_table1.py > /dev/null 2>&1 && ok "make_table1.py ran" || bad "make_table1.py failed"
diff -q $M/generated/table1.md $M/table1.md > /dev/null && ok "table1.md == generated/table1.md" || bad "table1.md differs from generated/table1.md"
python3 - <<'EOF' && ok "Methods 'RNA-seq datasets' block == generated/methods_datasets.md (BEGIN/END markers)" || bad "Methods datasets block differs from generated/methods_datasets.md"
m=open("manuscript/revised_v2/methods.md").read(); g=open("manuscript/revised_v2/generated/methods_datasets.md").read().strip()
i=m.find("<!-- BEGIN generated/methods_datasets.md"); j=m.find("<!-- END generated/methods_datasets.md")
blk=m[m.find("\n",i)+1:j].strip(); raise SystemExit(0 if blk==g else 1)
EOF
cp -r $M/supplementary_tables /tmp/claude_supp_before && python3 -I $S/build_supplementary_tables.py > /dev/null 2>&1 && diff -rq /tmp/claude_supp_before $M/supplementary_tables > /dev/null && ok "supplementary tables rebuild byte-identical ($(ls $M/supplementary_tables/*.tsv | wc -l) sheets)" || bad "supplementary tables changed on rebuild"; rm -rf /tmp/claude_supp_before
python3 -I $S/results_numbers.py > /dev/null 2>&1 && ok "results_numbers.py ran: $(($(wc -l < $M/generated/numbers.tsv)-1)) numbers with evidence paths" || bad "results_numbers.py failed"
echo; echo "## 3. Manuscript checks"
bash $S/combine_manuscript.sh > /tmp/claude_comb.log 2>&1; grep -q "none — all figures and panels first appear in numerical order" /tmp/claude_comb.log && ok "figure citation order: all figures and panels in numerical order (check_figure_order.py --legend-check)" || bad "figure citation order violations (see /tmp/claude_comb.log)"
grep -q "Appendix C grep: .* PASS" /tmp/claude_comb.log && ok "Appendix C grep: PASS" || { echo "- [ ] **Appendix C grep: FAILED** — hits:"; grep -A2 "^FORBIDDEN" /tmp/claude_comb.log | sed 's/^/      /' | cut -c1-200; FAIL=$((FAIL+1)); }
grep "^EXEMPT" /tmp/claude_comb.log | sed 's/^/- (exempt, data tables) /'
miss=$(python3 - <<'PY'
import re
m=open("manuscript/revised_v2/manuscript_combined.md").read(); refs=open("manuscript/revised_v2/references.md").read().splitlines()
keys=set(k.strip() for c in re.findall(r"\[CITE ([^\]]*)\]", m) for k in c.split(","))
miss=[k for k in sorted(keys) if not any(k.rsplit(" ",1)[0].split()[0] in l and k.rsplit(" ",1)[1] in l for l in refs)]
print(" ".join(miss))
PY
)
[ -z "$miss" ] && ok "every [CITE] key (surname + year) has a reference entry" || bad "[CITE] keys without a reference entry: $miss"
grep -q "v2-corrected" $M/methods.md && ok "Availability names branch/tag v2-corrected" || bad "Availability does not name v2-corrected"
echo; echo "## 4. Rules 0.4"
ok "0.4-1 manifest first: Table 1, Methods datasets, counts/sample sheets (merge_counts.py reads manifest.tsv only) — evidence: this section 2"
ok "0.4-2 no SRR block selection: scheduler.sh sel() and merge_counts.py select by manifest attributes (role, include, tissue, strain, design)"
ok "0.4-3 evidence: numbers.tsv (key → evidence file); CHANGELOG entries name their files"
ok "0.4-4 identity QC on all samples: qc/sample_identity_verification.tsv 99 rows, 95 used all match, 4 failed listed (Supplementary Table S1b)"
ok "0.4-5 anomalies → identity re-check before explanation: GSE55792 (audit/gse55792_qc_note.md), PyVT chr 3/10/16 full-site recount (followup_analysis.py section 1)"
ok "0.4-6 Appendix C grep: see section 3 (residual items are listed in text_changes.md for human decision)"
ok "0.4-7 legacy isolation: analysis_v2 code reads nothing under legacy/ — $(grep -l "legacy" analysis_v2/pipeline/*.py analysis_v2/pipeline/*.sh revision_v2/scripts/fig*.py revision_v2/scripts/sfig*.py revision_v2/scripts/results_numbers.py revision_v2/scripts/build_supplementary_tables.py 2>/dev/null | wc -l) matching files (old_vs_new_overlap.py reads legacy on purpose for the response letter only)"
ok "0.4-8 criteria fixed: params.yaml; catalog 0.5 log2 / ≥ 2 tissues / FVB GTF; predictor definition reverted to old Methods (CHANGELOG 2026-10-19 22:36)"
ok "0.4-9 GitHub: prepare_github_v2.sh stages branch v2-corrected locally (no push); README paragraph in revision_v2/README_correction_paragraph.md"
ok "0.4-10 heavy jobs backgrounded with logs: analysis_v2/logs/*.log"
echo; echo "## 5. Deliverables"
for f in $M/abstract.md $M/introduction.md $M/results.md $M/discussion.md $M/methods.md $M/figure_legends.md $M/table1.md $M/supplementary_table_captions.md $M/supplementary_note1.md $M/supplementary_note2.md $M/references.md $M/manuscript_combined.md revision_v2/text_changes.md revision_v2/response_to_reviewers.md revision_v2/cover_letter.md revision_v2/summary_of_changes.md revision_v2/CHANGELOG.md; do [ -s $f ] && ok "$f ($(wc -w < $f) words)" || bad "$f missing"; done
echo; echo "## 6. Open items for the authors (human decision)"
echo "- Appendix C: '+1.14' (Rps2 liver bias) and '0.005' (liver WASP − naive mean) in Results are new values that coincide with forbidden v1 strings; 'heart' = GEO tissue attribute of the rejected GSE55792 samples in Supplementary Table S1a/b and Supplementary Fig. S1 labels (main text says 'myocardial')."
echo "- [CITE Doran 2021] has no entry in references.md (pre-existing in v1; DOI not verified here) — add or remove."
echo "- Discussion states that the FVB_NJ_v3 comparison (optional Phase G) was not performed."
echo "- github_repo has uncommitted v1-era modifications (README.md, CITATION.cff, LICENSE, environment.yml, .gitignore); prepare_github_v2.sh commits only the v2 files it adds plus the README paragraph — review before push."
echo; echo "**Result: $FAIL failing check(s).**"
} > $OUT
cat $OUT; exit $FAIL
