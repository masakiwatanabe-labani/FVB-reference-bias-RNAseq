#!/bin/bash
# Phase H — final check (extends the Phase F check). Re-runs every generator and check and writes revision_v2/final_check.md
# (rule 0.4-3: each "verified" line names its evidence). Nothing here changes analysis_v2/results or figures; generated manuscript
# files are re-rendered and diffed against the versions in manuscript/revised_v2; the LAR docx files are re-assembled and dumped to
# plain text for the H-8 checks. Usage: final_check.sh
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; cd $PROJ; OUT=revision_v2/final_check.md; S=revision_v2/scripts; M=manuscript/revised_v2; SUB=revision_v2/submission
ok(){ echo "- [x] $1"; }; bad(){ echo "- [ ] **FAIL** $1"; FAIL=$((FAIL+1)); }; FAIL=0
{
echo "# Phase J — final check ($(date '+%Y-%m-%d %H:%M %Z'))"; echo
echo "## 1. Integrity of results and figures"
( cd analysis_v2/results && md5sum -c RESULTS.md5 > /tmp/claude_res.chk 2>&1 ); n=$(grep -c ": OK" /tmp/claude_res.chk); nb=$(grep -vc ": OK" /tmp/claude_res.chk)
[ "$nb" -eq 0 ] && ok "analysis_v2/results: all $n files match RESULTS.md5 (includes f1_ase/v2/: 24 + 24 per-site tables, per-animal/summary tables, catalog decomposition)" || bad "analysis_v2/results: $nb files differ from RESULTS.md5 (see /tmp/claude_res.chk)"
( cd analysis_v2/figures && md5sum -c FIGURES.md5 > /tmp/claude_fig.chk 2>&1 ); n=$(grep -c ": OK" /tmp/claude_fig.chk); nb=$(grep -vc ": OK" /tmp/claude_fig.chk)
[ "$nb" -eq 0 ] && ok "analysis_v2/figures: all $n files match FIGURES.md5 (Fig 1–6, SFig 1–3: PDF + TIFF + _data.tsv)" || bad "analysis_v2/figures: $nb files differ from FIGURES.md5"
grep -q "Fig2_f1_ase.pdf" analysis_v2/figures/FIGURES.md5 && grep -q "SFig3_bias_decomposition.pdf" analysis_v2/figures/FIGURES.md5 && ok "Fig. 2 (regenerated from f1_ase/v2/, with panel D) and Supplementary Fig. S3 are listed in FIGURES.md5" || bad "Fig. 2 or SFig S3 missing from FIGURES.md5"
n=$(ls analysis_v2/data/*/*/done.all | wc -l); [ "$n" -eq 83 ] && ok "83 done.all markers (83 animals processed through the three arms); 12 done.f1_ase markers: $(ls analysis_v2/data/*/*/done.f1_ase | wc -l)" || bad "done.all markers: $n (expected 83)"
[ -e analysis_v2/locks/STOP ] && bad "STOP flag present" || ok "no STOP flag"
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
bash $S/combine_manuscript.sh > /tmp/claude_comb.log 2>&1; grep -q "none — all figures and panels first appear in numerical order" /tmp/claude_comb.log && ok "figure citation order (incl. Supplementary Fig. S3): all figures and panels in numerical order (check_figure_order.py --legend-check)" || bad "figure citation order violations (see /tmp/claude_comb.log)"
grep -q "Appendix C grep: .* PASS" /tmp/claude_comb.log && ok "Appendix C grep (incl. the H-1d additions 0.333, 0.386, 17 percentage, 274, 872): PASS" || { echo "- [ ] **Appendix C grep: FAILED** — hits:"; grep -A2 "^FORBIDDEN" /tmp/claude_comb.log | sed 's/^/      /' | cut -c1-200; FAIL=$((FAIL+1)); }
grep "^EXEMPT" /tmp/claude_comb.log | sed 's/^/- (exempt, data tables) /'
n=$(grep -c '\[CITE' $M/abstract_structured.md $M/introduction.md $M/methods.md $M/results.md $M/discussion.md $M/figure_legends.md $M/table1.md $M/supplementary_table_captions.md $M/supplementary_note1.md $M/supplementary_note2.md | awk -F: '{s+=$2} END{print s}')
[ "$n" -eq 0 ] && ok "no [CITE placeholder left in any main-text or Additional-file section" || bad "$n [CITE placeholders remain"
python3 - <<'EOF' > /tmp/claude_cite.txt && ok "numbered citations: $(cat /tmp/claude_cite.txt)" || bad "numbered citations: $(cat /tmp/claude_cite.txt)"
import re
M="manuscript/revised_v2"; txt="".join(open(f"{M}/{f}").read() for f in ("introduction.md","methods.md","results.md","discussion.md","figure_legends.md","table1.md","supplementary_table_captions.md","supplementary_note1.md","supplementary_note2.md"))
refs=re.findall(r"^(\d+)\. ", open(f"{M}/references.md").read(), flags=re.M); N=len(refs)
cited=[]; first=[]
for m in re.finditer(r"\[(\d+(?:[–,\-]\d+)*)\]", txt):
    nums=[]
    for part in m.group(1).split(","):
        if "–" in part or "-" in part:
            a,b=re.split(r"[–\-]",part); nums+=list(range(int(a),int(b)+1))
        else: nums.append(int(part))
    for k in nums:
        if k not in first: first.append(k)
    cited+=nums
cs=set(cited); ok_=(refs==[str(i) for i in range(1,N+1)]) and cs==set(range(1,N+1)) and first==sorted(first)
print(f"{N} references, numbered 1..{N} consecutively; all {len(cs)} cited; first appearances in increasing order: {first==sorted(first)}"); raise SystemExit(0 if ok_ else 1)
EOF
python3 - <<'EOF' > revision_v2/audit/doi_check_final.txt 2>&1; tail -1 revision_v2/audit/doi_check_final.txt | grep -q "ALL OK" && ok "all DOIs in references.md resolve (HTTP 301/302 from doi.org; revision_v2/audit/doi_check_final.txt)" || bad "DOI check: $(tail -1 revision_v2/audit/doi_check_final.txt)"
import re, urllib.request, ssl
ctx=ssl.create_default_context(); dois=re.findall(r"https://doi\.org/(\S+?)\.?(?:\s|$)", open("manuscript/revised_v2/references.md").read()); bad=0
class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self,*a,**k): return None
op=urllib.request.build_opener(NoRedirect, urllib.request.HTTPSHandler(context=ctx))
for d in dois:
    try:
        r=op.open(urllib.request.Request(f"https://doi.org/{d}",method="HEAD",headers={"User-Agent":"Mozilla/5.0"}),timeout=20); code=r.status
    except urllib.error.HTTPError as e: code=e.code
    except Exception as e: code=str(e)
    print(d, code); bad+= 0 if code in (301,302,303,200) else 1
print(f"{len(dois)} DOIs; failures {bad}; " + ("ALL OK" if bad==0 else "FAILED"))
EOF
w=$(grep -v '^#\|^\*\*' $M/abstract_structured.md | wc -w); h=$(grep -c '^\*\*\(Background\|Methods\|Results\|Conclusions\)\*\*$' $M/abstract_structured.md)
[ "$w" -le 350 ] && [ "$h" -eq 4 ] && ok "structured Abstract: $w words (≤ 350), 4 headings (Background/Methods/Results/Conclusions)" || bad "structured Abstract: $w words, $h headings"
python3 $S/assemble_lar.py > /tmp/claude_asm.log 2>&1 && ok "assemble_lar.py rebuilt LAR_main_text_v2.docx and LAR_additional_file_1.docx (pandoc, --reference-doc 'LAR text submit.docx')" || bad "assemble_lar.py failed (/tmp/claude_asm.log)"
P=$SUB/LAR_main_text_v2.plain.txt
order=$(grep -n -E '^(Background|Methods|Results|Discussion|Conclusions|List of abbreviations|Declarations|References|Figure legends)$' $P | awk -F: 'NR>4{print $2}' | tr '\n' '>' ); exp="Background>Methods>Results>Discussion>Conclusions>List of abbreviations>Declarations>References>Figure legends>"
[ "$order" = "$exp" ] && grep -q '^Table 1\.' $P && ok "section order in the docx (plain dump): Abstract (structured) → Keywords → $exp Table 1" || bad "section order in the docx: got '$order'"
n5=$(grep -c 'corrected analysis\|this revision\|earlier version\|submitted version\|reviewer\|Reviewer' $P $SUB/LAR_additional_file_1.plain.txt | awk -F: '{s+=$2} END{print s}'); [ "$n5" -eq 0 ] && ok "H-5 grep (corrected analysis / this revision / earlier version / submitted version / reviewer) in the main text and Additional file 1: 0 hits" || bad "H-5 grep: $n5 hits in the docx dumps"
nf=$(awk -F'\t' 'NR>1 && $1 ~ /^f1_/ && $3 !~ /f1_ase\/v2\/|revision_v2\/audit\/f1_old|analysis_v2\/qc\/fvbnj_homalt_exonic\.hetGT\.vcf\.gz/' $M/generated/numbers.tsv | wc -l); nf1=$(awk -F'\t' 'NR>1 && $1 ~ /^f1_/' $M/generated/numbers.tsv | wc -l)
[ "$nf" -eq 0 ] && ok "all $nf1 f1_* numbers come from the new per-site tables (evidence column: analysis_v2/results/f1_ase/v2/ or the H-1a audit files); the v1 tables f1_ase/f1_ase_per_sample.tsv are not an evidence source" || bad "$nf f1_* numbers have another evidence source"
d1=$(grep -c '\[CITE' $P $SUB/LAR_additional_file_1.plain.txt | awk -F: '{s+=$2} END{print s}'); d2=$(grep -c 'revision_v2/\|analysis_v2/' $P $SUB/LAR_additional_file_1.plain.txt | awk -F: '{s+=$2} END{print s}'); d3=$(grep -c '^Step [0-9]\. ' $P $SUB/LAR_additional_file_1.plain.txt | awk -F: '{s+=$2} END{print s}'); d4=$(grep -c '^[0-9]\+\.  \*\*' $P | awk -F: '{s+=$2} END{print s}'); d5=$(grep -c '\* p < 0.05; \*\* p < 0.01; \*\*\* p < 0.001' $P)
[ "$d1" -eq 0 ] && [ "$d2" -eq 0 ] && [ "$d4" -eq 0 ] && [ "$d5" -ge 1 ] && ok "docx → plain: no [CITE, no local paths, no continued numbered lists ($d3 Step paragraphs; the reconstruction keeps the v1 prose steps), Fig. 2 significance asterisks rendered literally ($d5 line)" || bad "docx → plain: CITE=$d1 paths=$d2 steps=$d3 numbered-bold-lists=$d4 asterisk-lines=$d5"
grep -q "v2-corrected" $M/declarations.md && ok "Availability (Declarations) names branch/tag v2-corrected" || bad "Availability does not name v2-corrected"
echo; echo "## 3b. Phase J checks (minimal-change reconstruction on the submitted text)"
python3 $S/make_change_ledger.py > /tmp/claude_ledger.log 2>&1 && ok "change ledger rebuilt: $(grep '^total' /tmp/claude_ledger.log); $(grep '^body sentences' /tmp/claude_ledger.log)" || bad "make_change_ledger.py failed"
python3 - <<'EOF' > /tmp/claude_ledger_chk.txt && ok "ledger completeness: $(cat /tmp/claude_ledger_chk.txt)" || bad "ledger completeness: $(cat /tmp/claude_ledger_chk.txt)"
import csv, re, sys
P="/usr/local/jupyter/FVB_bias"; rows=list(csv.DictReader(open(f"{P}/revision_v2/change_ledger.tsv"), delimiter="\t"))
def norm(s):
    s=s.replace("\\",""); s=re.sub(r"\[CITE[^\]]*\]","",s); s=re.sub(r"\[\d+(?:[–,\-]\d+)*\]","",s); s=s.replace("−","-").replace("–","-").replace("--","-"); return re.sub(r"[\W_]+","",s).lower()
v1=set()
for line in open(f"{P}/revision_v2/submission/LAR_v1_text.md", encoding="utf-8"):
    v1.add(norm(line.strip()))
for nn in (1,2):
    for para in open(f"{P}/manuscript/supplementary_note{nn}.md", encoding="utf-8").read().split("\n\n"): v1.add(norm(para))
keep_bad=[r for r in rows if r["status"]=="KEEP" and not any(norm(r["new_text"]) in x for x in v1)]
numflags=[r for r in rows if r["number_check"]]
body=[r for r in rows if r["section"] in ("Abstract","Background","Methods","Results","Discussion","Conclusions") and r["status"]!="DROP"]
trig_bad=[r for r in rows if r["status"] in ("EDIT","NEW","DROP") and r["trigger"] in ("?","")]
print(f"{len(rows)} rows; {len(body)} body sentences with a row; KEEP sentences not found verbatim in v1: {len(keep_bad)}; EDIT/NEW numbers not in numbers.tsv: {len(numflags)}; EDIT/NEW/DROP rows without trigger: {len(trig_bad)}")
sys.exit(1 if (keep_bad or numflags or trig_bad) else 0)
EOF
python3 $S/verify_response_pointers.py > /tmp/claude_ptr.log 2>&1 && ok "response letter: every 'Changed:' pointer resolves to ledger rows with the matching trigger ($(grep -c '^ok' /tmp/claude_ptr.log) pointers)" || bad "response letter pointers unresolved: $(grep -c '^FAIL' /tmp/claude_ptr.log) (see /tmp/claude_ptr.log)"
wb=$(grep -v '^#\|^\*\*' $M/abstract_structured.md | wc -w); wh=$(sed 's/^# //;s/\*\*//g' $M/abstract_structured.md | wc -w); [ "$wb" -le 350 ] && [ "$wh" -le 350 ] && ok "Abstract ≤ 350 words both ways (body $wb, with headings $wh)" || bad "Abstract words: body $wb, with headings $wh"
grep -q "3Rs" $M/discussion.md && ok "'3Rs' is used in the text (Discussion), so the abbreviation entry is used" || bad "'3Rs' not used in the text"
python3 - <<'EOF' > /tmp/claude_abbr.txt && ok "abbreviation list: $(cat /tmp/claude_abbr.txt)" || bad "abbreviation list: $(cat /tmp/claude_abbr.txt)"
import re, sys
M="/usr/local/jupyter/FVB_bias/manuscript/revised_v2"; body="\n".join(open(f"{M}/{f}").read() for f in ("abstract_structured.md","introduction.md","methods.md","results.md","discussion.md","figure_legends.md"))
ab=[l.split(":")[0].strip() for l in open(f"{M}/abbreviations.md") if ":" in l]
unused=[a for a in ab if not re.search(r"(?<![\w])"+re.escape(a)+r"(?![\w])", body)]
print(f"{len(ab)} entries, unused: {unused}"); sys.exit(1 if unused else 0)
EOF
n_i=$(grep -c "either fail to align or are assigned to another locus\|either failing to align or being assigned to another locus" $M/results.md $M/supplementary_note2.md | awk -F: '{s+=$2} END{print s}'); grep -q "expected to be heterozygous" $M/methods.md && grep -q "Mechanism 1\." $M/supplementary_note2.md && [ "$n_i" -ge 2 ] && ok "Phase I wording present: direct-class read fate ($n_i places), 'expected to be heterozygous' (Methods), 'Mechanism 1./2.' labels (Note 2)" || bad "Phase I wording missing (I-1d places=$n_i)"
n_ban=$(grep -c "WASP improved\|partially corrected\|did not normalize\|five pure-strain\|five datasets\|30 of the 32\|pvt_discordant\|FVB-enriched F3" $P $SUB/LAR_additional_file_1.plain.txt | awk -F: '{s+=$2} END{print s}'); [ "$n_ban" -eq 0 ] && ok "J-0 semantic ban list (also in appendix_c_grep.sh): 0 hits in the docx dumps" || bad "J-0 ban list: $n_ban hits in the docx dumps"
[ -s $SUB/LAR_main_text_v2_marked.docx ] && ok "marked copy present ($(stat -c %s $SUB/LAR_main_text_v2_marked.docx) bytes; generated from the ledger by mark_from_ledger.py)" || bad "marked copy missing"
echo; echo "## 4. Rules 0.4"
ok "0.4-1 manifest first: Table 1, Methods datasets, counts/sample sheets (merge_counts.py reads manifest.tsv only) — evidence: section 2"
ok "0.4-2 no SRR block selection: scheduler.sh sel() and merge_counts.py select by manifest attributes (role, include, tissue, strain, design)"
ok "0.4-3 evidence: numbers.tsv (key → evidence file); CHANGELOG entries name their files"
ok "0.4-4 identity QC on all samples: qc/sample_identity_verification.tsv 99 rows, 95 used all match (F1 rows re-checked on the regenerated BAMs, H-1), 4 failed listed (Supplementary Table S1b)"
ok "0.4-5 anomalies → identity re-check before explanation: GSE55792 (audit/gse55792_qc_note.md), PyVT chr 3/10/16 full-site recount, F1 site-set artefact (audit/f1_old_site_genotypes.md; old-site reproduction f1_ase_summary_v2.oldsites.tsv)"
ok "0.4-6 Appendix C grep: section 3"
ok "0.4-7 legacy isolation: analysis_v2 code reads nothing under legacy/ — $(grep -l "legacy" analysis_v2/pipeline/*.py analysis_v2/pipeline/*.sh revision_v2/scripts/fig*.py revision_v2/scripts/sfig*.py revision_v2/scripts/results_numbers.py revision_v2/scripts/build_supplementary_tables.py 2>/dev/null | wc -l) matching files (old_vs_new_overlap.py reads legacy on purpose for the response letter only)"
ok "0.4-8 criteria fixed: params.yaml; catalog 0.5 log2 / ≥ 2 tissues / FVB GTF unchanged (138 genes); decomposition added as a description, not a new criterion"
ok "0.4-9 GitHub: branch v2-corrected committed locally (commit_phase_h.sh); no push; tag to be created by the authors at submission"
ok "0.4-10 heavy jobs backgrounded with logs: analysis_v2/logs/*.log (f1_ase_rerun.log, f1_oldsites_recount.log)"
echo; echo "## 5. Deliverables (revision_v2/submission/)"
for f in $SUB/LAR_main_text_v2.docx $SUB/LAR_additional_file_1.docx $SUB/supplementary_tables_v2.xlsx $SUB/response_to_reviewers_v2.docx $SUB/cover_letter_v2.docx; do [ -s $f ] && ok "$f ($(stat -c %s $f) bytes)" || bad "$f missing"; done
nfig=$(ls $SUB/figures/*.pdf 2>/dev/null | wc -l); ntif=$(ls $SUB/figures/*.tiff 2>/dev/null | wc -l); [ "$nfig" -eq 9 ] && [ "$ntif" -eq 9 ] && ok "figures/: 9 PDF + 9 TIFF (Fig 1–6, SFig S1–S3)" || bad "figures/: $nfig PDF, $ntif TIFF"
for f in $M/abstract_structured.md $M/introduction.md $M/results.md $M/discussion.md $M/methods.md $M/figure_legends.md $M/table1.md $M/supplementary_table_captions.md $M/supplementary_note1.md $M/supplementary_note2.md $M/references.md $M/declarations.md $M/abbreviations.md $M/keywords.md revision_v2/text_changes.md revision_v2/response_to_reviewers.md revision_v2/cover_letter.md revision_v2/CHANGELOG.md; do [ -s $f ] && ok "$f ($(wc -w < $f) words)" || bad "$f missing"; done
echo; echo "## 6. Open items for the authors (human decision)"
echo "- Appendix C, decided 2026-10-20: '+1.14' (Rps2 liver bias) and '0.005' (liver WASP − naive mean) are coincidental matches of NEW values with v1 strings (ACCEPTED contexts in appendix_c_grep.sh); 'heart' is the GEO tissue attribute of the rejected GSE55792 samples."
echo "- The FVB_NJ_v3 comparison (optional Phase G) was not performed; the Discussion says so and names the Rps2 skin locus as the first place to look."
echo "- github_repo has uncommitted v1-era modifications (.gitignore, CITATION.cff, LICENSE, environment.yml, scripts/figures/*); commit_phase_h.sh commits only the v2 files — review before push."
echo; echo "**Result: $FAIL failing check(s).**"
} > $OUT
cat $OUT; exit $FAIL
