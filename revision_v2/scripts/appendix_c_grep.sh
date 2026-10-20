#!/bin/bash
# Appendix C (v2 instruction document): strings from the frozen v1 manuscript that must not appear in any new text, figure or table.
# Scope (decided 2026-10-20, Phase D): PROSE = every .md under manuscript/revised_v2 (+ generated/*.md) and the Phase E letters in
# revision_v2/ -> all patterns. DATA = supplementary-table TSVs and figure *_data.tsv -> tissue words, verification sentences and the
# conditional gene names only; numeric patterns are not meaningful inside data tables (e.g. a mapping rate of 15.12 % or a
# per-sample value 0.005 are measurements, not the old manuscript's numbers). Two DATA sources are listed as EXEMPT and reported but
# do not fail the check: TableS5 (GEO series titles are external text) and the GEO tissue attribute "heart" of the excluded
# GSE55792 samples in TableS1a/S1b and in the Supplementary Fig. S1 row labels (true attribute of a rejected dataset; human decision).
# Exit 1 on any non-exempt hit. Usage: appendix_c_grep.sh
PROJ=/usr/local/jupyter/FVB_bias; cd $PROJ
PROSE=$(ls manuscript/revised_v2/*.md manuscript/revised_v2/generated/*.md 2>/dev/null)
LETTERS=$(ls revision_v2/response_to_reviewers.md revision_v2/cover_letter.md revision_v2/summary_of_changes.md 2>/dev/null)   # must quote the withdrawn v1 statements and numbers (disclosure): checked for tissue words and gene names only
DATA=$(ls manuscript/revised_v2/generated/*.tsv manuscript/revised_v2/supplementary_tables/*.tsv analysis_v2/figures/*_data.tsv 2>/dev/null)
TISSUE=("kidney" "cardiac" "cardiac muscle" "heart" "n = 10" "n = 12" "n = 7" "10 vs 12" "7 vs 6" "mesenteric adipose tissue samples. It included FVB/NJ mice (n = 6)")
NUM=("32-gene" "32 genes" "0.14%" "22,548" "144 unique genes" "0.64%" "+3.67" "+3.60" "−0.07" "+1.14" "+1.39" "+0.91" "+0.93" "+1.33" "+0.94" "+0.67" "+5.50" "+3.77" "+3.01" "7.89 SNPs/kb" "15.12" "60,818" "0.17%" "0.35%" "0.2% to 2.3%" "0.005")
VERIFY=("confirmed from GEO metadata" "recorded in the GEO sample attributes" "Verified from SRA")
GENES_CONDITIONAL=("Rsph3a" "Fau" "Pgam1" "Gm25939" "Atp5g2" "Mir6236" "Gm8909" "H2-Q5" "Gm55594" "Gm47283")
CAT=analysis_v2/results/catalog/catalog.primary.tsv
rc=0; nex=0
hit(){ # $1 pattern, $2 file list, $3 mode (FAIL|EXEMPT)
  local p="$1" files="$2" mode="$3" out
  out=$(grep -n -F -w -- "$p" $files 2>/dev/null | cut -c1-160)
  [ -z "$out" ] && return 0
  if [ "$mode" = FAIL ]; then echo "FORBIDDEN [$p]:"; echo "$out" | head -4; rc=1; else echo "EXEMPT (data) [$p]: $(echo "$out" | wc -l) line(s): $(echo "$out" | cut -d: -f1 | sort -u | tr '\n' ' ')"; nex=$((nex+1)); fi
}
allowed_gene(){ cut -f2 $CAT | grep -qx "$1"; }
# prose: everything
for p in "${TISSUE[@]}" "${NUM[@]}" "${VERIFY[@]}"; do hit "$p" "$PROSE" FAIL; done
for g in "${GENES_CONDITIONAL[@]}"; do allowed_gene "$g" || hit "$g" "$PROSE" FAIL; done
# letters: tissue words and gene names only
for p in "${TISSUE[@]}"; do hit "$p" "$LETTERS" FAIL; done
for g in "${GENES_CONDITIONAL[@]}"; do allowed_gene "$g" || hit "$g" "$LETTERS" FAIL; done
# data tables: tissue words / verification sentences / gene names; exempt sources reported separately
EXEMPT_FILES=$(ls manuscript/revised_v2/supplementary_tables/TableS5_survey_of_practice.tsv manuscript/revised_v2/supplementary_tables/TableS1a_manifest.tsv manuscript/revised_v2/supplementary_tables/TableS1b_identity_verification.tsv analysis_v2/figures/SFig1_GSE55792_haplotype_blocks_data.tsv 2>/dev/null)
OTHER_DATA=$(echo "$DATA" | tr " " "\n" | grep -vxF -f <(echo "$EXEMPT_FILES" | tr " " "\n"))
for p in "${TISSUE[@]}" "${VERIFY[@]}"; do hit "$p" "$OTHER_DATA" FAIL; [ -n "$EXEMPT_FILES" ] && hit "$p" "$EXEMPT_FILES" EXEMPT; done
for g in "${GENES_CONDITIONAL[@]}"; do allowed_gene "$g" || { hit "$g" "$OTHER_DATA" FAIL; hit "$g" "$EXEMPT_FILES" EXEMPT; }; done
echo "Appendix C grep: prose files $(echo $PROSE | wc -w), letters $(echo $LETTERS | wc -w), data files $(echo $DATA | wc -w); exempt hits $nex; $([ $rc -eq 0 ] && echo PASS || echo FAILED)"
exit $rc
