#!/bin/bash
# Assemble manuscript/revised_v2/manuscript_combined.md from the section files (same order as the v1 combined file) and run the
# two Phase F checks: figure citation order (revision/scripts/check_figure_order.py) and the Appendix C forbidden-string grep
# (revision_v2/scripts/appendix_c_grep.sh). Usage: combine_manuscript.sh
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; M=$PROJ/manuscript/revised_v2; cd $M
{
  for f in abstract_structured introduction results discussion methods references figure_legends table1 supplementary_table_captions supplementary_note1 supplementary_note2; do
    cat $f.md; echo; echo; echo "---"; echo; echo
  done
} > manuscript_combined.md
echo "combined: $(wc -l < manuscript_combined.md) lines, $(wc -w < manuscript_combined.md) words"
# Phase T: the combined file keeps the v1 section order (Methods after Discussion); the citation order is checked in the order of the
# LAR article (Background → Methods → Results → Discussion, as assembled by assemble_lar.py)
echo "== figure citation order (main text, LAR section order)"; python3 $PROJ/revision/scripts/check_figure_order.py manuscript_combined.md --legend-check --section-order "Abstract,Background,Methods,Results,Discussion"
echo "== Appendix C grep"; bash $PROJ/revision_v2/scripts/appendix_c_grep.sh
