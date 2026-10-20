#!/bin/bash
# Build all v2 figures from analysis_v2/results (or RES_DIR) into analysis_v2/figures (or FIG_DIR).
# Usage: make_figures.sh [RES_DIR] [FIG_DIR]   e.g. a partial-data test: make_figures.sh scratch/dryrun_partial scratch/figs_partial
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; S=$PROJ/revision_v2/scripts
[ $# -ge 1 ] && export RES_DIR=$1; [ $# -ge 2 ] && export FIG_DIR=$2
cd $S; rc=0
LIST="fig1_survey fig2_f1_ase fig3_strain_comp fig4_predictors fig5_catalog fig6_worked_example sfig1_multimap"; [ -z "${FIG_DIR:-}" ] && LIST="$LIST sfig_gse55792_blocks"   # the GSE55792 figure does not depend on RES_DIR
for f in $LIST; do
  echo "== $f"; python3 $f.py || { echo "!! $f failed"; rc=1; }
done
exit $rc
