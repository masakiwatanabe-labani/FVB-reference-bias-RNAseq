#!/bin/bash
# Phase C final run (decision 2026-10-16 #5: analysis_v2/results and analysis_v2/figures are written only by this single run,
# after all 83 animals finished; per-tissue dry runs validated the same scripts on scratch copies).
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; cd $PROJ
log(){ echo "[$(date '+%F %T')] final_run $*"; }
unset RES_DIR FIG_DIR
n=$(ls analysis_v2/data/*/*/done.all | wc -l); log "start: $n done.all markers (expected 83)"; [ "$n" -eq 83 ] || { log "abort: not all 83 done"; exit 1; }
[ -e analysis_v2/locks/STOP ] && { log "abort: STOP flag"; exit 1; }
log "merge_counts.py"; python3 analysis_v2/pipeline/merge_counts.py || { log "merge_counts FAILED"; exit 1; }
log "bias_analysis.py"; python3 -I analysis_v2/pipeline/bias_analysis.py || { log "bias_analysis FAILED"; exit 1; }
log "make_figures.sh"; bash revision_v2/scripts/make_figures.sh || { log "make_figures FAILED"; exit 1; }
( cd analysis_v2/results && find . -type f -name "*.tsv" | sort | xargs md5sum ) > analysis_v2/results/RESULTS.md5
( cd analysis_v2/figures && ls *.pdf *.tiff *_data.tsv 2>/dev/null | sort | xargs md5sum ) > analysis_v2/figures/FIGURES.md5
log "done"
