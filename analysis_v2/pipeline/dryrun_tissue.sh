#!/bin/bash
# Per-tissue dry run (decision 2026-10-16 #5): once every include=Y animal of a tissue is done, run merge_counts.py and
# bias_analysis.py with RES_DIR pointing to a scratch directory. analysis_v2/results/ is untouched; the final run is
# done once, after all tissues. Usage: dryrun_tissue.sh <tag>  (tag names the scratch dir, e.g. hypothalamus)
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; TAG=$1
OUT=/tmp/claude-0/-usr-local-jupyter-FVB-bias/e900eed5-e9ea-4c20-a6d8-90c422bc83a4/scratchpad/dryrun_$TAG; rm -rf $OUT; mkdir -p $OUT
export RES_DIR=$OUT
cd $PROJ
echo "[$(date '+%F %T')] dry run $TAG -> $OUT"
python3 analysis_v2/pipeline/merge_counts.py && python3 -I analysis_v2/pipeline/bias_analysis.py
echo "[$(date '+%F %T')] dry run $TAG rc=$?"
