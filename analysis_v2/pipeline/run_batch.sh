#!/bin/bash
# Run run_sample.sh sequentially over manifest samples selected by role/design/tissue.
# Usage: run_batch.sh <batch>   batch = main_water | fructose | heart | pyvt | <sample_id list file>
# Stops immediately on a strain-identity mismatch (exit 2) or any other failure.
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; B=$1
case $B in
  main_water) SEL='$16=="strain_comparison" && $13=="Y" && $17=="Y" && $7!="heart"'; ORDER="adipose hypothalamus liver" ;;
  heart)      SEL='$16=="strain_comparison" && $13=="Y" && $7=="heart"'; ORDER="heart" ;;
  fructose)   SEL='$16=="strain_comparison" && $13=="Y" && $17=="N" && $18=="Y"'; ORDER="adipose hypothalamus liver" ;;
  skin)       SEL='$16=="strain_comparison" && $13=="Y" && $7=="skin"'; ORDER="skin" ;;
  pyvt)       SEL='$16=="negative_control_candidate" && $13=="Y"'; ORDER="mammary_tumour_PyVT" ;;
  *)          SEL=""; ;;
esac
if [ -n "$SEL" ]; then
  LIST=$(for t in $ORDER; do awk -F'\t' -v t=$t "NR>1 && $SEL && \$7==t {print \$1}" $PROJ/revision_v2/manifest/manifest.tsv; done)
else LIST=$(cat $B); fi
echo "[$(date '+%F %T')] batch $B: $(echo "$LIST" | wc -w) samples"
for s in $LIST; do
  d=$(awk -F'\t' -v s=$s 'NR>1 && $1==s{print $7}' $PROJ/revision_v2/manifest/manifest.tsv)
  [ -e $PROJ/analysis_v2/data/$d/$s/done.all ] && { echo "[$(date '+%F %T')] $s already done"; continue; }
  bash $PROJ/analysis_v2/pipeline/run_sample.sh $s; rc=$?
  if [ $rc -eq 2 ]; then echo "[$(date '+%F %T')] STOP: strain identity mismatch in $s (rule 0.4-4). Batch halted."; exit 2; fi
  if [ $rc -ne 0 ]; then echo "[$(date '+%F %T')] STOP: $s failed rc=$rc. Batch halted."; exit $rc; fi
  echo "[$(date '+%F %T')] free disk: $(df -h $PROJ | tail -1 | awk '{print $4}')"
done
echo "[$(date '+%F %T')] batch $B COMPLETE"
