#!/bin/bash
# Speed-up Step 2 switch-over at an animal boundary (run only after 06_shm_retest_b6gtf_index.md says ADOPT):
# graceful STOP -> wait for every worker and the scheduler to exit -> record which animals used the old index ->
# clear the flag -> restart scheduler.sh with STAR_MODE=shm (star_server.sh + rebuilt FVB index).
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; cd $PROJ; L=analysis_v2/locks
log(){ echo "[$(date '+%F %T')] switch $*"; }
grep -q '^\*\*ADOPT\*\*' revision_v2/speedup/06_shm_retest_b6gtf_index.md || { log "re-test verdict is not ADOPT; nothing done"; exit 1; }
log "touching STOP (workers finish their current step, scheduler reports HALTED)"; touch $L/STOP
while pgrep -f "[s]cheduler.sh" >/dev/null || pgrep -f "[r]un_sample.sh" >/dev/null || pgrep -x STAR >/dev/null; do sleep 60; done
log "scheduler and workers exited: $(tail -n 1 analysis_v2/logs/scheduler.log)"
# evidence: which fvb alignments were made with the old index + on-the-fly GTF (all fvb Log.out files so far)
{ echo -e "sample\tfvb_genomeDir\tmapping_time_GTF\tdate"; for f in analysis_v2/data/*/*/bam/*.fvb.Log.out; do s=$(basename $f .fvb.Log.out); g=$(grep -m1 '^genomeDir' $f | awk '{print $2}' | xargs basename); gt=$(grep -m1 '^sjdbGTFfile' $f | awk '{print $2}'); d=$(grep -m1 'started STAR run' $f | sed 's/ \.\.\.\.\. started STAR run//'); echo -e "$s\t$g\t${gt:--}\t$d"; done; } > revision_v2/speedup/08_fvb_index_per_animal.tsv
log "per-animal index record: $(tail -n +2 revision_v2/speedup/08_fvb_index_per_animal.tsv | wc -l) fvb alignments so far (revision_v2/speedup/08_fvb_index_per_animal.tsv)"
rm -f $L/STOP; log "STOP cleared; starting scheduler in STAR_MODE=shm"
STAR_MODE=shm nohup bash analysis_v2/pipeline/scheduler.sh >> analysis_v2/logs/scheduler.log 2>&1 &
sleep 5; log "scheduler pid $! ; $(tail -n 2 analysis_v2/logs/scheduler.log | tr '\n' ' ')"
