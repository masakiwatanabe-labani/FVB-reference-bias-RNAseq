#!/bin/bash
# Stage-pipelined scheduler (speed-up Step 1). Up to MAXJ animals are "in progress" at once, each handled end to end by
# run_sample.sh; the two STAR runs per animal are serialised through flock on analysis_v2/locks/star.lock (inside
# run_sample.sh), so at most ONE STAR process exists at any time and everything else (download, fastp, identity QC,
# WASP filter, featureCounts, cleanup) overlaps with it. Same sample order as run_batch.sh main_water -> fructose ->
# pyvt -> skin, same manifest selection, same commands and parameters (only thread counts of non-STAR tools differ:
# fastp 2, samtools 2, featureCounts 2 x 2 concurrent; GATK is single-threaded).
# Stop rules: a strain-identity mismatch (rc 2) or any other failure creates analysis_v2/locks/STOP; running workers
# finish their current step and exit (rc 3), no new work starts, and the scheduler reports HALTED. `touch
# analysis_v2/locks/STOP` by hand gives the same graceful stop (used for the Step 2 switch-over). Remove the flag and
# restart to resume; done-markers make every step resumable.
# Guards: new workers start only when free disk >= 60 GB and MemAvailable >= 6 GB (run_sample.sh also waits for RAM
# before each non-STAR step).
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; MAN=$PROJ/revision_v2/manifest/manifest.tsv
export STAR_MODE=${STAR_MODE:-lock}   # shm = speed-up Step 2: star_server.sh keeps one genome in shared memory and maps 2 animals x 4 threads
MAXJ=${MAXJ:-$([ "$STAR_MODE" = shm ] && echo 5 || echo 3)}   # shm: 5 in progress so that naive/fvb requests form batches of 2-3 per genome switch
L=$PROJ/analysis_v2/locks; STOP=$L/STOP; mkdir -p $L
export FASTP_T=2 SAM_T=2 FC_T=2 FC_PAR=2 THREADS=8
log(){ echo "[$(date '+%F %T')] scheduler $*"; }
sel(){ case $1 in
  main_water)   S='$16=="strain_comparison" && $13=="Y" && $17=="Y" && $7!="heart" && $7!="skin"'; O="adipose hypothalamus liver" ;;
  fructose)     S='$16=="strain_comparison" && $13=="Y" && $17=="N" && $18=="Y" && $7!="skin"'; O="adipose hypothalamus liver" ;;
  skin_control) S='$16=="strain_comparison" && $13=="Y" && $7=="skin" && $17=="Y"'; O="skin" ;;
  skin_imq)     S='$16=="strain_comparison" && $13=="Y" && $7=="skin" && $17=="N"'; O="skin" ;;
  pyvt)         S='$16=="negative_control_candidate" && $13=="Y"'; O="mammary_tumour_PyVT" ;;
  esac
  for t in $O; do awk -F'\t' -v t=$t "NR>1 && $S && \$7==t {print \$1}" $MAN; done; }
LIST=$(for b in main_water fructose skin_control skin_imq pyvt; do sel $b; done)   # order decided 2026-10-16: main analysis first, sensitivity (IMQ) and negative control (PyVT) last
ddir(){ echo $PROJ/analysis_v2/data/$(awk -F'\t' -v s=$1 'NR>1 && $1==s{print $7}' $MAN)/$1; }
free_gb(){ df -BG $PROJ | tail -1 | awk '{print $4+0}'; }
avail_gb(){ awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo; }
declare -A PID=()   # "=()" required: under set -u a declared-but-empty array is "unbound" (bug found 2026-10-16 18:30)
reap(){ local s rc; for s in "${!PID[@]}"; do   # "local s": without it reap overwrote the caller's loop variable (bug found 2026-10-16 18:33, duplicate workers)
  if ! kill -0 ${PID[$s]} 2>/dev/null || [ "$(awk '/^State/{print $2}' /proc/${PID[$s]}/status 2>/dev/null)" = "Z" ]; then wait ${PID[$s]}; rc=$?; unset PID[$s]
    case $rc in
      0) log "$s finished rc=0 (free disk $(free_gb)G)";;
      2) log "STOP: strain identity mismatch in $s (rule 0.4-4); other workers finish their current step"; touch $STOP;;
      3) log "$s paused by STOP flag";;
      *) log "STOP: $s failed rc=$rc; other workers finish their current step"; touch $STOP;;
    esac
  fi; done; }
# do not start while the last animal of the old sequential batch (its STAR is not under the lock) is still running
while pgrep -f "[r]un_sample.sh" >/dev/null || pgrep -x STAR >/dev/null; do sleep 60; done
log "start: $(echo "$LIST" | wc -w) samples in list, MAXJ=$MAXJ, STAR_MODE=$STAR_MODE, $(for s in $LIST; do [ -e $(ddir $s)/done.all ] && echo 1; done | wc -l) already done"
SRV=""
if [ "$STAR_MODE" = shm ]; then
  [ -e $PROJ/ref/star_fvb_b6gtf/done.build ] || { log "STAR_MODE=shm but ref/star_fvb_b6gtf is not built; abort"; exit 1; }
  rm -f $L/starq/*.req $L/starq/*.done $L/starq/*.fail $L/starq/*.rc $L/starq/QUIT 2>/dev/null
  bash $PROJ/analysis_v2/pipeline/star_server.sh >> $PROJ/analysis_v2/logs/star_server.log 2>&1 & SRV=$!; log "star_server.sh started (pid $SRV)"
fi
for s in $LIST; do
  D=$(ddir $s)
  [ -e $D/done.all ] && continue
  [ -n "${PID[$s]+x}" ] && { log "BUG: $s already running (pid ${PID[$s]}); skipped"; continue; }
  while :; do
    reap
    [ -e $STOP ] && break
    if [ ${#PID[@]} -lt $MAXJ ] && [ "$(free_gb)" -ge 60 ] && [ "$(avail_gb)" -ge 6 ]; then break; fi
    sleep 60
  done
  [ -e $STOP ] && break
  mkdir -p $D/logs
  bash $PROJ/analysis_v2/pipeline/run_sample.sh $s all >> $D/logs/worker.log 2>&1 &
  PID[$s]=$!; log "$s started (worker pid $!, in progress: ${#PID[@]}, free disk $(free_gb)G, RAM avail $(avail_gb)G)"
  sleep 5
done
while [ ${#PID[@]} -gt 0 ]; do reap; sleep 30; done
if [ -n "$SRV" ]; then touch $L/starq/QUIT; while kill -0 $SRV 2>/dev/null; do sleep 10; done; log "star_server.sh exited (genome unloaded)"; fi
if [ -e $STOP ]; then log "HALTED (STOP flag at $STOP). Fix, remove the flag and rerun scheduler.sh to resume."; exit 2; fi
log "COMPLETE: all samples in list done"
