#!/bin/bash
# Shared-memory STAR server (speed-up Step 2, approved 2026-10-17 with conditions). ONE instance runs beside scheduler.sh.
# It owns the STAR lock for its whole lifetime (so no other STAR genome can be resident at the same time), keeps ONE
# genome (B6 or rebuilt FVB) in shared memory, and executes mapping requests that run_sample.sh (STAR_MODE=shm) drops
# into analysis_v2/locks/starq/<sid>.<arm>.req. Up to CONC requests of the resident arm run concurrently with 4 threads
# each; the genome is switched only when nothing is running, the resident arm's queue is empty and the other arm's
# queue is not (stage batches: B6 for a few animals' naive arm -> FVB for the same animals' FVB arm -> ...).
# Command line and outputs are exactly those of revision_v2/speedup/06_shm_retest_b6gtf_index.md (byte-identical to
# the old flock mode): LoadAndKeep, no mapping-time --sjdbGTFfile, BAM Unsorted + samtools sort -@4 -m 1G, written as
# <prefix>Aligned.sortedByCoord.out.bam so that run_sample.sh continues unchanged.
# Guards kept (condition 4): STOP flag -> pending requests are failed with "STOP" (workers pause, rc 3), running
# mappings finish, genome unloaded, server exits; a mapping starts only with MemAvailable >= RAM_JOB GB, a genome is
# loaded only with MemAvailable >= RAM_LOAD GB. Disk guard stays in scheduler.sh. QUIT file -> same orderly exit.
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; L=${STAR_LOCKDIR:-$PROJ/analysis_v2/locks}; Q=$L/starq; STOP=$L/STOP; QUIT=$Q/QUIT; mkdir -p $Q   # STAR_LOCKDIR/STAR_BIN/SAMTOOLS_BIN: overridable only for the stub simulation (revision_v2/speedup/07_server_simulation.log)
CONC=${CONC:-2}; RAM_JOB=${RAM_JOB:-6}; RAM_LOAD=${RAM_LOAD:-27}
STAR=${STAR_BIN:-/usr/local/bin/STAR}; SAMTOOLS=${SAMTOOLS_BIN:-/usr/local/bin/samtools}
VCF=$PROJ/ref/fvb_vcf/FVB_NJ.snps.GRCm39.vcf
declare -A GDIR=([naive]=$PROJ/ref/star_b6 [fvb]=$PROJ/ref/star_fvb_b6gtf)
declare -A EXTRA=([naive]="--waspOutputMode SAMtag --varVCFfile $VCF --outSAMattributes NH HI AS NM MD vW" [fvb]="")
log(){ echo "[$(date '+%F %T')] star_server $*"; }
avail_gb(){ awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo; }
G=none                      # resident genome arm: none|naive|fvb
declare -A RUN=()           # req name -> pid of the mapping+sort job
declare -A RUNARM=()        # req name -> arm
pending(){ ls $Q/*.$1.req 2>/dev/null | wc -l; }   # $1 arm
first_pending(){ ls -tr $Q/*.$1.req 2>/dev/null | head -1; }
load(){ local arm=$1 t0=$(date +%s)
  while [ "$(avail_gb)" -lt $RAM_LOAD ]; do log "waiting for RAM to load $arm genome (avail $(avail_gb) G < $RAM_LOAD G)"; sleep 60; done
  $STAR --genomeLoad LoadAndExit --genomeDir ${GDIR[$arm]} --outFileNamePrefix $Q/load.$arm. > $Q/load.$arm.log 2>&1 || { log "LoadAndExit $arm FAILED (see $Q/load.$arm.log) -> STOP"; touch $STOP; return 1; }
  G=$arm; log "$arm genome loaded in $(( $(date +%s)-t0 )) s (RAM avail $(avail_gb) G)"; }
unload(){ [ $G = none ] && return 0; $STAR --genomeLoad Remove --genomeDir ${GDIR[$G]} --outFileNamePrefix $Q/remove.$G. > $Q/remove.$G.log 2>&1; log "$G genome removed from shared memory"; G=none; }
run_job(){ # $1 <name>.running file (foreground inside a background subshell): STAR (4 threads, shared genome) + samtools sort; writes .done or .fail (+ .rc)
  local run=$1 base=${1%.running} sid arm reads prefix slog rc
  eval "$(grep -E '^(sid|arm|reads|prefix|slog)=' $run)"
  $STAR --runThreadN 4 --genomeLoad LoadAndKeep --genomeDir ${GDIR[$arm]} --readFilesIn $reads --readFilesCommand zcat \
    --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical ${EXTRA[$arm]} \
    --outSAMattrRGline ID:$sid SM:$sid PL:ILLUMINA --outFileNamePrefix $prefix > $slog 2>&1; rc=$?
  [ $rc -eq 0 ] || { echo "STAR rc=$rc" > $base.rc; rm -f $run; echo "STAR rc=$rc" > $base.fail; return; }
  $SAMTOOLS sort -@ 4 -m 1G -o ${prefix}Aligned.sortedByCoord.out.bam ${prefix}Aligned.out.bam >> $slog 2>&1; rc=$?
  [ $rc -eq 0 ] || { echo "samtools sort rc=$rc" > $base.rc; rm -f $run; echo "samtools sort rc=$rc" > $base.fail; return; }
  rm -f ${prefix}Aligned.out.bam; echo ok > $base.rc; rm -f $run; touch $base.done
}
reap(){ local r; for r in "${!RUN[@]}"; do
  if ! kill -0 ${RUN[$r]} 2>/dev/null; then wait ${RUN[$r]}; log "job $r (${RUNARM[$r]}) finished: $(cat $Q/$r.rc 2>/dev/null || echo 'no rc file')"; rm -f $Q/$r.rc; unset RUN[$r] RUNARM[$r]; fi; done; }
fail_pending(){ local f; for f in $Q/*.req; do [ -e "$f" ] || continue; echo STOP > ${f%.req}.fail; rm -f $f; log "pending $(basename $f) failed with STOP (STOP flag or QUIT: worker pauses, rc 3)"; done; }

exec 9>$L/star.lock
log "waiting for STAR lock (pid $$, CONC=$CONC, RAM_JOB=$RAM_JOB G, RAM_LOAD=$RAM_LOAD G)"
/usr/bin/flock 9; log "STAR lock acquired; serving $Q"
rm -f $QUIT
for f in $Q/*.running; do [ -e "$f" ] && { mv $f ${f%.running}.req; log "stale $(basename $f) from a previous server re-queued"; }; done
while :; do
  reap
  if [ -e $STOP ] || [ -e $QUIT ]; then
    fail_pending
    if [ ${#RUN[@]} -eq 0 ]; then unload; log "exit ($([ -e $STOP ] && echo STOP || echo QUIT))"; exit 0; fi
    sleep 20; continue
  fi
  other=$([ $G = naive ] && echo fvb || echo naive)
  if [ ${#RUN[@]} -eq 0 ]; then
    if [ $G = none ]; then
      [ "$(pending naive)" -gt 0 ] && load naive || { [ "$(pending fvb)" -gt 0 ] && load fvb; }
    elif [ "$(pending $G)" -eq 0 ] && [ "$(pending $other)" -gt 0 ]; then
      unload; load $other
    fi
  fi
  if [ $G != none ] && [ ${#RUN[@]} -lt $CONC ] && [ "$(pending $G)" -gt 0 ] && [ "$(avail_gb)" -ge $RAM_JOB ]; then
    req=$(first_pending $G); name=$(basename ${req%.req})
    if [ -n "$req" ] && [ -z "${RUN[$name]+x}" ]; then
      mv $req $Q/$name.running   # a request is pending only while its .req exists; .running marks it as taken (simulation bug 2026-10-17 10:49: the running job was re-picked / failed by STOP)
      ( run_job $Q/$name.running ) & RUN[$name]=$!; RUNARM[$name]=$G; log "job $name started (pid $!, running ${#RUN[@]}, RAM avail $(avail_gb) G)"
      sleep 5; continue
    fi
  fi
  sleep 15
done
