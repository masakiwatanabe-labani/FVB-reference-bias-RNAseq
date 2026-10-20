#!/bin/bash
# Download FASTQ for upcoming samples while run_sample.sh is computing on the current one (I/O overlap only).
# Same batch order as run_batch.sh (main_water -> fructose -> pyvt -> skin), same ENA FTP URLs and md5 values from
# the manifest. Each file goes to <sample>/fastq/.pf/, is md5-checked, and is moved into <sample>/fastq/ only if
# run_sample.sh has not started that sample yet. run_sample.sh then skips the download and repeats its own md5 check.
# At most AHEAD prefetched-but-not-started samples exist at any time (disk bound).
set -uo pipefail
PROJ=/usr/local/jupyter/FVB_bias; MAN=$PROJ/revision_v2/manifest/manifest.tsv; AHEAD=${AHEAD:-2}
log(){ echo "[$(date '+%F %T')] prefetch $*"; }
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
started(){ [ -e $(ddir $1)/logs/run_sample.log ]; }
waiting(){ n=0; for s in $LIST; do d=$(ddir $s); [ -e $d/fastq/prefetched ] && ! started $s && n=$((n+1)); done; echo $n; }
for s in $LIST; do
  D=$(ddir $s)
  started $s && continue
  [ -e $D/fastq/prefetched ] && continue
  while [ "$(waiting)" -ge $AHEAD ]; do sleep 300; done
  started $s && continue
  mkdir -p $D/fastq/.pf
  ftp=$(awk -F'\t' -v s=$s 'NR>1 && $1==s{print $22}' $MAN); md5=$(awk -F'\t' -v s=$s 'NR>1 && $1==s{print $23}' $MAN)
  IFS=';' read -ra U <<< "$ftp"; IFS=';' read -ra M <<< "$md5"
  ok=1; log "$s start (${#U[@]} files)"
  for i in "${!U[@]}"; do
    f=$(basename ${U[$i]})
    good=0
    for try in 1 2 3; do
      started $s && break 2
      rm -f $D/fastq/.pf/$f $D/fastq/.pf/$f.aria2
      /usr/bin/aria2c -q -x 4 -s 4 --retry-wait=10 --max-tries=8 -d $D/fastq/.pf -o $f "ftp://${U[$i]}" \
        && [ "$(md5sum $D/fastq/.pf/$f | cut -d' ' -f1)" = "${M[$i]}" ] && { good=1; break; }
      log "$s $f attempt $try failed (download or md5)"; sleep 120
    done
    [ $good -eq 1 ] || { ok=0; break; }
  done
  if [ $ok -eq 1 ] && ! started $s; then
    mv $D/fastq/.pf/*.fastq.gz $D/fastq/ && rmdir $D/fastq/.pf && touch $D/fastq/prefetched && log "$s ready ($(du -sh $D/fastq | cut -f1))"
  else
    rm -rf $D/fastq/.pf; log "$s skipped (run_sample started or download failed; run_sample.sh will download it itself)"
  fi
done
log "all samples handled"
