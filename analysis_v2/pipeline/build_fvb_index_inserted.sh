#!/bin/bash
# Step 2 variant (2026-10-17 18:12): materialise EXACTLY the genome production uses for the FVB arm — the old FVB index
# (ref/star_fvb: g2gtools FVB/NJ FASTA, SAindexNbases 14, overhang 100) with the B6 Ensembl 110 GTF inserted on the fly —
# by running STAR once with --sjdbInsertSave All, which writes the inserted genome as a complete genomeDir
# (_STARgenome). Motivation: a from-scratch rebuild with the same GTF (ref/star_fvb_b6gtf) reproduced every alignment,
# flag, CIGAR and count but renumbered the HI tag of 129 chrY multi-mapper records (06_shm_retest_b6gtf_index.md);
# the inserted genome keeps production's junction ordering and is expected to reproduce HI as well.
set -euo pipefail
PROJ=/usr/local/jupyter/FVB_bias; cd $PROJ; OUT=ref/star_fvb_b6gtf_ins; LOCK=analysis_v2/locks/star.lock
W=/tmp/claude-0/-usr-local-jupyter-FVB-bias/e900eed5-e9ea-4c20-a6d8-90c422bc83a4/scratchpad/step2
log(){ echo "[$(date '+%F %T')] build_ins $*"; }
for f in ref/g2g/GRCm39.FVB_NJ.fa ref/b6_gtf/GRCm39.110.gtf; do have=$(md5sum $f | cut -d' ' -f1); want=$(grep -F "$f" analysis_v2/references.md5 | awk '{print $1}'); [ "$have" = "$want" ] || { log "MD5 MISMATCH $f"; exit 3; }; log "md5 ok $f $have"; done
( cd ref/star_fvb && md5sum Genome SA SAindex genomeParameters.txt ) | diff - <(grep star_fvb/ analysis_v2/references.md5 | sed 's| ref/star_fvb/| |') > /dev/null && log "old index files match references.md5" || { log "old index md5 differs from references.md5"; exit 3; }
# BUG found 2026-10-18 09:05: under `set -o pipefail`, `head` closing the pipe gives zcat SIGPIPE -> non-zero pipeline -> `set -e` killed the
# script silently here at 18:33 on 2026-10-17 (no log line, tiny_2.fq never written, chained re-test never ran). Fixed with a subshell
# that tolerates SIGPIPE; the script was NOT re-run (cutoff 2026-10-18 09:00 JST passed, Step 2 abandoned per condition 5).
mkdir -p $OUT/tmp; (set +o pipefail; zcat $W/r1.fq.gz | head -n 4000 > $OUT/tmp/tiny_1.fq; zcat $W/r2.fq.gz | head -n 4000 > $OUT/tmp/tiny_2.fq)
log "waiting for STAR lock"
( /usr/bin/flock 9
  log "STAR lock acquired; inserting B6 GTF into ref/star_fvb with --sjdbInsertSave All (RAM avail $(free -g | awk 'NR==2{print $7}') G)"
  t0=$(date +%s)
  /usr/local/bin/STAR --runThreadN 8 --genomeDir $PROJ/ref/star_fvb --sjdbGTFfile $PROJ/ref/b6_gtf/GRCm39.110.gtf --sjdbInsertSave All \
    --readFilesIn $OUT/tmp/tiny_1.fq $OUT/tmp/tiny_2.fq --outSAMtype None --outFileNamePrefix $PROJ/$OUT/tmp/ins. > $OUT/tmp/ins.stdout 2>&1
  log "STAR insert run rc=$? in $(( ($(date +%s)-t0)/60 )) min"
) 9>$LOCK
ls $OUT/tmp/ins._STARgenome/ | tr '\n' ' ' | cut -c1-300; echo
for f in Genome SA SAindex chrName.txt chrLength.txt chrStart.txt chrNameLength.txt genomeParameters.txt sjdbList.out.tab sjdbInfo.txt exonInfo.tab geneInfo.tab transcriptInfo.tab exonGeTrInfo.tab sjdbList.fromGTF.out.tab; do [ -e $OUT/tmp/ins._STARgenome/$f ] && mv $OUT/tmp/ins._STARgenome/$f $OUT/; done
# the saved genomeParameters.txt lacks the Log of the insert run; keep the insert run's Log.out beside it
cp $OUT/tmp/ins.Log.out $OUT/insert_run.Log.out
log "junctions: $(wc -l < $OUT/sjdbList.out.tab) (old 286280, from-scratch rebuild 290071); chromosomes $(wc -l < $OUT/chrName.txt)"
( cd $OUT && md5sum Genome SA SAindex sjdbList.out.tab chrName.txt genomeParameters.txt > BUILD.md5 )
rm -rf $OUT/tmp; touch $OUT/done.build; log "done"
