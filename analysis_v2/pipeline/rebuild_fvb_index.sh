#!/bin/bash
# Speed-up Step 2 (approved 2026-10-17 09:00 with conditions): rebuild the FVB/NJ STAR index with the B6 Ensembl 110 GTF
# so that the index alone carries the junction set production obtains by mapping-time insertion (02_shm_equivalence.md).
# Condition 1: runs only while holding the STAR lock (genomeGenerate ~32 GB RAM beside a 27 GB STAR does not fit in 58 GB);
# --genomeSAindexNbases 14, --sjdbOverhang 100; input FASTA = current g2gtools FVB/NJ FASTA (md5 verified against
# analysis_v2/references.md5 before building); GTF = B6 Ensembl 110. The old index ref/star_fvb is not touched.
set -euo pipefail
PROJ=/usr/local/jupyter/FVB_bias; cd $PROJ
FA=ref/g2g/GRCm39.FVB_NJ.fa; GTF=ref/b6_gtf/GRCm39.110.gtf; OUT=ref/star_fvb_b6gtf; LOCK=analysis_v2/locks/star.lock
log(){ echo "[$(date '+%F %T')] rebuild $*"; }
log "md5 check of inputs"
for f in $FA $GTF; do
  have=$(md5sum $f | cut -d' ' -f1); want=$(grep -F "$f" analysis_v2/references.md5 | awk '{print $1}' | head -1)
  [ -n "$want" ] || { log "no recorded md5 for $f in analysis_v2/references.md5"; exit 2; }
  [ "$have" = "$want" ] || { log "MD5 MISMATCH for $f: have $have want $want -> abort"; exit 3; }
  log "md5 ok $f $have"
done
mkdir -p $OUT
log "waiting for STAR lock"
( /usr/bin/flock 9
  log "STAR lock acquired; genomeGenerate start (RAM avail $(free -g | awk 'NR==2{print $7}') G)"
  t0=$(date +%s)
  /usr/local/bin/STAR --runMode genomeGenerate --runThreadN 8 --genomeDir $PROJ/$OUT --genomeFastaFiles $PROJ/$FA \
    --genomeSAindexNbases 14 --sjdbGTFfile $PROJ/$GTF --sjdbOverhang 100 --outFileNamePrefix $PROJ/$OUT/ > $OUT/genomeGenerate.stdout 2>&1
  rc=$?; log "genomeGenerate rc=$rc in $(( ($(date +%s)-t0)/60 )) min"
  [ $rc -eq 0 ] || exit $rc
) 9>$LOCK
log "junctions: $(wc -l < $OUT/sjdbList.out.tab) (old index $(wc -l < ref/star_fvb/sjdbList.out.tab), B6 index $(wc -l < ref/star_b6/sjdbList.out.tab))"
log "chromosomes: $(wc -l < $OUT/chrName.txt) (old $(wc -l < ref/star_fvb/chrName.txt))"
( cd $OUT && md5sum Genome SA SAindex sjdbList.out.tab chrName.txt genomeParameters.txt > BUILD.md5 )
log "BUILD.md5 written; done"
touch $OUT/done.build
