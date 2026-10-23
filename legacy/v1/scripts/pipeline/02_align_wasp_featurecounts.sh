#!/usr/bin/env bash
# Pipeline step 2: Three-arm alignment (naive / WASP / FVB personalised reference)
#                  followed by featureCounts (unique-only; no -M flag)
#
# For each dataset, this script:
#   1. Downloads FASTQs from ENA FTP (aria2c) or uses pre-downloaded files
#   2. Trims with fastp
#   3. Aligns to B6 GRCm39 (naive arm, WASP-tagged)
#   4. Filters WASP vW=2 reads to produce WASP arm
#   5. Aligns to FVB/NJ personalised reference (FVB arm)
#   6. Runs featureCounts (unique-only, -s 0) for all three arms
#
# Primary featureCounts quantification uses DEFAULT settings (no -M flag).
# Multi-mapping reads appear in Unassigned_MultiMapping (~7-10M per sample).
# The -M --fraction mode is ONLY used in step 5 (multimapping_sensitivity.py).
#
# Usage:
#   export FVB_PROJ=/path/to/working/directory
#   export THREADS=8
#   bash 02_align_wasp_featurecounts.sh
#
# Note: Download FASTQs from SRA before running, or let the script download them.
# SRA Toolkit: prefetch + fasterq-dump; ENA FTP: aria2c (URLs constructed below).

set -euo pipefail
PROJ="${FVB_PROJ:-$(cd "$(dirname "$0")/../.." && pwd)}"
THREADS="${THREADS:-8}"
REF="$PROJ/ref"
RESULTS="$PROJ/results/strain_comp"
FASTQ_DIR="$PROJ/data/fastq"
mkdir -p "$RESULTS" "$FASTQ_DIR" "$PROJ/logs"

STAR_B6="$REF/star_b6"
STAR_FVB="$REF/star_fvb"
GTF_B6="$REF/b6_gtf/GRCm39.110.gtf"
FVB_VCF_SNP="$REF/fvb_vcf/FVB_NJ.snps.GRCm39.vcf"   # uncompressed; STAR requires plain VCF

echo "[$(date)] Step 2 START"

# ── Sample accession lists ────────────────────────────────────────────────────
# Verified from SRA for each GEO series

declare -A FVB_SRRS=(
    ["GSE123875"]="SRR8321484 SRR8321485 SRR8321486 SRR8321487 SRR8321488 SRR8321489"
    ["GSE123893"]="SRR8321791 SRR8321792 SRR8321793 SRR8321794 SRR8321795 SRR8321796"
    ["GSE123894"]="SRR8322104 SRR8322105 SRR8322106 SRR8322107 SRR8322108 SRR8322109"
    ["GSE135230"]="SRR9879444 SRR9879445 SRR9879446 SRR9879447 SRR9879448 SRR9879449 SRR9879450 SRR9879451 SRR9879452 SRR9879453"
    ["GSE175625"]="SRR14664845 SRR14664846 SRR14664847 SRR14664848 SRR14664849 SRR14664850 SRR14664851"
)
declare -A B6_SRRS=(
    ["GSE123875"]="SRR8321439 SRR8321440 SRR8321441 SRR8321442 SRR8321443 SRR8321444"
    ["GSE123893"]="SRR8321732 SRR8321733 SRR8321734 SRR8321735 SRR8321736 SRR8321737"
    ["GSE123894"]="SRR8322030 SRR8322031 SRR8322032 SRR8322033 SRR8322034 SRR8322035"
    ["GSE135230"]="SRR9879454 SRR9879455 SRR9879456 SRR9879457 SRR9879458 SRR9879459 SRR9879460 SRR9879461 SRR9879462 SRR9879463 SRR9879464 SRR9879465"
    ["GSE175625"]="SRR14664793 SRR14664794 SRR14664795 SRR14664796 SRR14664797 SRR14664798"
)
# PE = paired-end; SE = single-end
declare -A LIBRARY=(
    ["GSE123875"]="PE" ["GSE123893"]="PE" ["GSE123894"]="PE"
    ["GSE135230"]="SE" ["GSE175625"]="SE"
)

# ── Helper: download one SRR from ENA FTP ────────────────────────────────────
download_srr() {
    local srr=$1 outdir=$2 pe=$3
    local prefix="${srr:0:6}"
    local url_base="https://ftp.sra.ebi.ac.uk/vol1/fastq/$prefix"
    [[ ${#srr} -eq 10 ]] && url_base="$url_base/00${srr: -1}/$srr" || url_base="$url_base/$srr"
    mkdir -p "$outdir"
    if [[ "$pe" == "PE" ]]; then
        [[ -f "$outdir/${srr}_1.fastq.gz" ]] && return 0
        aria2c -x 4 -s 4 --retry-wait=10 --max-tries=5 \
            "${url_base}/${srr}_1.fastq.gz" "${url_base}/${srr}_2.fastq.gz" -d "$outdir"
    else
        [[ -f "$outdir/${srr}.fastq.gz" ]] && return 0
        aria2c -x 4 -s 4 --retry-wait=10 --max-tries=5 \
            "${url_base}/${srr}.fastq.gz" -d "$outdir"
    fi
}

# ── Process each dataset ──────────────────────────────────────────────────────
for GSE in GSE123875 GSE123893 GSE123894 GSE135230 GSE175625; do
    LIB="${LIBRARY[$GSE]}"
    PE_FLAG=""; [[ "$LIB" == "PE" ]] && PE_FLAG="-p"
    echo "[$(date)] === $GSE ($LIB) ==="

    ALL_SRRS="${FVB_SRRS[$GSE]} ${B6_SRRS[$GSE]}"
    BAM_LIST_NAIVE=(); BAM_LIST_FVB=()

    for SRR in $ALL_SRRS; do
        SDIR="$FASTQ_DIR/$GSE"
        download_srr "$SRR" "$SDIR" "$LIB"

        ADIR="$RESULTS/bam/${GSE}/${SRR}"
        mkdir -p "$ADIR"

        # ── Trim ──
        TRIM_DIR="$PROJ/tmp/trim/${GSE}/${SRR}"
        mkdir -p "$TRIM_DIR"
        if [[ "$LIB" == "PE" ]]; then
            TRIM_R1="$TRIM_DIR/${SRR}_1.fq.gz"
            TRIM_R2="$TRIM_DIR/${SRR}_2.fq.gz"
            [[ -f "$TRIM_R1" ]] || fastp \
                -i "$SDIR/${SRR}_1.fastq.gz" -I "$SDIR/${SRR}_2.fastq.gz" \
                -o "$TRIM_R1" -O "$TRIM_R2" \
                --detect_adapter_for_pe -w 4 -j /dev/null -h /dev/null
            READS="$TRIM_R1 $TRIM_R2"
        else
            TRIM_R1="$TRIM_DIR/${SRR}.fq.gz"
            [[ -f "$TRIM_R1" ]] || fastp \
                -i "$SDIR/${SRR}.fastq.gz" -o "$TRIM_R1" \
                -w 4 -j /dev/null -h /dev/null
            READS="$TRIM_R1"
        fi

        # ── Naive alignment (B6 ref, WASP-tagged) ──
        NAIVE_BAM="$ADIR/${SRR}.naive.bam"
        if [[ ! -f "$NAIVE_BAM" ]]; then
            STAR --runThreadN "$THREADS" \
                --genomeDir "$STAR_B6" --sjdbGTFfile "$GTF_B6" \
                --readFilesIn $READS --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMstrandField intronMotif \
                --outFilterIntronMotifs RemoveNoncanonical \
                --waspOutputMode SAMtag --varVCFfile "$FVB_VCF_SNP" \
                --outFileNamePrefix "$ADIR/${SRR}.naive." \
                --outSAMattributes NH HI AS NM MD vW
            mv "$ADIR/${SRR}.naive.Aligned.sortedByCoord.out.bam" "$NAIVE_BAM"
            samtools index "$NAIVE_BAM"
        fi

        # ── WASP arm: exclude vW=2 reads ──
        WASP_BAM="$ADIR/${SRR}.wasp.bam"
        [[ -f "$WASP_BAM" ]] || {
            samtools view -b -e '![vW] || [vW]==1' "$NAIVE_BAM" > "$WASP_BAM"
            samtools index "$WASP_BAM"
        }

        # ── FVB reference alignment ──
        FVB_BAM="$ADIR/${SRR}.fvb.bam"
        if [[ ! -f "$FVB_BAM" ]]; then
            STAR --runThreadN "$THREADS" \
                --genomeDir "$STAR_FVB" --sjdbGTFfile "$GTF_B6" \
                --readFilesIn $READS --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMstrandField intronMotif \
                --outFilterIntronMotifs RemoveNoncanonical \
                --outFileNamePrefix "$ADIR/${SRR}.fvb."
            mv "$ADIR/${SRR}.fvb.Aligned.sortedByCoord.out.bam" "$FVB_BAM"
            samtools index "$FVB_BAM"
        fi

        BAM_LIST_NAIVE+=("$NAIVE_BAM")
        BAM_LIST_FVB+=("$FVB_BAM")

        # Delete trimmed FASTQs to save disk space
        rm -rf "$TRIM_DIR"
    done

    # ── featureCounts (unique-only; default: no -M flag) ──
    # Three outputs: naive, wasp, fvb_b6gtf
    echo "[$(date)] featureCounts $GSE..."

    # Collect WASP BAMs
    BAM_LIST_WASP=()
    for SRR in $ALL_SRRS; do
        BAM_LIST_WASP+=("$RESULTS/bam/${GSE}/${SRR}/${SRR}.wasp.bam")
    done

    for ARM in naive wasp fvb_b6gtf; do
        OUT_COUNTS="$RESULTS/${GSE}.${ARM}.counts.txt"
        [[ -f "$OUT_COUNTS" ]] && continue
        if [[ "$ARM" == "naive" ]]; then
            BAMS=("${BAM_LIST_NAIVE[@]}")
        elif [[ "$ARM" == "wasp" ]]; then
            BAMS=("${BAM_LIST_WASP[@]}")
        else
            BAMS=("${BAM_LIST_FVB[@]}")
        fi
        featureCounts -a "$GTF_B6" -o "$OUT_COUNTS" \
            -s 0 -T "$THREADS" $PE_FLAG \
            "${BAMS[@]}"
    done

    echo "[$(date)] $GSE DONE — deleting raw FASTQs"
    rm -rf "$FASTQ_DIR/$GSE"
done

echo "[$(date)] Step 2 COMPLETE"
echo "Count files: $RESULTS/"
