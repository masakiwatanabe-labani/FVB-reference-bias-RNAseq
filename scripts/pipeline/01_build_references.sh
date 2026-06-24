#!/usr/bin/env bash
# Pipeline step 1: Build reference genomes and personalised FVB/NJ reference
#
# Prerequisites:
#   - conda environment activated (see environment.yml)
#   - g2gtools Python 3 patch: in the g2gtools installation directory,
#     edit chain.py and replace all ".next()" calls with "next()" (Python 3).
#     Typically: sed -i 's/\.next()$/next(self)/g' $(python -c "import g2gtools; import os; print(os.path.dirname(g2gtools.__file__))")/chain.py
#
# Usage:
#   export FVB_PROJ=/path/to/working/directory
#   bash 01_build_references.sh
#
# Note: The MGP REL-2021 indel VCF carries FILTER="." for all records (not FILTER="PASS").
# The --pass flag in g2gtools vcf2vci therefore silently excludes all indels; the resulting
# personalised reference encodes SNP substitutions only. All bias estimates in the manuscript
# reflect this SNP-only correction.
# Stage 2: Download references and build STAR indices
# GRCm39 unified; B6 FASTA+GTF, FVB VCF, STAR indices (B6 + FVB personalized)
set -euo pipefail

PROJ="${FVB_PROJ:-$(cd "$(dirname "$0")/../.." && pwd)}"
# Tools are expected in PATH (see environment.yml for exact versions)
MM=""
REF=$PROJ/ref
THREADS=${THREADS:-8}


exec > >(tee -a "$LOG") 2>&1
echo "[$(date)] Stage 2 START"

mkdir -p $REF/{b6_fasta,b6_gtf,fvb_vcf,star_b6,star_fvb,g2g}

# ── 2.1 B6 GRCm39 primary FASTA (Ensembl release 110) ────────────────────────
FASTA_URL="https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
FASTA_GZ="$REF/b6_fasta/GRCm39.primary.fa.gz"
FASTA="$REF/b6_fasta/GRCm39.primary.fa"
if [ ! -f "$FASTA" ]; then
    echo "[$(date)] Downloading B6 FASTA..."
    aria2c -x 8 -s 8 --retry-wait=5 --max-tries=5 -d "$REF/b6_fasta" -o "GRCm39.primary.fa.gz" "$FASTA_URL"
    bgzip -d "$FASTA_GZ"
    echo "[$(date)] B6 FASTA downloaded: $FASTA"
else
    echo "[$(date)] B6 FASTA already exists, skipping."
fi

# ── 2.2 B6 GRCm39 GTF (Ensembl release 110) ──────────────────────────────────
GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz"
GTF_GZ="$REF/b6_gtf/GRCm39.110.gtf.gz"
GTF="$REF/b6_gtf/GRCm39.110.gtf"
if [ ! -f "$GTF" ]; then
    echo "[$(date)] Downloading GTF..."
    aria2c -x 4 -s 4 --retry-wait=5 --max-tries=5 -d "$REF/b6_gtf" -o "GRCm39.110.gtf.gz" "$GTF_URL"
    bgzip -d "$GTF_GZ"
    echo "[$(date)] GTF downloaded: $GTF"
else
    echo "[$(date)] GTF already exists, skipping."
fi

# ── 2.3 FVB/NJ VCF (Ensembl Mouse Genomes Project, GRCm39) ───────────────────
# Ensembl Mouse Genomes: mgp.v5 or sanger
# Primary source: https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/
FVB_VCF_SNP_URL="https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz"
FVB_VCF_INDEL_URL="https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_indels.vcf.gz"
FVB_SNP_ALL="$REF/fvb_vcf/mgp_REL2021_snps.vcf.gz"
FVB_INDEL_ALL="$REF/fvb_vcf/mgp_REL2021_indels.vcf.gz"
FVB_SNP="$REF/fvb_vcf/FVB_NJ.snps.GRCm39.vcf.gz"
FVB_INDEL="$REF/fvb_vcf/FVB_NJ.indels.GRCm39.vcf.gz"
FVB_COMBINED="$REF/fvb_vcf/FVB_NJ.combined.GRCm39.vcf.gz"

if [ ! -f "$FVB_COMBINED" ]; then
    echo "[$(date)] Downloading MGP VCF (all strains, SNPs)..."
    if [ ! -f "$FVB_SNP_ALL" ]; then
        aria2c -x 4 -s 4 --retry-wait=5 --max-tries=5 \
            -d "$REF/fvb_vcf" -o "mgp_REL2021_snps.vcf.gz" "$FVB_VCF_SNP_URL"
        aria2c -x 4 -s 4 --retry-wait=5 --max-tries=5 \
            -d "$REF/fvb_vcf" -o "mgp_REL2021_indels.vcf.gz" "$FVB_VCF_INDEL_URL"
        # Download index files too
        aria2c -x 2 -s 2 -d "$REF/fvb_vcf" -o "mgp_REL2021_snps.vcf.gz.tbi" "${FVB_VCF_SNP_URL}.tbi" 2>/dev/null || \
            bcftools index -t "$FVB_SNP_ALL"
        aria2c -x 2 -s 2 -d "$REF/fvb_vcf" -o "mgp_REL2021_indels.vcf.gz.tbi" "${FVB_VCF_INDEL_URL}.tbi" 2>/dev/null || \
            bcftools index -t "$FVB_INDEL_ALL"
    fi
    echo "[$(date)] Subsetting FVB/NJ from MGP VCF..."
    # The strain is labeled FVB_NJ in MGP
    bcftools view -s "FVB_NJ" -f "PASS" --min-ac 1 "$FVB_SNP_ALL" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bgzip > "$FVB_SNP"
    bcftools index -t "$FVB_SNP"

    bcftools view -s "FVB_NJ" -f "PASS" --min-ac 1 "$FVB_INDEL_ALL" | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
        bgzip > "$FVB_INDEL"
    bcftools index -t "$FVB_INDEL"

    # Combine SNP+indel for WASP and g2gtools
    bcftools concat --allow-overlaps -a "$FVB_SNP" "$FVB_INDEL" | \
        bcftools sort | bgzip > "$FVB_COMBINED"
    bcftools index -t "$FVB_COMBINED"
    echo "[$(date)] FVB VCF created: $FVB_COMBINED"
    # Count variants
    echo "[$(date)] FVB SNPs: $(bcftools view -H "$FVB_SNP" | wc -l)"
    echo "[$(date)] FVB indels: $(bcftools view -H "$FVB_INDEL" | wc -l)"
else
    echo "[$(date)] FVB VCF already exists, skipping."
fi

# ── 2.4 B6 STAR index (naive / WASP) ─────────────────────────────────────────
STAR_B6=$REF/star_b6
if [ ! -f "$STAR_B6/Genome" ]; then
    echo "[$(date)] Building B6 STAR index (~30 GB, will take ~1 hour)..."
    STAR --runMode genomeGenerate \
        --genomeDir "$STAR_B6" \
        --genomeFastaFiles "$FASTA" \
        --sjdbGTFfile "$GTF" \
        --genomeSAindexNbases 14 \
        --runThreadN $THREADS \
        --outFileNamePrefix "$STAR_B6/"
    echo "[$(date)] B6 STAR index complete."
else
    echo "[$(date)] B6 STAR index already exists, skipping."
fi

# ── 2.5 FVB personalized genome via g2gtools ─────────────────────────────────
G2G=$REF/g2g
FVB_VCI="$G2G/FVB_NJ.vci.gz"
FVB_FASTA="$G2G/FVB_NJ.fa"
FVB_GTF="$G2G/FVB_NJ.gtf"
STAR_FVB=$REF/star_fvb

if [ ! -f "$FVB_FASTA" ]; then
    echo "[$(date)] Building FVB personalized genome via g2gtools..."

    # Check g2gtools available
    if ! g2gtools --version &>/dev/null; then
        echo "[$(date)] g2gtools not in conda env, installing via pip..."
        pip install g2gtools --quiet
    fi

    # Step 1: vcf2vci
    echo "[$(date)] g2gtools vcf2vci..."
    g2gtools vcf2vci \
        -p $THREADS \
        --vcf "$FVB_COMBINED" \
        --fasta "$FASTA" \
        --strain FVB_NJ \
        -o "$FVB_VCI"

    # Step 2: patch (apply SNPs to B6 FASTA → FVB FASTA)
    echo "[$(date)] g2gtools patch..."
    g2gtools patch \
        -i "$FASTA" \
        -c "$FVB_VCI" \
        -o "${FVB_FASTA%.fa}_snp.fa"

    # Step 3: transform (apply indels)
    echo "[$(date)] g2gtools transform..."
    g2gtools transform \
        -i "${FVB_FASTA%.fa}_snp.fa" \
        -c "$FVB_VCI" \
        -o "$FVB_FASTA"

    echo "[$(date)] FVB personalized genome: $FVB_FASTA"
else
    echo "[$(date)] FVB genome already exists, skipping patch/transform."
fi

# Step 4: convert GTF to FVB coordinates
if [ ! -f "$FVB_GTF" ]; then
    echo "[$(date)] g2gtools convert GTF to FVB coords..."
    g2gtools convert \
        -i "$GTF" \
        -c "$FVB_VCI" \
        -o "$FVB_GTF"
    echo "[$(date)] FVB GTF: $FVB_GTF"
else
    echo "[$(date)] FVB GTF already exists."
fi

# Step 5: FVB STAR index
if [ ! -f "$STAR_FVB/Genome" ]; then
    echo "[$(date)] Building FVB STAR index..."
    mkdir -p "$STAR_FVB"
    STAR --runMode genomeGenerate \
        --genomeDir "$STAR_FVB" \
        --genomeFastaFiles "$FVB_FASTA" \
        --sjdbGTFfile "$FVB_GTF" \
        --genomeSAindexNbases 14 \
        --runThreadN $THREADS \
        --outFileNamePrefix "$STAR_FVB/"
    echo "[$(date)] FVB STAR index complete."
else
    echo "[$(date)] FVB STAR index already exists, skipping."
fi

echo "[$(date)] Stage 2 COMPLETE"
echo "  B6 FASTA: $FASTA"
echo "  B6 GTF:   $GTF"
echo "  FVB VCF:  $FVB_COMBINED"
echo "  B6 STAR:  $STAR_B6"
echo "  FVB STAR: $STAR_FVB"
echo "  FVB VCI:  $FVB_VCI"
