#!/bin/bash
set -euo pipefail

# Input/Output paths
RMSK_TXT="/work/InternalMedicine/s184335/sc/workspace/sisi/scTE_genome_hsa/rmsk.txt"
OUT_DIR="/archive/InternalMedicine/Chung_lab/shared/sc/genome"
GTF_OUT="${OUT_DIR}/hg38_rmsk_TE.gtf"
LOC_OUT="${OUT_DIR}/hg38_rmsk_TE.gtf.locInd"
TELOCAL_SIF="/archive/InternalMedicine/Chung_lab/shared/sc/sig_image/telocal_latest.sif"

echo "=== 1. Converting rmsk.txt to GTF format ==="
awk 'BEGIN {FS="\t"; OFS="\t"} 
{
    chr = $6;
    start = $7 + 1; # Convert 0-based start to 1-based GTF coordinate
    end = $8;
    strand = $10;
    repName = $11;
    repClass = $12;
    repFamily = $13;

    # Filter out low-complexity and simple repeats
    if (repClass != "Simple_repeat" && repClass != "Low_complexity" && repClass != "Unknown") {
        attr = sprintf("gene_id \"%s\"; transcript_id \"%s\"; family_id \"%s\"; class_id \"%s\";", repName, repName, repFamily, repClass);
        print chr, "repeatmasker", "exon", start, end, ".", strand, ".", attr;
    }
}' "$RMSK_TXT" > "$GTF_OUT"

echo "[SUCCESS] Created GTF: $GTF_OUT"

echo "=== 2. Building TElocal index (.locInd) via Singularity ==="
singularity exec \
    --bind /archive/InternalMedicine/Chung_lab/shared/sc \
    "$TELOCAL_SIF" \
    TElocal_indexer --afile "$GTF_OUT" --itype TE --out "$LOC_OUT"

echo "[SUCCESS] Created TElocal Index: $LOC_OUT"
