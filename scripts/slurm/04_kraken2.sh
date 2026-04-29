#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --output=logs/kraken2_%A_%a.out
#SBATCH --error=logs/kraken2_%A_%a.err
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --partition=highmem

# ============================================================
# 04_kraken2.sh
# k-mer based taxonomic classification with Kraken2.
# Uses the standard Kraken2 database (RefSeq archaea, bacteria,
# viral, plasmid, human, UniVec_Core).
# ============================================================

# ---------- user-configurable paths ----------
CLEAN_DIR="/scratch/$USER/metagenomics/decontaminated"
OUTDIR="/scratch/$USER/metagenomics/kraken2"
DB="/scratch/$USER/databases/kraken2_standard_20230314"
SAMPLE_LIST="/scratch/$USER/metagenomics/metadata/sample_list.txt"
CONFIDENCE=0.1    # classification confidence threshold
THREADS=$SLURM_CPUS_PER_TASK
# ---------------------------------------------

module purge
module load kraken2/2.1.3

mkdir -p "$OUTDIR"
mkdir -p logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "[$(date)] Kraken2 starting for sample: $SAMPLE"

kraken2 \
    --db            "$DB" \
    --threads       "$THREADS" \
    --confidence    "$CONFIDENCE" \
    --paired        \
    --gzip-compressed \
    --output        "${OUTDIR}/${SAMPLE}.kraken2" \
    --report        "${OUTDIR}/${SAMPLE}.k2report" \
    --report-minimizer-data \
    "${CLEAN_DIR}/${SAMPLE}_R1.fastq.gz" \
    "${CLEAN_DIR}/${SAMPLE}_R2.fastq.gz"

echo "[$(date)] Kraken2 finished for sample: $SAMPLE"
echo "  Classification rate:"
grep "sequences classified" "${OUTDIR}/${SAMPLE}.k2report" | head -2
