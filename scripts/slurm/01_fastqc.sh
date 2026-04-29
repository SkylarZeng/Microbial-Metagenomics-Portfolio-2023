#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --error=logs/fastqc_%A_%a.err
#SBATCH --array=1-24          # one task per sample
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=standard

# ============================================================
# 01_fastqc.sh
# Run FastQC on raw paired-end FASTQ files (array job).
# Each SLURM array task processes one sample.
# ============================================================

# ---------- user-configurable paths ----------
FASTQ_DIR="/scratch/$USER/metagenomics/raw_reads"
OUTDIR="/scratch/$USER/metagenomics/qc/fastqc_raw"
SAMPLE_LIST="/scratch/$USER/metagenomics/metadata/sample_list.txt"
THREADS=$SLURM_CPUS_PER_TASK
# ---------------------------------------------

module purge
module load fastqc/0.11.9

mkdir -p "$OUTDIR"
mkdir -p logs

# Get the sample name for this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "[$(date)] FastQC starting for sample: $SAMPLE"

fastqc \
    --threads "$THREADS" \
    --outdir  "$OUTDIR" \
    "${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" \
    "${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

echo "[$(date)] FastQC finished for sample: $SAMPLE"
