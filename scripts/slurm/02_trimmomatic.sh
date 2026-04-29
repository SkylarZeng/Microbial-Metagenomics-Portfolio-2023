#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --output=logs/trimmomatic_%A_%a.out
#SBATCH --error=logs/trimmomatic_%A_%a.err
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=standard

# ============================================================
# 02_trimmomatic.sh
# Adapter trimming and quality filtering with Trimmomatic.
# Removes Illumina TruSeq adapters; trims low-quality bases;
# discards reads shorter than 50 bp.
# ============================================================

# ---------- user-configurable paths ----------
FASTQ_DIR="/scratch/$USER/metagenomics/raw_reads"
OUTDIR="/scratch/$USER/metagenomics/trimmed"
SAMPLE_LIST="/scratch/$USER/metagenomics/metadata/sample_list.txt"
ADAPTERS="/opt/software/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa"
THREADS=$SLURM_CPUS_PER_TASK
# ---------------------------------------------

module purge
module load trimmomatic/0.39
module load fastqc/0.11.9

mkdir -p "$OUTDIR"
mkdir -p logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "[$(date)] Trimmomatic starting for sample: $SAMPLE"

trimmomatic PE \
    -threads "$THREADS" \
    -phred33 \
    "${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" \
    "${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz" \
    "${OUTDIR}/${SAMPLE}_R1_paired.fastq.gz"   \
    "${OUTDIR}/${SAMPLE}_R1_unpaired.fastq.gz" \
    "${OUTDIR}/${SAMPLE}_R2_paired.fastq.gz"   \
    "${OUTDIR}/${SAMPLE}_R2_unpaired.fastq.gz" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10:8:true \
    LEADING:3        \
    TRAILING:3       \
    SLIDINGWINDOW:4:20 \
    MINLEN:50

# Post-trim QC
fastqc \
    --threads "$THREADS" \
    --outdir  "${OUTDIR}/fastqc_post" \
    "${OUTDIR}/${SAMPLE}_R1_paired.fastq.gz" \
    "${OUTDIR}/${SAMPLE}_R2_paired.fastq.gz"

echo "[$(date)] Trimmomatic finished for sample: $SAMPLE"
