#!/bin/bash
#SBATCH --job-name=host_decontam
#SBATCH --output=logs/decontam_%A_%a.out
#SBATCH --error=logs/decontam_%A_%a.err
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=standard

# ============================================================
# 03_host_decontam.sh
# Remove host (plant/eukaryote) reads by mapping against the
# host reference genome with Bowtie2.  Only unmapped read pairs
# (putative microbial reads) are retained.
# ============================================================

# ---------- user-configurable paths ----------
TRIMMED_DIR="/scratch/$USER/metagenomics/trimmed"
OUTDIR="/scratch/$USER/metagenomics/decontaminated"
HOST_INDEX="/scratch/$USER/databases/host_genome/host_bt2_index"
SAMPLE_LIST="/scratch/$USER/metagenomics/metadata/sample_list.txt"
THREADS=$SLURM_CPUS_PER_TASK
# ---------------------------------------------

module purge
module load bowtie2/2.4.5
module load samtools/1.17

mkdir -p "$OUTDIR"
mkdir -p logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "[$(date)] Host decontamination starting for sample: $SAMPLE"

# Map to host; capture unmapped pairs (--un-conc-gz)
bowtie2 \
    -x  "$HOST_INDEX" \
    -1  "${TRIMMED_DIR}/${SAMPLE}_R1_paired.fastq.gz" \
    -2  "${TRIMMED_DIR}/${SAMPLE}_R2_paired.fastq.gz" \
    --threads "$THREADS" \
    --sensitive \
    --un-conc-gz "${OUTDIR}/${SAMPLE}_clean_R%.fastq.gz" \
    -S /dev/null \
    2> "${OUTDIR}/${SAMPLE}_bowtie2_host.log"

# Rename Bowtie2 output (replaces % with 1 and 2)
mv "${OUTDIR}/${SAMPLE}_clean_R1.fastq.gz" \
   "${OUTDIR}/${SAMPLE}_R1.fastq.gz"
mv "${OUTDIR}/${SAMPLE}_clean_R2.fastq.gz" \
   "${OUTDIR}/${SAMPLE}_R2.fastq.gz"

echo "[$(date)] Host decontamination finished for sample: $SAMPLE"
echo "  Alignment summary:"
grep -E "reads|aligned" "${OUTDIR}/${SAMPLE}_bowtie2_host.log"
