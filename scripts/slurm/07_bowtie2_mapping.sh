#!/bin/bash
#SBATCH --job-name=bowtie2_map
#SBATCH --output=logs/mapping_%A_%a.out
#SBATCH --error=logs/mapping_%A_%a.err
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --dependency=afterok:MEGAHIT_JOBID

# ============================================================
# 07_bowtie2_mapping.sh
# Map trimmed reads back to per-sample assemblies to generate
# coverage profiles required for depth-aware binning.
# ============================================================

# ---------- user-configurable paths ----------
CLEAN_DIR="/scratch/$USER/metagenomics/decontaminated"
ASSEMBLY_DIR="/scratch/$USER/metagenomics/assembly"
OUTDIR="/scratch/$USER/metagenomics/mapping"
SAMPLE_LIST="/scratch/$USER/metagenomics/metadata/sample_list.txt"
THREADS=$SLURM_CPUS_PER_TASK
# ---------------------------------------------

module purge
module load bowtie2/2.4.5
module load samtools/1.17

mkdir -p "$OUTDIR"
mkdir -p logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
CONTIG_FA="${ASSEMBLY_DIR}/${SAMPLE}/${SAMPLE}.contigs.fa"
INDEX_PREFIX="${OUTDIR}/${SAMPLE}_index"

echo "[$(date)] Indexing assembly for sample: $SAMPLE"
bowtie2-build --threads "$THREADS" "$CONTIG_FA" "$INDEX_PREFIX"

echo "[$(date)] Mapping reads for sample: $SAMPLE"
bowtie2 \
    -x "$INDEX_PREFIX" \
    -1 "${CLEAN_DIR}/${SAMPLE}_R1.fastq.gz" \
    -2 "${CLEAN_DIR}/${SAMPLE}_R2.fastq.gz" \
    --threads "$THREADS" \
    --no-unal \
    -S "${OUTDIR}/${SAMPLE}.sam" \
    2> "${OUTDIR}/${SAMPLE}_mapping.log"

echo "[$(date)] Converting and sorting BAM"
samtools view  -bS -@ "$THREADS" "${OUTDIR}/${SAMPLE}.sam" | \
samtools sort      -@ "$THREADS" -o "${OUTDIR}/${SAMPLE}_sorted.bam"
samtools index "${OUTDIR}/${SAMPLE}_sorted.bam"

# Remove large SAM to save disk space
rm "${OUTDIR}/${SAMPLE}.sam"

echo "[$(date)] Mapping finished for sample: $SAMPLE"
samtools flagstat "${OUTDIR}/${SAMPLE}_sorted.bam"
