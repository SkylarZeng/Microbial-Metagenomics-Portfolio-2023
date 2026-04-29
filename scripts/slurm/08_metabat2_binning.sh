#!/bin/bash
#SBATCH --job-name=metabat2
#SBATCH --output=logs/metabat2_%A_%a.out
#SBATCH --error=logs/metabat2_%A_%a.err
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=standard
#SBATCH --dependency=afterok:MAPPING_JOBID

# ============================================================
# 08_metabat2_binning.sh
# Genome binning with MetaBAT2 using per-contig coverage
# depth profiles computed from sorted BAM files.
# ============================================================

# ---------- user-configurable paths ----------
ASSEMBLY_DIR="/scratch/$USER/metagenomics/assembly"
MAPPING_DIR="/scratch/$USER/metagenomics/mapping"
OUTDIR="/scratch/$USER/metagenomics/bins/metabat2"
SAMPLE_LIST="/scratch/$USER/metagenomics/metadata/sample_list.txt"
MIN_CONTIG_LEN=2500    # MetaBAT2 default; raise to 3000 for cleaner bins
THREADS=$SLURM_CPUS_PER_TASK
# ---------------------------------------------

module purge
module load metabat2/2.15

mkdir -p "$OUTDIR"
mkdir -p logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
CONTIG_FA="${ASSEMBLY_DIR}/${SAMPLE}/${SAMPLE}.contigs.fa"
BAM="${MAPPING_DIR}/${SAMPLE}_sorted.bam"
DEPTH_FILE="${OUTDIR}/${SAMPLE}_depth.txt"
BIN_PREFIX="${OUTDIR}/${SAMPLE}/${SAMPLE}_bin"

echo "[$(date)] Computing contig depth for sample: $SAMPLE"
jgi_summarize_bam_contig_depths \
    --outputDepth "$DEPTH_FILE" \
    "$BAM"

echo "[$(date)] MetaBAT2 binning for sample: $SAMPLE"
mkdir -p "${OUTDIR}/${SAMPLE}"

metabat2 \
    -i "$CONTIG_FA" \
    -a "$DEPTH_FILE" \
    -o "$BIN_PREFIX" \
    -m "$MIN_CONTIG_LEN" \
    -t "$THREADS" \
    --unbinned \
    --saveCls \
    --verbose

echo "[$(date)] Binning finished for sample: $SAMPLE"
echo "  Bin count: $(ls ${OUTDIR}/${SAMPLE}/*.fa 2>/dev/null | wc -l)"
