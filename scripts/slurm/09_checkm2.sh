#!/bin/bash
#SBATCH --job-name=checkm2
#SBATCH --output=logs/checkm2_%A_%a.out
#SBATCH --error=logs/checkm2_%A_%a.err
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --dependency=afterok:METABAT2_JOBID

# ============================================================
# 09_checkm2.sh
# Assess completeness and contamination of MAGs with CheckM2.
# Uses a neural-network approach with the DIAMOND database.
# ============================================================

# ---------- user-configurable paths ----------
BINS_DIR="/scratch/$USER/metagenomics/bins/metabat2"
OUTDIR="/scratch/$USER/metagenomics/checkm2"
CHECKM2_DB="/scratch/$USER/databases/checkm2_db/uniref100.KO.1.dmnd"
SAMPLE_LIST="/scratch/$USER/metagenomics/metadata/sample_list.txt"
THREADS=$SLURM_CPUS_PER_TASK
# ---------------------------------------------

module purge
module load checkm2/1.0.1

mkdir -p "$OUTDIR"
mkdir -p logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
SAMPLE_BINS="${BINS_DIR}/${SAMPLE}"

echo "[$(date)] CheckM2 starting for sample: $SAMPLE"
echo "  Input bins directory: $SAMPLE_BINS"
echo "  Bin count: $(ls ${SAMPLE_BINS}/*.fa 2>/dev/null | wc -l)"

checkm2 predict \
    --input    "$SAMPLE_BINS" \
    --output-directory "${OUTDIR}/${SAMPLE}" \
    --extension fa \
    --database-path "$CHECKM2_DB" \
    --threads   "$THREADS" \
    --force

echo "[$(date)] CheckM2 finished for sample: $SAMPLE"

# Quick quality summary
echo "  High-quality bins (≥90% completeness, ≤5% contamination):"
awk -F'\t' 'NR>1 && $2>=90 && $3<=5 {print $1, $2, $3}' \
    "${OUTDIR}/${SAMPLE}/quality_report.tsv" | wc -l

echo "  Medium-quality bins (≥70% completeness, ≤10% contamination):"
awk -F'\t' 'NR>1 && $2>=70 && $3<=10 {print $1, $2, $3}' \
    "${OUTDIR}/${SAMPLE}/quality_report.tsv" | wc -l
