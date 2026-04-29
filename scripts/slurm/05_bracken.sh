#!/bin/bash
#SBATCH --job-name=bracken
#SBATCH --output=logs/bracken_%A_%a.out
#SBATCH --error=logs/bracken_%A_%a.err
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=standard
#SBATCH --dependency=afterok:KRAKEN2_JOBID   # replace with actual job ID

# ============================================================
# 05_bracken.sh
# Bayesian abundance re-estimation from Kraken2 reports.
# Run separately at genus and species level.
# ============================================================

# ---------- user-configurable paths ----------
KRAKEN_DIR="/scratch/$USER/metagenomics/kraken2"
OUTDIR="/scratch/$USER/metagenomics/bracken"
DB="/scratch/$USER/databases/kraken2_standard_20230314"
SAMPLE_LIST="/scratch/$USER/metagenomics/metadata/sample_list.txt"
READ_LEN=150      # sequencing read length
# ---------------------------------------------

module purge
module load bracken/2.8

mkdir -p "${OUTDIR}/species"
mkdir -p "${OUTDIR}/genus"
mkdir -p logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

echo "[$(date)] Bracken starting for sample: $SAMPLE"

# Species-level
bracken \
    -d  "$DB" \
    -i  "${KRAKEN_DIR}/${SAMPLE}.k2report" \
    -o  "${OUTDIR}/species/${SAMPLE}_bracken_species.txt" \
    -w  "${OUTDIR}/species/${SAMPLE}_bracken_species.report" \
    -r  "$READ_LEN" \
    -l  S \
    -t  10

# Genus-level
bracken \
    -d  "$DB" \
    -i  "${KRAKEN_DIR}/${SAMPLE}.k2report" \
    -o  "${OUTDIR}/genus/${SAMPLE}_bracken_genus.txt" \
    -w  "${OUTDIR}/genus/${SAMPLE}_bracken_genus.report" \
    -r  "$READ_LEN" \
    -l  G \
    -t  10

echo "[$(date)] Bracken finished for sample: $SAMPLE"
