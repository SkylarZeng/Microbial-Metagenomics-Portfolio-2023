#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --output=logs/gtdbtk_%j.out
#SBATCH --error=logs/gtdbtk_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --partition=highmem
#SBATCH --dependency=afterok:CHECKM2_JOBID

# ============================================================
# 10_gtdbtk.sh
# Phylogenomic taxonomy placement of all medium/high-quality
# MAGs using GTDB-Tk (GTDB r214).
# Runs on ALL passing bins at once (not array job).
# ============================================================

# ---------- user-configurable paths ----------
BINS_DIR="/scratch/$USER/metagenomics/bins/metabat2"
CHECKM_DIR="/scratch/$USER/metagenomics/checkm2"
OUTDIR="/scratch/$USER/metagenomics/gtdbtk"
GTDBTK_DATA="/scratch/$USER/databases/gtdb_r214"
PASSING_BINS_DIR="/scratch/$USER/metagenomics/bins/passing"
MIN_COMPLETENESS=70
MAX_CONTAMINATION=10
THREADS=$SLURM_CPUS_PER_TASK
# ---------------------------------------------

module purge
module load gtdbtk/2.3.2

export GTDBTK_DATA_PATH="$GTDBTK_DATA"

mkdir -p "$PASSING_BINS_DIR"
mkdir -p "$OUTDIR"
mkdir -p logs

echo "[$(date)] Collecting passing MAGs (completeness >= ${MIN_COMPLETENESS}%, contamination <= ${MAX_CONTAMINATION}%)"

TOTAL_PASSING=0
for SAMPLE_CHECKM in "${CHECKM_DIR}"/*/quality_report.tsv; do
    SAMPLE=$(basename "$(dirname "$SAMPLE_CHECKM")")
    while IFS=$'\t' read -r BIN COMPLETENESS CONTAMINATION _REST; do
        [[ "$BIN" == "Name" ]] && continue
        if (( $(echo "$COMPLETENESS >= $MIN_COMPLETENESS" | bc -l) )) && \
           (( $(echo "$CONTAMINATION <= $MAX_CONTAMINATION" | bc -l) )); then
            cp "${BINS_DIR}/${SAMPLE}/${BIN}.fa" \
               "${PASSING_BINS_DIR}/${SAMPLE}__${BIN}.fa"
            TOTAL_PASSING=$((TOTAL_PASSING + 1))
        fi
    done < "$SAMPLE_CHECKM"
done

echo "  Total passing MAGs collected: $TOTAL_PASSING"

echo "[$(date)] Running GTDB-Tk classify_wf"
gtdbtk classify_wf \
    --genome_dir  "$PASSING_BINS_DIR" \
    --out_dir     "$OUTDIR" \
    --extension   fa \
    --cpus        "$THREADS" \
    --skip_ani_screen

echo "[$(date)] GTDB-Tk finished"
echo "  Results: ${OUTDIR}/classify/gtdbtk.bac120.summary.tsv"
