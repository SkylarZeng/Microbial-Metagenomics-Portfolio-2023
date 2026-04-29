#!/bin/bash
#SBATCH --job-name=megahit
#SBATCH --output=logs/megahit_%A_%a.out
#SBATCH --error=logs/megahit_%A_%a.err
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=standard

# ============================================================
# 06_megahit_assembly.sh
# De novo metagenomic assembly with MEGAHIT.
# Assembles each sample independently (non-co-assembly).
# ============================================================

# ---------- user-configurable paths ----------
CLEAN_DIR="/scratch/$USER/metagenomics/decontaminated"
OUTDIR="/scratch/$USER/metagenomics/assembly"
SAMPLE_LIST="/scratch/$USER/metagenomics/metadata/sample_list.txt"
THREADS=$SLURM_CPUS_PER_TASK
MEM_BYTES=120000000000   # 120 GB (leave headroom below #SBATCH --mem)
# ---------------------------------------------

module purge
module load megahit/1.2.9

mkdir -p "$OUTDIR"
mkdir -p logs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
ASSEMBLY_DIR="${OUTDIR}/${SAMPLE}"

echo "[$(date)] MEGAHIT assembly starting for sample: $SAMPLE"

# Remove previous partial run (MEGAHIT refuses to overwrite)
rm -rf "$ASSEMBLY_DIR"

megahit \
    -1 "${CLEAN_DIR}/${SAMPLE}_R1.fastq.gz" \
    -2 "${CLEAN_DIR}/${SAMPLE}_R2.fastq.gz" \
    --out-dir       "$ASSEMBLY_DIR" \
    --out-prefix    "$SAMPLE" \
    --num-cpu-threads "$THREADS" \
    --memory        "$MEM_BYTES" \
    --k-list        21,29,39,59,79,99,119,141 \
    --min-contig-len 1000 \
    --verbose

# Summarise assembly statistics
echo "[$(date)] Assembly finished for sample: $SAMPLE"
echo "  Contig count: $(grep -c '>' ${ASSEMBLY_DIR}/${SAMPLE}.contigs.fa)"
echo "  Largest contig:"
awk '/^>/{if(seq) print length(seq); seq=""} !/^>/{seq=seq$0} END{if(seq) print length(seq)}' \
    "${ASSEMBLY_DIR}/${SAMPLE}.contigs.fa" | sort -n | tail -1
