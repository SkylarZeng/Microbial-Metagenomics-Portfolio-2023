#!/bin/bash
#SBATCH -J Coverm
#SBATCH -p cpunode
#SBATCH -o coverm-%j.out
#SBATCH -e coverm-%j.err
#SBATCH -n 2

export LC_ALL="en_US.UTF-8"

ml load coverm-musl-0.6.1
ml load samtools-1.14-gcc-8.5.0-mxflewt
ml load minimap2-2.24-r1122
ml load bwa-0.7.17-gcc-8.5.0-oo4toav
ml load parallel-20210922-gcc-8.5.0-v7ilxoy

cat initial_samples.txt | parallel -j 1 \
  coverm genome \
    --coupled \
    /data/bio/zengyg/Coverm_Samples/read1/{}_repaired_paired_1.fastq.gz \
    /data/bio/zengyg/Coverm_Samples/read2/{}_repaired_paired_2.fastq.gz \
    --genome-fasta-files /data/bio/zengyg/Drep/drep_g_agrococcus/dereplicated_genomes/SRR12893903_bin5.fa \
    -o output.tsv

