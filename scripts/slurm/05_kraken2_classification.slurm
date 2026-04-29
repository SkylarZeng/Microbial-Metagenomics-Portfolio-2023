#!/bin/bash
#SBATCH -J Kraken2
#SBATCH -p cpunode
#SBATCH -o kraken2-%j.out
#SBATCH -e kraken2-%j.err
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yige.zeng19@student.xjtlu.edu.cn

export LC_ALL="en_US.UTF-8"

ml load kraken2-2.1.1-gcc-8.5.0-5vbredn
ml load parallel-20210922-gcc-8.5.0-v7ilxoy

cat test1.txt | parallel -j 1 \
  kraken2 \
    --db /data/bio/zengyg/minikraken2_v1_8GB \
    --paired \
    /data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq \
    /data/bio/zengyg/Sample_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq \
    --use-names \
    --report /data/bio/zengyg/kraken2_report/{}_report2.txt \
    --report-zero-counts \
    --output /data/bio/zengyg/kraken2_output/kraken2.out

