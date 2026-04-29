#!/bin/bash
#SBATCH -J Binning
#SBATCH -p cpunode
#SBATCH -o binning-%j.out
#SBATCH -e binning-%j.err
#SBATCH -n 8

export LC_ALL="en_US.UTF-8"

ml load metawrap-1.3.2
ml load parallel-20210922-gcc-8.5.0-v7ilxoy
ml load metabat-2.15-gcc-8.5.0-vo32min
ml load bwa-0.7.17-gcc-8.5.0-oo4toav
ml load samtools-1.14-gcc-8.5.0-uvqyufe
ml load perl-5.34.1-gcc-8.5.0-qvpxqn2
ml load perl-net-ssleay-1.85-gcc-8.5.0-xcuifrw
ml load fraggenescan-1.31-gcc-8.5.0-gm2572v
ml load hmmer-3.3.2-gcc-8.5.0-r4vrjwg
ml load maxbin2-2.2.7
ml load concoct-1.1.0

cat samples.txt | parallel -j 8 \
  metawrap binning \
    -o /data/bio/zengyg/Initial_Binning/{} \
    -t 2 \
    -a /data/bio/zengyg/Sample_Assembly/{}_assembly/final_assembly.fasta \
    --metabat2 --maxbin2 --concoct \
    /data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq \
    /data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq

