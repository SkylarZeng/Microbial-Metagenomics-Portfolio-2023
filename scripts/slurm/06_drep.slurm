#!/bin/bash
#SBATCH -J dRep
#SBATCH -p cpunode
#SBATCH -o dReptestplus-%j.out
#SBATCH -e dReptestplus-%j.err
#SBATCH -n 2

export LC_ALL="en_US.UTF-8"
export CHECKM_DATA_PATH=/data/public_data/checkm

ml load parallel-20210922-gcc-8.5.0-v7ilxoy
ml load metabat-2.15-gcc-8.5.0-vo32min
ml load bwa-0.7.17-gcc-8.5.0-oo4toav
ml load samtools-1.14-gcc-8.5.0-uvqyufe
ml load perl-5.34.1-gcc-8.5.0-qvpxqn2
ml load perl-net-ssleay-1.85-gcc-8.5.0-xcuifrw
ml load fraggenescan-1.31-gcc-8.5.0-gm2572v
ml load hmmer-3.3.2-gcc-8.5.0-r4vrjwg
ml load prodigal-2.6.3-gcc-8.5.0-ixjx62u
ml load py-matplotlib-3.5.2-gcc-8.5.0-2b4naql
ml load py-seaborn-0.9.0-gcc-8.5.0-nvtwd25
ml load python-3.9.12-gcc-8.5.0-ews4x7x
ml load pplacer-1.1.alpha19-gcc-8.5.0-a6jmrxw
ml load mummer-3.23-gcc-8.5.0-i6dbyy5
ml load gcc-11.3.0-gcc-8.5.0-bzzlmoa

dRep dereplicate dreptestplus/ \
  -comp 50 \
  -g /data/bio/zengyg/sample_genome_assemblies_genome_fasta/*.fa

