#!/bin/bash
#SBATCH -J KneadData
#SBATCH -p cpunode
#SBATCH -o kneaddata_control.out
#SBATCH -e kneaddata_control.err
#SBATCH -N 1 -n 12

export LC_ALL="en_US.UTF-8"

ml load kneaddata-0.12.0
ml load parallel-20210922-gcc-8.5.0-v7ilxoy

KNEADDIR=/data/bio/zengyg/AdapterRemoval/Controls
BOWTIE2PATH=/opt/apps/spack/opt/spack/linux-centos7-skylake_avx512/gcc-8.5.0/bowtie2-2.4.2-kgezdpzdg2eaevx5wstksrjpz5bhhajh/bin/
TRIMMOMATIC=/opt/apps/spack/opt/spack/linux-centos7-skylake_avx512/gcc-8.5.0/trimmomatic-0.39-b262mrifwct67kpk3ss7m7jmfqlxnjhz/bin/

cat fastq12.txt | parallel -j 12 \
  kneaddata --bypass-trf \
    -i1 "$KNEADDIR/{}_adapterremove_1.fastq" \
    -i2 "$KNEADDIR/{}_adapterremove_2.fastq" \
    -db /data/bio/zengyg/Human_Genome_hg37 \
    -o /data/bio/zengyg/kneaddata/controls/{} \
    --remove-intermediate-output \
    --bowtie2 "$BOWTIE2PATH" \
    --trimmomatic "$TRIMMOMATIC"

