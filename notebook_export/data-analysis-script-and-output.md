# Data Analysis Script and Output

## Literature Review

[Nature综述：宏基因组测序研究耐药基因的方法和资源_刘永鑫Adam的博客-CSDN博客](https://blog.csdn.net/woodcorpse/article/details/115744071)

[Sequencing-based methods and resources to study antimicrobial resistance - Nature Reviews Genetics](https://www.nature.com/articles/s41576-019-0108-4)

# Data Preprocessing

## Fastqc

### Samples

```jsx
[zengyg@lnxbio fastQC]$ cat fastqc.slurm 
#!/bin/bash
#SBATCH -J fastqc
#SBATCH -p cpunode
#SBATCH -o fastqc1.out
#SBATCH -e fastqc1.err
#SBATCH -N 1 -n 4

ml load fastqc-0.11.9-gcc-8.5.0-couiege

FASTQDIR=/data/ncbi_data/fastq

for fastq12 in SRR12893917 SRR12893916  SRR12893913  SRR12893912  SRR12893911  SRR12893910  SRR12893909 SRR12893908 SRR12893907 SRR12893907 SRR12893906  SRR12893905 SRR12893904  SRR12893903 SRR12893902 SRR12893901 SRR12893901 SRR12893900 SRR12893900 SRR12893899 SRR12893899 SRR12893898 SRR12893898 SRR12893897 SRR12893897 SRR12893896 SRR12893896 SRR12893895 SRR12893895 SRR12893894 SRR12893894 SRR12893884 SRR12893884 SRR12893883 SRR12893883 SRR12893881 SRR12893881 SRR12893880 SRR12893880 SRR12893879 SRR12893879 SRR12893859 SRR12893859 SRR12893854 SRR12893854 SRR12893852 SRR12893852 SRR12893849 SRR12893849 SRR12893845 SRR12893845 SRR12893844 SRR12893842 SRR12893842 SRR12893841 SRR12893841 SRR12893834 SRR12893834 SRR12893831 SRR12893831 SRR12893830 SRR12893830 SRR12893829 SRR12893829 SRR12893826 SRR12893826  SRR12893824  SRR12893823  SRR12893822 SRR12893821  SRR12893820  SRR12893818 SRR12893816 SRR12893815  SRR12893814 SRR12893812   SRR12893811 SRR12893810   SRR12893804  SRR12893793 SRR12893782  SRR12893771  SRR12893694  SRR12893972  SRR12893950

do

fastqc $FASTQDIR/${fastq12}_{1,2}.fastq  -o  /data/bio/zengyg/Sample_Fastqc/

done
```

### [linux 复制多个文件夹下的文件到一个文件夹下面](https://www.cnblogs.com/awinlei/archive/2013/02/05/2893292.html)

`for i in $(find ./ -name *.gif);do cp -vf $i ./images/;done`

### Controls

```jsx
[zengyg@lnxbio fastQC]$ cat fastqc_control.slurm 
#!/bin/bash
#SBATCH -J fastqc
#SBATCH -p cpunode
#SBATCH -o fastqc2.out
#SBATCH -e fastqc2.err
#SBATCH -N 1 -n 4

ml load fastqc-0.11.9-gcc-8.5.0-couiege

FASTQDIR=/data/ncbi_data/fastq

for fastq12 in SRR12893993 SRR12893992 SRR12893991 SRR12893990 SRR12893989 SRR12893988 SRR12893987 SRR12893985 SRR12893984 SRR12893983 SRR12893982 SRR12893981
do

fastqc $FASTQDIR/${fastq12}_{1,2}.fastq  -o  /data/bio/zengyg/Control_Fastqc/

done
```

### gunzip

```jsx
#!/bin/bash
#SBATCH -J gunzip
#SBATCH -p cpunode
#SBATCH -o gunzip-%j.out
#SBATCH -e gunzip-%j.err
#SBATCH -N 1 -n 4

ml load parallel-20210922-gcc-8.5.0-v7ilxoy

cat Samples.txt | parallel -j 4 gunzip /data/bio/jiangmj/X101SC22124970-Z01-J001/00.CleanData/{}/*.fq.gz

```

## MultiQC

### Sample

```jsx
[zengyg@lnxbio multiQC]$ cat multiqc.slurm 
#!/bin/bash
#SBATCH -J multiqc
#SBATCH -p cpunode
#SBATCH -o multiqc1.out
#SBATCH -e multiqc1.err
#SBATCH -N 1 -n 1

ml load py-multiqc-1.7-gcc-8.5.0-cpm2hiw

MULTIQCDIR=/data/bio/zengyg/Sample_Fastqc

multiqc $MULTIQCDIR/*_fastqc.zip
```

### Control

```jsx
#!/bin/bash
#SBATCH -J multiqc
#SBATCH -p cpunode
#SBATCH -o multiqc1.out
#SBATCH -e multiqc1.err
#SBATCH -N 1 -n 1

ml load py-multiqc-1.7-gcc-8.5.0-cpm2hiw

MULTIQCDIR=/data/bio/zengyg/Control_Fastqc

multiqc $MULTIQCDIR/*_fastqc.zip

[zengyg@lnxbio multiQC]$
```

## AdapterRemoval

### Samples

```jsx
[zengyg@lnxbio zengyg]$ cat adapterremoval1.slurm 
#!/bin/bash
#SBATCH -J AdapterRemoval
#SBATCH -p cpunode
#SBATCH -o adapterremoval_sample.out
#SBATCH -e adapterremoval_sample.err
#SBATCH -N 1 -n 4

ml load AdapterRemoval-2.3.3 

FASTQDIR=/data/ncbi_data/fastq

for fastq12 in SRR12893917 SRR12893916  SRR12893913  SRR12893912  SRR12893911  SRR12893910  SRR12893909 SRR12893908  SRR12893907 SRR12893907 SRR12893906  SRR12893905 SRR12893904  SRR12893903 SRR12893902 SRR12893901 SRR12893901 SRR12893900 SRR12893900 SRR12893899 SRR12893899 SRR12893898 SRR12893898 SRR12893897 SRR12893897 SRR12893896 SRR12893896 SRR12893895 SRR12893895 SRR12893894 SRR12893894 SRR12893884 SRR12893884 SRR12893883 SRR12893883 SRR12893881 SRR12893881 SRR12893880 SRR12893880 SRR12893879 SRR12893879 SRR12893859 SRR12893859 SRR12893854 SRR12893854 SRR12893852 SRR12893852 SRR12893849 SRR12893849 SRR12893845 SRR12893845 SRR12893844 SRR12893842 SRR12893842 SRR12893841 SRR12893841 SRR12893834 SRR12893834 SRR12893831 SRR12893831 SRR12893830 SRR12893830 SRR12893829 SRR12893829 SRR12893826 SRR12893826  SRR12893824  SRR12893823  SRR12893822 SRR12893821  SRR12893820  SRR12893818 SRR12893816 SRR12893815  SRR12893814 SRR12893812   SRR12893811 SRR12893810   SRR12893804  SRR12893793 SRR12893782  SRR12893771  SRR12893694  SRR12893972  SRR12893950

do

AdapterRemoval --file1 $FASTQDIR/${fastq12}_1.fastq --file2 $FASTQDIR/${fastq12}_2.fastq --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --output1 /data/bio/zengyg/Sample_AdapterRemoval/${fastq12}_adapterremove_1.fastq --output2 /data/bio/zengyg/Sample_AdapterRemoval/${fastq12}_adapterremove_2.fastq

done
```

### Controls

```jsx
[zengyg@lnxbio zengyg]$ cat adapterremoval2.slurm
#!/bin/bash
#SBATCH -J AdapterRemoval
#SBATCH -p cpunode
#SBATCH -o adapterremoval_control.out
#SBATCH -e adapterremoval_control.err
#SBATCH -N 1 -n 4

ml load AdapterRemoval-2.3.3 

FASTQDIR=/data/ncbi_data/fastq

for fastq12 in SRR12893993 SRR12893992 SRR12893991 SRR12893990 SRR12893989 SRR12893988 SRR12893987 SRR12893985 SRR12893984 SRR12893983 SRR12893982 SRR12893981

do

AdapterRemoval --file1 $FASTQDIR/${fastq12}_1.fastq --file2 $FASTQDIR/${fastq12}_2.fastq --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --output1 /data/bio/zengyg/Control_AdapterRemoval/${fastq12}_adapterremove_1.fastq --output2 /data/bio/zengyg/Control_AdapterRemoval/${fastq12}_adapterremove_2.fastq

done
```

## KneadData (https://github.com/biobakery/kneaddata)

### Controls

```jsx
[zengyg@lnxbio zengyg]$ cat KneadData.slurm 
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

cat fastq12.txt | parallel -j 12 kneaddata --bypass-trf  -i1 /data/bio/zengyg/AdapterRemoval/Controls/{}_adapterremove_1.fastq -i2 /data/bio/zengyg/AdapterRemoval/Controls/{}_adapterremove_2.fastq -db /data/bio/zengyg/Human_Genome_hg37 -o /data/bio/zengyg/kneaddata/controls/{} --remove-intermediate-output --bowtie2 $BOWTIE2PATH  --trimmomatic $TRIMMOMATIC
```

### Build Bowtie-2 Database

[Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer)

```jsx
ml load kneaddata-0.12.0
bowtie2-build final_assembly.fasta Controls
```

### Samples

```jsx
[zengyg@lnxbio tmp1sst1]$ cat kneaddata_sample_1103.slurm 
#!/bin/bash
#SBATCH -J KneadData
#SBATCH -p cpunode
#SBATCH -o kneaddata_sample_controls.out
#SBATCH -e kneaddata_sample_controls.err
#SBATCH --ntasks-per-node=8

export LC_ALL="en_US.UTF-8"

ml load kneaddata-0.12.0
ml load parallel-20210922-gcc-8.5.0-v7ilxoy

BOWTIE2PATH=/opt/apps/spack/opt/spack/linux-centos7-skylake_avx512/gcc-8.5.0/bowtie2-2.4.2-kgezdpzdg2eaevx5wstksrjpz5bhhajh/bin/
TRIMMOMATIC=/opt/apps/spack/opt/spack/linux-centos7-skylake_avx512/gcc-8.5.0/trimmomatic-0.39-b262mrifwct67kpk3ss7m7jmfqlxnjhz/bin/

cat samples.txt | parallel -j 8 kneaddata  --bypass-trf -i1 /data/bio/tmp1sst1/AdapterRemoval/Samples/{}_adapterremove_1.fastq -i2 /data/bio/tmp1sst1/AdapterRemoval/Samples/{}_adapterremove_2.fastq  -db /data/bio/tmp1sst1/Human_Genome_hg37/hg37dec_v0.1 -t 2 -p 4 --remove-intermediate-output  -o /data/bio/tmp1sst1/kneaddata1103/{}  --bowtie2 $BOWTIE2PATH --trimmomatic $TRIMMOMATIC

```

```jsx
[zengyg@lnxbio tmp1sst1]$ cat kneaddata_sample_1103.slurm 
#!/bin/bash
#SBATCH -J KneadData
#SBATCH -p cpunode
#SBATCH -o kneaddata_sample_controls.out
#SBATCH -e kneaddata_sample_controls.err
#SBATCH --ntasks-per-node=8

export LC_ALL="en_US.UTF-8"

ml load kneaddata-0.12.0
ml load parallel-20210922-gcc-8.5.0-v7ilxoy

BOWTIE2PATH=/opt/apps/spack/opt/spack/linux-centos7-skylake_avx512/gcc-8.5.0/bowtie2-2.4.2-kgezdpzdg2eaevx5wstksrjpz5bhhajh/bin/
TRIMMOMATIC=/opt/apps/spack/opt/spack/linux-centos7-skylake_avx512/gcc-8.5.0/trimmomatic-0.39-b262mrifwct67kpk3ss7m7jmfqlxnjhz/bin/

cat samples.txt | parallel -j 8 kneaddata  --bypass-trf -i1 /data/bio/tmp1sst1/kneaddata1103/{}/{}_adapterremove_1_kneaddata_paired_1.fastq -i2 /data/bio/tmp1sst1/kneaddata1103/{}/{}_adapterremove_1_kneaddata_paired_2.fastq -db /data/bio/tmp1sst1/Controls_reference/Controls -t 2 -p 4 --remove-intermediate-output -o /data/bio/tmp1sst1/kneaddata1103/{}  --bowtie2 $BOWTIE2PATH --trimmomatic $TRIMMOMATIC

```

### 合并文件到目录

```jsx
find ./ -name "*paired_clean_*.fastq*" | xargs -I "{}" mv "{}" ./clean_1_2/
```

## MetaWRAP（https://github.com/bxlab/metaWRAP)

### Control_co-assembly

```jsx
cat CLEAN_READS/ERR*_1.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/ERR*_2.fastq > CLEAN_READS/ALL_READS_2.fastq

[zengyg@lnxbio zengyg]$ cat co-assembly_control.slurm 
#!/bin/bash
#SBATCH -J metaWRAP
#SBATCH -p cpunode
#SBATCH -o metawrap1.out
#SBATCH -e metawrap1.err
#SBATCH -n 4

ml load metawrap-1.3.2
ml load megahit-1.1.4-gcc-8.5.0-7uwl3gd
ml load py-quast-4.6.3-gcc-8.5.0-bsftrpy 
ml load bwa-0.7.17-gcc-8.5.0-oo4toav

METAWRAPDIR=/data/bio/zengyg/co_assembly_data
 
metawrap assembly -1 $METAWRAPDIR/ALL_READS_1.fastq -2 $METAWRAPDIR/ALL_READS_2.fastq -t 4 -o /data/bio/zengyg/Co_assembly_control
```

```jsx
[zengyg@lnxbio metawarp_sample_slurm_log]$ cat SRR12893811-3694.slurm 
#!/bin/bash
#SBATCH -J metaWRAP
#SBATCH -p cpunode
#SBATCH -o metawrap_3811-3694.out
#SBATCH -e metawrap_3811-3694.err
#SBATCH -n 8

ml load metawrap-1.3.2
ml load megahit-1.1.4-gcc-8.5.0-7uwl3gd
ml load py-quast-4.6.3-gcc-8.5.0-k7uacp6 
ml load bwa-0.7.17-gcc-8.5.0-oo4toav

metawrap assembly -1 /data/bio/zengyg/Samples_Clean/SRR12893811_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq -2 /data/bio/zengyg/Samples_Clean/SRR12893811_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq -t 8 -o /data/bio/zengyg/SRR12893811_assembly
```

## Binning

[MetaWRAP分箱流程实战和结果解读_刘永鑫Adam的博客-CSDN博客](https://blog.csdn.net/woodcorpse/article/details/106066754)

### Fastq Repair

```jsx
ml load bbmap-39.01
ml load parallel-20210922-gcc-8.5.0-v7ilxoy
 
 
cat samples.txt | parallel -j 30 repair.sh in=/data/bio/zengyg/Samples_Clean/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq in2=/data/bio/zengyg/Samples_Clean/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq out=/data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq out2=/data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq
```

### Rename

```jsx
rename adapterremove_1_kneaddata_paired_1_kneaddata repaired *
```

### Binning

```jsx
[zengyg@lnxbio zengyg]$ cat binning.slurm 
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

cat samples.txt | parallel -j 8 metawrap binning -o /data/bio/zengyg/Initial_Binning/{} -t 2 -a /data/bio/zengyg/Sample_Assembly/{}_assembly/final_assembly.fasta --metabat2 --maxbin2 --concoct /data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq /data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq
```

### Bin-Refinement

```jsx
[zengyg@lnxbio zengyg]$ cat binning.slurm 
#!/bin/bash
#SBATCH -J Bin-refinement
#SBATCH -p cpunode
#SBATCH -o bin-refinement-%j.out
#SBATCH -e bin-refinement-%j.err
#SBATCH -n 2

export LC_ALL="en_US.UTF-8"
export CHECKM_DATA_PATH=/data/public_data/checkm

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
ml load py-biopython-1.73-gcc-4.8.5-qfbrlp5
ml load py-checkm-genome-1.0.13-gcc-8.5.0-sdtawxc
ml load prodigal-2.6.3-gcc-8.5.0-ixjx62u

cat samples-binning.txt | parallel -j 8 metawrap bin_refinement -o /data/bio/zengyg/Bin_Refinement_1SRR_plot_numpy/SRR12893880 -t 2 -A /data/bio/zengyg/Initial_Binning/SRR12893880_backup2/metabat2_bins -B /data/bio/zengyg/Initial_Binning/SRR12893880/concoct_bins -c 50 -x 10
```

# Data Analysis

## Kraken2

[Kraken2 tutorial](https://bioinformaticsworkbook.org/dataAnalysis/Metagenomics/Kraken.html)

[Kraken使用手册](https://www.jianshu.com/p/df3f923451f2)

### Download Database

```jsx
wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v1_8GB_201904.tgz
```

### Build Database

```jsx
#!/bin/bash
#SBATCH -J Kraken2
#SBATCH -p cpunode
#SBATCH -o kraken2-%j.out
#SBATCH -e kraken2-%j.err
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Yige.zeng19@student.xjtlu.edu.cn

export LC_ALL="en_US.UTF-8"
export dataRoot=/data/public_data/checkm

ml load kraken2-2.1.1-gcc-8.5.0-5vbredn

kraken2-build --standard -t 8 -db /data/bio/zengyg/minikraken2_v1_8GB --use-ftp
```

## Kraken2

### Download Database

```jsx
#下载数据库
wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v1_8GB_201904.tgz
#解压
tar zxvf minikraken2_v1_8GB_201904_UPDATE.tgz -C ./ #绝对路径为～/db/minikraken2_v1_8GB/ 
```

### Bracken数据库构建

[Kraken2+Bracken](https://www.jianshu.com/p/dd8182e861cb)

```jsx
# Bracken数据库构建
[zengyg@lnxbio zengyg]$ cat yigeplus.slurm.sh 
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

cat test1.txt | parallel -j 1 kraken2 --db /data/bio/zengyg/minikraken2_v1_8GB  --paired /data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq /data/bio/zengyg/Sample_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq  --use-names --report /data/bio/zengyg/kraken2_report/{}_report2.txt --report-zero-counts --output /data/bio/zengyg/kraken2_output/kraken2.out

# Bracken参数解释：
-d，Kraken数据库路径；
-t，线程数；
-k，k-mer长度，Kraken1默认为31，Kraken2默认为32；
-l，reads读长。
```

### Kraken2

```jsx
[zengyg@lnxbio zengyg]$ cat kraken2-test.slurm 
#!/bin/bash
#SBATCH -J Kraken2
#SBATCH -p cpunode
#SBATCH -o kraken2-%j.out
#SBATCH -e kraken2-%j.err
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Yige.zeng19@student.xjtlu.edu.cn

export LC_ALL="en_US.UTF-8"

ml load kraken2-2.1.1-gcc-8.5.0-5vbredn
ml load parallel-20210922-gcc-8.5.0-v7ilxoy

kraken2 --db /data/bio/zengyg/minikraken2_v1_8GB_Update --report /data/bio/zengyg/kraken2_report/SRR12893972_report --paired /data/bio/zengyg/Samples_Clean_repair/SRR12893972_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq /data/bio/zengyg/Sample_Clean_repair/SRR12893972_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq --use-names --report-zero-counts --output /data/bio/zengyg/kraken2_output
```

### Find software

```jsx
spack find | grep prokka

prokka@1.14.6
```

### bin重命名

```jsx
[zengyg@lnxbio zengyg]$ cat rename_mv_bins.sh 
#!/bin/bash

for SRR in SRR12893771 SRR12893782 SRR12893793 SRR12893804 SRR12893810 SRR12893811 SRR12893812 SRR12893814 SRR12893815 SRR12893816 SRR12893818 SRR12893820 SRR12893821 SRR12893822 SRR12893823 SRR12893824 SRR12893826 SRR12893829 SRR12893830 SRR12893831 SRR12893834 SRR12893841 SRR12893842 SRR12893844 SRR12893845 SRR12893849 SRR12893852 SRR12893854 SRR12893859 SRR12893881 SRR12893883 SRR12893884 SRR12893894 SRR12893895 SRR12893896 SRR12893897 SRR12893898 SRR12893899 SRR12893900 SRR12893901 SRR12893902 SRR12893903 SRR12893904 SRR12893905 SRR12893906 SRR12893907 SRR12893908 SRR12893909 SRR12893910 SRR12893911 SRR12893912 SRR12893913 SRR12893916 SRR12893917 SRR12893950 SRR12893972; do
   for BINSN in {1..9}; do
	if [[ -f "$SRR"/metawrap_50_10_bins/bin."$BINSN".fa ]];then
	mv  "$SRR"/metawrap_50_10_bins/bin."$BINSN".fa  "$SRR"/metawrap_50_10_bins/"$SRR"_bin"$BINSN"
	fi
   done

   find "$SRR/metawrap_50_10_bins/" -name "$SRR"_bin"?" | xargs -i mv {} ./SRRBIN/

done
```

### 批量加后缀

```jsx
ls | xargs -t -i mv {} {}.fa
```

## CheckM

[科学网-CheckM评估基因组完整度 - 涂波的博文](https://blog.sciencenet.cn/blog-2379401-1068993.html)

[拿到Bin以后，会用才是高级分析！](https://zhuanlan.zhihu.com/p/328454924)

```jsx
[zengyg@lnxbio zengyg]$ cat checkm.slurm 
#!/bin/bash
#SBATCH -J CheckM
#SBATCH -p cpunode
#SBATCH -o checkm-%j.out
#SBATCH -e checkm-%j.err
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yige.zeng19@student.xjtlu.edu.cn

export LC_ALL="en_US.UTF-8"

ml load python-3.9.12-gcc-8.5.0-ews4x7x
ml load py-checkm-genome-1.0.13-gcc-8.5.0-sdtawxc
ml load hmmer-3.3.2-gcc-8.5.0-r4vrjwg
ml load prodigal-2.6.3-gcc-8.5.0-ixjx62u
ml load pplacer-1.1.alpha19-gcc-8.5.0-a6jmrxw
ml load fasttree-2.1.10-gcc-8.5.0-jzumr6k
export GTDBTK_DATA_PATH=/data/public_data/gtdbtk-data/release207_v2

checkm lineage_wf -x fa SRRBIN/ checkm/
```

### TSV转Excel

[宏基因组分箱CheckM评估结果的提取_EmmettPeng的博客-CSDN博客_checkm使用](https://blog.csdn.net/Emmett_Bioinfo/article/details/116891391)

```jsx
[zengyg@lnxbio checkm_2_excel]$ cat readme 
# run

python3 checmsum.py bin_stats_ext.tsv bin_stats_ext.txt

# input bin_stats_ext.txt TO excel
```

```jsx
[zengyg@lnxbio checkm_2_excel]$ cat checkmsum.py 
#!/bin/python3

'''
This program is designed to process CheckM output file (e.g. bin_stats_ext.tsv) to easy-to-read text format (.txt)
Usage: python3 checkm_summary.py <inputfile> <outputfile>
Written by: Emmett Peng
'''
import json
import sys

with open(sys.argv[1] ,'r') as f:
	Load = {}
	for line in f:
		line = line.replace('\'','\"')
		line = line.split('\t')
		line[0] = 'bin_' + line[0].replace('.bin.', '_')
		#print(line[0])
		#line[0]:Bin Id; line[1]:Bin information
		#exec(f"line[0] = json.loads(line[1])")
		Load[line[0]] = json.loads(line[1])
		print(Load)
		#print(text)

with open(sys.argv[2], 'w+') as output:
	output.write('Bin Id\tMarker Lineage\tGenomes\tMarkers\tMarker Sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tGC\tGC std\tGenome Size\tAmbiguous Bases\tScaffolds\tContigs\tTranslation Table\tPredicted Genes\n')
	for key in Load:
		output.write(key + '\t')
		output.write(Load[key]['marker lineage'] + '\t')
		output.write(str(Load[key]['# genomes']) + '\t')
		output.write(str(Load[key]['# markers']) + '\t')
		output.write(str(Load[key]['# marker sets']) + '\t')
		output.write(str(Load[key]['0']) + '\t')
		output.write(str(Load[key]['1']) + '\t')
		output.write(str(Load[key]['2']) + '\t')
		output.write(str(Load[key]['3']) + '\t')
		output.write(str(Load[key]['4']) + '\t')
		output.write(str(Load[key]['5+']) + '\t')
		output.write(str(Load[key]['Completeness']) + '\t')
		output.write(str(Load[key]['Contamination']) + '\t')
		output.write(str(Load[key]['GC']) + '\t')
		output.write(str(Load[key]['GC std']) + '\t')
		output.write(str(Load[key]['Genome size']) + '\t')
		output.write(str(Load[key]['# ambiguous bases']) + '\t')
		output.write(str(Load[key]['# scaffolds']) + '\t')
		output.write(str(Load[key]['# contigs']) + '\t')
		output.write(str(Load[key]['Translation table']) + '\t')
		output.write(str(Load[key]['# predicted genes']) + '\t')
		output.write('\n')[zengyg@lnxbio checkm_2_excel]
```

## GTDB-tk

[GTDB：基因组分类数据库，物种注释和进化树构建工具GTDB-tk_刘永鑫Adam的博客-CSDN博客](https://blog.csdn.net/woodcorpse/article/details/108924563)

### 运行软件

```jsx
ml load python-3.9.12-gcc-8.5.0-ews4x7x   pplacer-1.1.alpha19-gcc-8.5.0-a6jmrxw  hmmer-3.3.2-gcc-8.5.0-r4vrjwg  fasttree-2.1.10-gcc-8.5.0-jzumr6k

export GTDBTK_DATA_PATH=/data/public_data/gtdbtk-data/release207_v2
```

### 检查数据库

```jsx
gtdbtk check_install
```

### 测试流程

```jsx
ml load python-3.9.12-gcc-8.5.0-ews4x7x   pplacer-1.1.alpha19-gcc-8.5.0-a6jmrxw  hmmer-3.3.2-gcc-8.5.0-r4vrjwg

ml load  fasttree-2.1.10-gcc-8.5.0-jzumr6k

gtdbtk test --out_dir gtdbtk_test
```

### GTDB-tk for 186genomes

```jsx
[zengyg@lnxbio zengyg]$ cat GTDB-tk.slurm 
#!/bin/bash
#SBATCH -J GTDB-tk
#SBATCH -p cpunode
#SBATCH -o gtdbtk-%j.out
#SBATCH -e gtdbtk-%j.err
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yige.zeng19@student.xjtlu.edu.cn

export LC_ALL="en_US.UTF-8"

ml load python-3.9.12-gcc-8.5.0-ews4x7x
ml load pplacer-1.1.alpha19-gcc-8.5.0-a6jmrxw
ml load hmmer-3.3.2-gcc-8.5.0-r4vrjwg
ml load fasttree-2.1.10-gcc-8.5.0-jzumr6k
export GTDBTK_DATA_PATH=/data/public_data/gtdbtk-data/release207_v2

gtdbtk classify_wf --genome_dir /data/bio/zengyg/SRRBIN \
    --extension fa \
    --out_dir /data/bio/zengyg/GTDB_classify_wf \
    --cpus 1
```

## R绘图

[基于R语言的微生物群落组成多样性分析--物种丰度计算及可视化](https://zhuanlan.zhihu.com/p/539464200#:~:text=%E9%97%A8%E6%B0%B4%E5%B9%B3%E7%9A%84%E7%89%A9%E7%A7%8D%E4%B8%B0%E5%BA%A6%E8%AE%A1%E7%AE%97%E5%8F%8A%E5%8F%AF%E8%A7%86%E5%8C%96%201%201%E3%80%81%E8%AF%BB%E5%8F%96%E6%95%B0%E6%8D%AE%20df1%20%3C-%20read.table%28file%3D%22Phylum.txt%22%2Csep%3D%22%5Ct%22%2Cheader%3DT%2Ccheck.names%3DFALSE%29%202%202%E3%80%81%E5%AF%B9%E6%95%B0%E6%8D%AE%E8%BF%9B%E8%A1%8C%E5%A4%84%E7%90%86,rownames%28df2%29%3Ddf2%24Tax%20df3%3Ddf2%5B%2C-1%5D%203%203%E3%80%81%E8%AE%A1%E7%AE%97%E7%89%A9%E7%A7%8D%E6%80%BB%E4%B8%B0%E5%BA%A6%E5%B9%B6%E8%BF%9B%E8%A1%8C%E9%99%8D%E5%BA%8F%E6%8E%92%E5%88%97%20...%204%204%E3%80%81%E7%BB%98%E5%9B%BE%20)

[跟着Nat Commun学作图 | 3.物种丰度堆积柱状图](https://zhuanlan.zhihu.com/p/433683673)

[R语言-使用ggplot2绘制物种丰度堆叠柱状图](https://www.jianshu.com/p/594ddf604159)

[Diversity statistics](https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html)

[Data Visualization with R](https://rkabacoff.github.io/datavis/Bivariate.html)

[Grouped, stacked and percent stacked barplot in ggplot2](https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html)

```jsx
install.packages("devtools")
install.packages("stringr")
library(devtools)
install_github('wilkox/wilkoxmisc')
install.packages("starsExtra")

library(wilkoxmisc)
library(reshape2)
library(dplyr)
library(readr)
library(ggplot2)
library(starsExtra)

OTU <- read.csv("kraken_report_all_species.csv")

OTU <- OTU %>%
gather(Sample, Count, 2:(ncol(OTU))) %>%
filter(Count > 0)

names(OTU)[1] <- c("Species")
OTU$Sample <- gsub("kraken.report.", "",OTU$Sample)
OTU$Sample <- gsub("_kraken_report", "",OTU$Sample)

Meta <- read_tsv("metadata_doorknob_skin.txt")
OTU <- OTU %>% filter(OTU$Sample %in% Meta$Sample)
#Tabulate read counts by genus
OTU <- ddply(OTU, .(Sample, Species), summarise, count = sum(Count))

#Convert count to relativeabundance and add column (require dplyr)
OTU <- ddply(OTU, .(Sample), mutate, RelativeAbundance = (count * 100) / sum(count))

write_tsv(OTU,"kraken_species_abundance.tidy.txt")

#Open taxonomy OTU table
OTU <- read_tsv("kraken_species_abundance.tidy.txt")
#Merge relative abundance table and metatable together with filter condition
Meta <- read_tsv("metadata_doorknob_skin.txt")
Meta2 <- Meta %>% filter(Location == "Residence 3")
OTUTable <- merge(OTU, Meta2, by = "Sample", all.x = TRUE)
OTU3 <- OTUTable %>% filter(Location == "Residence 3")

#collapse taxa table to only 5 or 8 top phyla, genus, family, etc (require reshape2).
Top12 <- collapse_taxon_table(OTU3, n = 12, Rank = "Species")

#Merge relative abundance table and metatable together 
Top12species <- merge(Top12, Meta2, by = "Sample", all.x = TRUE)
write_tsv(Top12species, "Top12species.txt")

#Plot
SpeciesAbundance <- read_tsv("Top12Species.txt")
Plot <- ggplot(SpeciesAbundance, aes(x = factor(Sampling_day,levels = c("Day 1 day","Day 1 night","Day 2 day","Day 2 night","Day 3 day","Day 3 night","Day 4 day","Day 4 night","Day 5 day","Day 5 night","Day 6 day","Day 6 night","Day 7 day","Day 7 night","Day 8 day","Day 8 night","Day 9 day","Day 9 night","Day 10 day","Day 10 night")), y = RelativeAbundance, fill = factor(Species, levels = c("Chryseobacterium taklimakanense", "Cutibacterium acnes", "Dietzia lutea","Dietzia sp. oral taxon 368", "Gordonia bronchialis", "Gordonia sp. KTR9", "Gordonia terrae","Janibacter indicus", "Kytococcus sedentarius", "Micrococcus luteus", "Moraxella osloensis",   "Minor/Unclassified"))))
Plot <- Plot + geom_bar(stat="identity")
Plot <- Plot + scale_fill_brewer(palette = "Paired")
Plot <- Plot + theme_classic()
Plot <- Plot + facet_grid(Sample_type~Location,scales = "free",space="free") +scale_y_continuous(expand = c(0,0))
Plot <- Plot + ylab(paste0("Relative Abundance (%)"))+xlab(paste0("Sampling day"))+labs(fill = "Species")
Plot <- Plot + theme(axis.text.x = element_text(hjust=0,angle=90))

#Export plot
ggsave("taxonomy_general.png", width=12, height=7)
ggsave("top12species_by_sample_type.pdf",dpi=1086)
```

### Results

![taxonomy_general.png](Data%20Analysis%20Script%20and%20Output/taxonomy_general.png)

![taxonomy_general.png](Data%20Analysis%20Script%20and%20Output/taxonomy_general%201.png)

## Bray-Curtis dissimilarity

### Density Plot

```jsx
install.packages("cultevo")

library(reshape2)
library(wilkoxmisc)
library(ggplot2)
library(dplyr)
library(cultevo)

tax_cast <- read.table("kraken_species_abundance.cast.txt",header=T,row.names=1)
#distance <- distance(tax_cast, method = "bray-curtis") in R package "ecodist
library(vegan)
distance <- vegdist(tax_cast, method="bray")
Matrix <- as.matrix(distance)
library(MASS)
write.matrix(Matrix, sep=" ", "bray_curtis_distance_matrix.txt")

#Load weighted UniFrac distances
UniFrac <- read_tsv("bray_curtis_distance_matrix.txt") # judge table in excel

#Melt
UniFrac <- melt(UniFrac, value.name = "Distance")
names(UniFrac)[1:2] <- c("Sample1", "Sample2")

#Add household and type
Samples <- metadata_doorknob_skin %>% select(Sample, Sample_type, Location, Sampling_day)
#Samples <- metadata_doorknob_skin[,c(1,5,10,11)] #Alternative 
names(Samples)[1] <- c("Sample2")
UniFrac <- merge(UniFrac, Samples, by = "Sample2", all.x = FALSE)
names(UniFrac)[4:5] <- c("Sample_type2", "Location2")
names(UniFrac)[6] <- c("Sampling_day2")
names(Samples)[1] <- c("Sample1")
UniFrac <- merge(UniFrac, Samples, by = "Sample1", all.x = FALSE)
names(UniFrac)[7:8] <- c("Sample_type1", "Location1")
names(UniFrac)[9] <- c("Sampling_day1")

time <- read_tsv("sampling_time_hour_conversion.txt")
names(time)[1] <- c("Sampling_day1")
merge <- left_join(UniFrac, time)
names(merge)[10] <- c("Time_hour1")
names(time)[1] <- c("Sampling_day2")
merge <- left_join(merge, time)
names(merge)[11] <- c("Time_hour2")

UniFrac <- merge
#Remove self-self samples
UniFrac <- UniFrac[which( ! UniFrac$Sample1 == UniFrac$Sample2), ]

UniFrac$HouseholdType <- ifelse(UniFrac$Location1 == UniFrac$Location2, "Same household", "Different household")

table <- UniFrac %>%
    filter(HouseholdType == "Same household") %>%
    filter(Sample_type2 == "Door knob") %>%
    filter(Sample_type1 %in% c("Right palm", "Left palm"))

table$Group <- paste0(table$Sample_type2," vs. ", table$Sample_type1)
table <- table %>%
    mutate(time_decay =table$Time_hour1-table$Time_hour2)
subtable1 <- table %>% filter(Location1 == "Residence 3")

devtools::install_github('Mikata-Project/ggthemr')
library(ggthemr)
Plot1 <- ggplot(subtable, aes(x=Location1, y = Distance, fill = Group))
#Plot1 <- Plot1 + geom_boxplot(width=0.5)
Plot1 <- Plot1 + theme_classic()
#Plot1 <- Plot1 + facet_wrap(~Location1)
#Plot1 <- Plot1 + geom_boxplot(width=0.5,outlier.shape = NA)
Plot1 <- Plot1 + geom_boxplot()
#Plot1 <- Plot1 + geom_jitter(size=1, alpha=0.3, position=position_jitter(0.2),aes(color=Location1))
#Plot1 <- Plot1 + scale_color_brewer(palette = "Set1")
Plot1 <- Plot1 + scale_fill_brewer(palette = "Paired")
#Plot1 <- Plot1 + facet_wrap(~Species,scales = "free_y")+ scale_y_continuous(expand = c(0,0))
Plot1 <- Plot1 + theme(legend.position="none")
Plot1 <- Plot1 + theme(axis.title = element_text(size=12))
Plot1 <- Plot1 + theme(axis.text.x = element_text(size = 10, hjust = 1, angle=45))
Plot1 <- Plot1 + ylab(paste0("Bray-Curtis dissimilarity")) + xlab(paste0("Residence"))
Plot1 <- Plot1 + theme(axis.title.x= element_blank())
Plot1 <- Plot1 + theme(legend.position = "none")
#Plot1 <- Plot1 + ggtitle("Deinococcus-Thermus") + theme(plot.title = element_text(hjust = 0.5))
Plot1

subtable <- table %>% filter(Location1 == "Residence 3")

wilcox.test(Distance~Group, data=subtable)
# residence 1
Wilcoxon rank sum test with continuity correction
data:  Distance by Group
W = 46123, p-value = 0.0003799
alternative hypothesis: true location shift is not equal to 0

# residence 2
Wilcoxon rank sum test with continuity correction
data:  Distance by Group
W = 40549, p-value = 0.005159
alternative hypothesis: true location shift is not equal to 0

# residence 3
Wilcoxon rank sum test with continuity correction
data:  Distance by Group
W = 72650, p-value = 0.2869
alternative hypothesis: true location shift is not equal to 0

# residence 4
Wilcoxon rank sum test with continuity correction
data:  Distance by Group
W = 13201, p-value = 0.07187
alternative hypothesis: true location shift is not equal to 0

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ddply")
library(ddply)
mu <- ddply(table, "Group", summarise, grp.mean=mean(Distance))
Group  grp.mean
1 Door knob vs. Left palm 0.4774842
2 Door knob vs. Right palm 0.4801904

kruskal.test(Distance~HouseholdComparison,data=UniFrac)
kruskalmc(Distance~HouseholdComparison,data=UniFrac)

cbPalette <- c("#fd017d", "#555def")
UniFrac$Comparisonf <- factor(UniFrac$Comparisonf, c("Within Individuals", "Between Individuals Within Households", "Between Households"))

P1 <- ggplot(subtable, aes(x = Distance, color = Group))
P1 <- P1 + geom_line(stat="density", size =2) + scale_y_continuous(expand=c(0,0), limits = c(0,3))
P1 <- P1 + theme_classic()
P1 <- P1 + facet_wrap(~Location1)
P1 <- P1 + theme(panel.background = element_rect(colour = "Black"))
P1 <- P1 + scale_colour_manual(values = cbPalette)
P1 <- P1 + xlab(paste0("Normalized Bray Curtis Dissimilarity")) + ylab(paste0("Density (%)"))
P1 <- P1 + theme(legend.title=element_blank())
P1 <- P1 + geom_segment(aes(x=0.480,y=0,xend=0.480,yend=1.45),size=1,linetype="dotdash",color="#fd017d")
P1 <- P1 + geom_segment(aes(x=0.477,y=0,xend=0.477,yend=1.2),size=1,linetype="dotdash",color="#455def")
P1 <- P1 + annotate("text", x=0.6, y=0.9, label="Mean = 0.480",angle=0, color="#fd017d")
P1 <- P1 + annotate("text", x=0.35, y=0.6, label="Mean = 0.477",angle=0, color="#455def")
P1 <- P1 + annotate("text", x=0.6, y=2.5, label = "P = 0.2869", size = 8,  fontface=3)
P1 <- P1 + theme(legend.position="bottom")
P1
ggsave(plot = P1, "doorkob_palm_comparison_Bray_Curtis_residence3.png", width=7, height=5)
write_tsv(table, "doorkob_palm_comparison_Bray_Curtis.txt")
write_tsv(subtable, "doorkob_palm_comparison_Bray_Curtis_residence3.txt")

UniFrac$AreaType <- ifelse(UniFrac$Area1 == UniFrac$Area2, "Same Area", "Different Area")
UniFrac$IndividualType <- ifelse(UniFrac$Individual1 == UniFrac$Individual2, "Same Individual","Different Individual")
UniFrac$AgeGroup <- ifelse(UniFrac$Age_Group1 == UniFrac$Age_Group2, "Same Age Group","Different Age Group")
UniFrac$GenderType <- ifelse(UniFrac$Gender1 == UniFrac$Gender2, "Same Gender","Different Gender")
UniFrac$AnatomyType <- ifelse(UniFrac$Anatomy1 == UniFrac$Anatomy2, "Same Anatomy","Different Anatomy")

UniFrac$Comparisonf <- ifelse(UniFrac$Individual1 == UniFrac$Individual2, "Within Individuals", "Between Individuals")
UniFrac$Comparisonf <- ifelse((UniFrac$Comparisonf == "Between Individuals") & (UniFrac$Household1 == UniFrac$Household2), "Between Individuals Within Households", UniFrac$Comparisonf)
UniFrac$Comparisonf <- ifelse(UniFrac$Comparisonf == "Between Individuals", "Between Households",UniFrac$Comparisonf)
```

### Line plot

```jsx
library(reshape2)
library(wilkoxmisc)
library(ggplot2)
library(dplyr)
library(cultevo)

tax_cast <- read.table("kraken_species_abundance.cast.txt",header=T,row.names=1)
#distance <- distance(tax_cast, method = "bray-curtis") in R package "ecodist
library(vegan)
distance <- vegdist(tax_cast, method="bray")
Matrix <- as.matrix(distance)
library(MASS)
write.matrix(Matrix, sep=" ", "bray_curtis_distance_matrix.txt")

#Load weighted UniFrac distances
UniFrac <- read.tsv("bray_curtis_distance_matrix.txt")

#Melt
UniFrac <- melt(UniFrac, value.name = "Distance")
names(UniFrac)[1:2] <- c("Sample1", "Sample2")

#Add household and type
Samples <- read_tsv("metadata_doorknob_skin.txt") %>%
    select(Sample, Sample_type, Location, Sampling_day)
names(Samples)[1] <- c("Sample2")
UniFrac <- merge(UniFrac, Samples, by = "Sample2", all.x = FALSE)
names(UniFrac)[4:5] <- c("Sample_type2", "Location2")
names(UniFrac)[6] <- c("Sampling_day2")
names(Samples)[1] <- c("Sample1")
UniFrac <- merge(UniFrac, Samples, by = "Sample1", all.x = FALSE)
names(UniFrac)[7:8] <- c("Sample_type1", "Location1")
names(UniFrac)[9] <- c("Sampling_day1")

time <- read_tsv("sampling_time_hour_conversion.txt")
names(time)[1] <- c("Sampling_day1")
merge <- left_join(UniFrac, time)
names(merge)[10] <- c("Time_hour1")
names(time)[1] <- c("Sampling_day2")
merge <- left_join(merge, time)
names(merge)[11] <- c("Time_hour2")

UniFrac <- merge
#Remove self-self samples
UniFrac <- UniFrac[which( ! UniFrac$Sample1 == UniFrac$Sample2), ]

UniFrac$HouseholdType <- ifelse(UniFrac$Location1 == UniFrac$Location2, "Same household", "Different household")
UniFrac <- UniFrac %>% filter(Location1 == "Residence 3")
UniFrac <- UniFrac %>% filter(Location2 == "Residence 3")
table <- UniFrac %>%
    filter(HouseholdType == "Same household") %>%
    filter(Sample_type2 == "Door knob") %>%
    filter(Sample_type1 %in% c("Right palm", "Left palm"))

table$Group <- paste0(table$Sample_type2," vs. ", table$Sample_type1)
table <- table %>%
    mutate(time_decay =table$Time_hour1-table$Time_hour2)

table2 <- table %>%
    filter(time_decay >= 0)

#cbPalette <- c("#fd017d", "#255def")
Plot1 <- ggplot(table, aes(x=time_decay, y = Distance))
Plot1 <- Plot1 + geom_point(color="#455def",alpha = 0.4)
Plot1 <- Plot1 + geom_smooth(color="#fd017d", fill="#fd017d")
Plot1 <- Plot1 + theme_bw()
#Plot1 <- Plot1 + scale_color_manual(values = cbPalette)
#Plot1 <- Plot1 + facet_grid(Group~Location1)
Plot1 <- Plot1 + geom_vline(xintercept=0, linetype="dotted", color = "black", size=0.5)
#Plot1 <- Plot1 + geom_boxplot(width=0.5,outlier.shape = NA)
#Plot1 <- Plot1 + geom_boxplot()
#Plot1 <- Plot1 + geom_jitter(size=1, alpha=0.3, position=position_jitter(0.2),aes(color=Location1))
#Plot1 <- Plot1 + scale_color_brewer(palette = "Set1")
#Plot1 <- Plot1 + scale_color_brewer(palette = "Set2")
#Plot1 <- Plot1 + facet_wrap(~Species,scales = "free_y")+ scale_y_continuous(expand = c(0,0))
#Plot1 <- Plot1 + theme(legend.position="none")
#Plot1 <- Plot1 + theme(axis.title = element_text(size=12))
#Plot1 <- Plot1 + theme(axis.text.x = element_text(size = 10, hjust = 1, angle=45))
Plot1 <- Plot1 + ylab(paste0("Bray-Curtis dissimilarity between doorkonb and palm")) + xlab(paste0("Time interval (hours)"))
#Plot1 <- Plot1 + theme(axis.title.x= element_blank())
#Plot1 <- Plot1 + ggtitle("Deinococcus-Thermus") + theme(plot.title = element_text(hjust = 0.5))
Plot1 <- Plot1 + theme(legend.position = "none")
Plot1
ggsave(plot = Plot1, "Bray-Curtis dissimilarity between doorknob and palm.png", width = 5, height = 5)

#library(ggpubr)
#library(gridExtra)
#Figure <- grid.arrange(Plot1, P1, ncol=2)
library(cowplot)
Figure <- plot_grid(Plot1, P1, labels = LETTERS[1:2], ncol = 2)
ggsave(plot = Figure, "Normalized_Bray_Curtis_Dissimilarity_combined_by_sample_type.png", width = 10, height = 5)

subtable <- table %>% filter(Location1 == "Residence 3")

wilcox.test(Distance~Group, data=table)

# residence 1
Wilcoxon rank sum test with continuity correction
data:  Distance by Group
W = 46123, p-value = 0.0003799
alternative hypothesis: true location shift is not equal to 0

# residence 2
Wilcoxon rank sum test with continuity correction
data:  Distance by Group
W = 40549, p-value = 0.005159
alternative hypothesis: true location shift is not equal to 0

# residence 3
Wilcoxon rank sum test with continuity correction
data:  Distance by Group
W = 72650, p-value = 0.2869
alternative hypothesis: true location shift is not equal to 0

# residence 4
Wilcoxon rank sum test with continuity correction
data:  Distance by Group
W = 13201, p-value = 0.07187
alternative hypothesis: true location shift is not equal to 0

library(ddply)
mu <- ddply(table, "Group", summarise, grp.mean=mean(Distance))
Group  grp.mean
1 Door knob vs. Left palm 0.4774842
2 Door knob vs. Right palm 0.4801904

kruskal.test(Distance~HouseholdComparison,data=UniFrac)
kruskalmc(Distance~HouseholdComparison,data=UniFrac)

#cbPalette <- c("#590A30", "#90AA3C", "#EF6125")
#UniFrac$Comparisonf <- factor(UniFrac$Comparisonf, c("Within Individuals", "Between Individuals Within Households", "Between Households"))

P1 <- ggplot(table, aes(x = Distance, color = Group))
P1 <- P1 + geom_line(stat="density", size =2) + scale_y_continuous(expand=c(0,0), limits = c(0,5))
P1 <- P1 + theme_classic()
P1 <- P1 + facet_wrap(~Location1)
P1 <- P1 + theme(panel.background = element_rect(colour = "Black"))
P1 <- P1 + scale_colour_manual(values = cbPalette)
P1 <- P1 + xlab(paste0("Normalized Bray Curtis Dissimilarity")) + ylab(paste0("Density (%)"))
P1 <- P1 + theme(legend.title=element_blank())
P1 <- P1 + geom_segment(aes(x=0.417,y=0,xend=0.417,yend=1.4),size=1,linetype="dotdash",color="#590A30")
P1 <- P1 + geom_segment(aes(x=0.484,y=0,xend=0.484,yend=0.8),size=1,linetype="dotdash",color="#90AA3C")
P1 <- P1 + geom_segment(aes(x=0.499,y=0,xend=0.499,yend=0.95),size=1,linetype="dotdash",color="#EF6125")
P1 <- P1 + annotate("text", x=0.355, y=1.25, label="Mean = 0.417",angle=315, color="#590A30")
P1 <- P1 + annotate("text", x=0.355, y=0.7, label="Mean = 0.484",angle=315, color="#90AA3C")
P1 <- P1 + annotate("text", x=0.52, y=1.2, label="Mean = 0.499",angle=315, color="#EF6125")
#P1 <- P1 + annotate("text", x=0.5, y=3.5, label = "P < 0.001", size = 8,  fontface=3)
P1 <- P1 + theme(legend.position="bottom")

write_tsv(table, "doorkob_palm_comparison_Bray_Curtis.txt")

UniFrac$AreaType <- ifelse(UniFrac$Area1 == UniFrac$Area2, "Same Area", "Different Area")
UniFrac$IndividualType <- ifelse(UniFrac$Individual1 == UniFrac$Individual2, "Same Individual","Different Individual")
UniFrac$AgeGroup <- ifelse(UniFrac$Age_Group1 == UniFrac$Age_Group2, "Same Age Group","Different Age Group")
UniFrac$GenderType <- ifelse(UniFrac$Gender1 == UniFrac$Gender2, "Same Gender","Different Gender")
UniFrac$AnatomyType <- ifelse(UniFrac$Anatomy1 == UniFrac$Anatomy2, "Same Anatomy","Different Anatomy")

UniFrac$Comparisonf <- ifelse(UniFrac$Individual1 == UniFrac$Individual2, "Within Individuals", "Between Individuals")
UniFrac$Comparisonf <- ifelse((UniFrac$Comparisonf == "Between Individuals") & (UniFrac$Household1 == UniFrac$Household2), "Between Individuals Within Households", UniFrac$Comparisonf)
UniFrac$Comparisonf <- ifelse(UniFrac$Comparisonf == "Between Individuals", "Between Households",UniFrac$Comparisonf)
```

### PCoA

```jsx
library(readr)
library(dplyr)
library(tidyr)
library(wilkoxmisc)
library(reshape2)
source("AICc.PERMANOVA2.r")

tax <- read_tsv("kraken_species_abundance.tidy.txt")
tax <- tax %>% filter(!Species == "unclassified")
tax <- as.data.frame(tax)
Meta <- read_tsv("metadata_doorknob_skin.txt")
tax2 <- merge(tax, Meta, by = "Sample", all.x = TRUE)
tax3 <- tax2 %>% filter(Location == "Residence 3")
Tax <- tax3[,1:4]

tax_cast <- dcast(Tax, Sample ~ Species, value.var = "count", fill = 0)
colnames(tax_cast) <- c("Sample", paste0("Species",1:4713))

write_tsv(tax_cast, "kraken_species_abundance_pcoa_rename_residence3.cast.txt")

#replace " " with "_" before loading the table

tax_cast <- kraken_species_abundance_pcoa_rename_residence3.cast
install.packages("ecodist")
library(ecodist)
distance <- distance(tax_cast, method = "bray-curtis") # in R package "ecodist
library(vegan)
distance1 <- vegdist(tax_cast, method="bray")
distance2 <- vegdist(tax_cast, method="jaccard")

UniFracMatrix1 <- as.matrix(distance1)
write.table(UniFracMatrix,"kraken_species_bray_cutis_matrix_residence3_bray.txt", sep="\t", col.names=TRUE, row.names=TRUE)

UniFracMatrix2 <- as.matrix(distance2)
write.table(UniFracMatrix,"kraken_species_bray_cutis_matrix_residence3_jaccard.txt", sep="\t", col.names=TRUE, row.names=TRUE)

Meta <- read_tsv("metadata_doorknob_skin.txt")
Meta2 <- Meta %>% filter(Location == "Residence 3")

#adonis <- adonis2(UniFracMatrix~Surface_material + Sample_type + Group + Location + Sample_type:Location,data=Meta)

PCoA_bray <- cmdscale(UniFracMatrix1, k = 2, eig = TRUE)
PCoA_jaccard <- cmdscale(UniFracMatrix2, k = 2, eig = TRUE)
DF1 <- data.frame(Sample = row.names(PCoA_bray$points), PCoA1 = PCoA_bray$points[,1], PCoA2 = PCoA_bray$points[,2], row.names = NULL)
DF2 <- data.frame(Sample = row.names(PCoA_jaccard$points), PCoA1 = PCoA_jaccard$points[,1], PCoA2 = PCoA_jaccard$points[,2], row.names = NULL)
PCoA_1 <- merge(DF1, Meta2, by = "Sample", all.x = TRUE)
PCoA_2 <- merge(DF2, Meta2, by = "Sample", all.x = TRUE)

Eigenvalues <- eigenvals(PCoA)
VarianceExplained <- Eigenvalues /sum(Eigenvalues)
VarianceExplained1 <- 100 * signif(VarianceExplained[1], 2)
VarianceExplained2 <- 100 * signif(VarianceExplained[2], 2)

#Merge Dataframe with metadata
AllSamples <-data.frame(Sample = row.names(UniFracMatrix1))
AllSamples <- merge(AllSamples, Meta, by = "Sample", all.x = TRUE)
subsample <- AllSamples %>% filter(Location == "Residence 3")
#Check that UniFrac matrix rows match samples table (for ANOSIM)
sum(row.names(UniFracMatrix)==AllSamples$Sample) == length(AllSamples$Sample)

install.packages("AICcPermanova")
library(AICcPermanova)
#Perform ANOSIM
BrayCurtis <- as.dist(distance)
ANOSIM <- anosim(BrayCurtis, grouping = AllSamples$Sample_type)
#See ANOSIM options
ls(ANOSIM)

#Load ANOSIM statistics and signif
ANOSIM$statistic
# [1] 0.06960386
ANOSIM$signif
# [1] 0.014

newdata <- na.omit(AllSamples)
#Perform Adonis(permanova)
source("AICc.PERMANOVA2.r")
Adonis <- adonis2(BrayCurtis~Sampling_day:Sample_type, data=AllSamples)

#best model for nested permanova (change variable)
Adonis <- adonis2(BrayCurtis~Surface_material + Location + Sampling_day:Sample_type,data=AllSamples)

AICc_permanova2(Adonis)$AIC
[1] -331.5818

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = BrayCurtis ~ Sampling_day:Sample_type, data = AllSamples)
Df SumOfSqs      R2        F Pr(>F)
adonis2(formula = BrayCurtis ~ Surface_material + Location + Sampling_day:Sample_type, data = AllSamples)
Df SumOfSqs      R2       F Pr(>F)    
Surface_material           1    1.208 0.02686 10.9671  0.001 ***
Location                   3   22.010 0.48947 66.6163  0.001 ***
Sampling_day:Sample_type  58    6.661 0.14813  1.0428  0.292    
Residual                 137   15.089 0.33554                   
Total                    199   44.968 1.00000                   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

dbrda <- dbrda(BrayCurtis ~ Sample_type+Sampling_day, data = AllSamples)
anova(dbrda,by="margin")

Adonis <- adonis(BrayCurtis~Location, data=AllSamples, permutations=999)
Call:
adonis(formula = BrayCurtis ~ Location, data = AllSamples, permutations = 999)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Location    8     7.715 0.96444  5.7046 0.21564  0.001 ***
Residuals 166    28.064 0.16906         0.78436
Total     174    35.780                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Adonis <- adonis(BrayCurtis~Sample_type, data=AllSamples, permutations=999)
Call:
adonis(formula = BrayCurtis ~ Sample_type, data = AllSamples,      permutations = 999)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Sample_type   3     8.388 2.79615  17.456 0.23445  0.001 ***
Residuals   171    27.391 0.16018         0.76555
Total       174    35.780                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Adonis <- adonis(BrayCurtis~Surface_material, data=AllSamples, permutations=999)
Call:
adonis(formula = BrayCurtis ~ Surface_material, data = AllSamples,      permutations = 999)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Surface_material   1     6.139  6.1389   35.83 0.17157  0.001 ***
Residuals        173    29.641  0.1713         0.82843
Total            174    35.780                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Adonis <- adonis(BrayCurtis~Group, data=AllSamples, permutations=999)
Call:
adonis(formula = BrayCurtis ~ Group, data = AllSamples, permutations = 999)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Group       1     2.892  2.8918  15.211 0.08082  0.001 ***
Residuals 173    32.888  0.1901         0.91918
Total     174    35.780                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Specify X and Y axis
#bray-curtis
xlab1 <- paste0("PCoA1 (Variance explained: 43%)")
ylab1 <- paste0("PCoA2 (Variance explained: 24%)")
#jaccard
xlab2 <- paste0("PCoA1 (Variance explained: 35%)")
ylab2 <- paste0("PCoA2 (Variance explained: 16%)")

#Draw plot by sample types
install.packages("ggthemes")
library(ggthemes)

PCoA_1$Sampling_day <- factor(PCoA_1$Sampling_day, levels = c("Day 1 day","Day 1 night","Day 2 day","Day 2 night","Day 3 day","Day 3 night","Day 4 day","Day 4 night","Day 5 day","Day 5 night","Day 6 day","Day 6 night","Day 7 day","Day 7 night","Day 8 day","Day 8 night","Day 9 day","Day 9 night","Day 10 day","Day 10 night"))
Plot1 <- ggplot(PCoA_1, aes(x = PCoA1, y = PCoA2, color=Sample_type))
Plot1 <- Plot1 + geom_point(size=4,alpha = 0.8)
Plot1 <- Plot1 + stat_ellipse(geom = "polygon",alpha = 0.05, aes(color = Sample_type))
cbPalette <- c("#2dbaff","#a75cfb","#fd289d")
Plot1 <- Plot1 + scale_color_manual(values = cbPalette)
Plot1 <- Plot1 + xlab(xlab1)
Plot1 <- Plot1 + ylab(ylab1)
Plot1 <- Plot1 + theme_bw()
#Plot <- Plot + scale_color_brewer(palette = "Paired")
Plot1 <- Plot1 + theme(axis.title=element_text(size=12))
Plot1 <- Plot1 + labs(title = "Bray Curtis PCoA by Sample Type")
Plot1 <- Plot1 + theme(axis.text=element_text(size=10))
ggsave("bray_curtis_pcoa_by_sample_type_ellipse_2.png", width = 7, height = 7)

PCoA_2$Sampling_day <- factor(PCoA_2$Sampling_day, levels = c("Day 1 day","Day 1 night","Day 2 day","Day 2 night","Day 3 day","Day 3 night","Day 4 day","Day 4 night","Day 5 day","Day 5 night","Day 6 day","Day 6 night","Day 7 day","Day 7 night","Day 8 day","Day 8 night","Day 9 day","Day 9 night","Day 10 day","Day 10 night"))
Plot2 <- ggplot(PCoA_2, aes(x = PCoA1, y = PCoA2, color=Sample_type))
Plot2 <- Plot2 + geom_point(size=4,alpha = 0.8)
Plot2 <- Plot2 + stat_ellipse(geom = "polygon",alpha = 0.05, aes(color = Sample_type))
cbPalette <- c("#2dbaff","#a75cfb","#fd289d")
Plot2 <- Plot2 + scale_color_manual(values = cbPalette)
Plot2 <- Plot2 + xlab(xlab2)
Plot2 <- Plot2 + ylab(ylab2)
Plot2 <- Plot2 + theme_bw()
#Plot <- Plot + scale_color_brewer(palette = "Paired")
Plot2 <- Plot2 + theme(axis.title=element_text(size=12))
Plot2 <- Plot2 + labs(title = "Jaccard PCoA by Sample Type")
Plot2 <- Plot2 + theme(axis.text=element_text(size=10))
ggsave("jaccard_pcoa_by_sample_type_ellipse_2.png", width = 7, height = 7)

library(ggpubr)
library(gridExtra)
install.packages("cowplot")
library(cowplot)
Figure <- plot_grid(Plot1, Plot2, labels = LETTERS[1:2], ncol = 2)
ggsave(plot = Figure, "pier_samples_kraken_subsamp_shannon_boxplot_combined_by_sample_type_location.png", width = 14, height = 7)
```

## Sloan Neutral Model

### Calculatation

```jsx
library(readr)
library(Hmisc)
library(bbmle)
library(wilkoxmisc)
library(reshape2)

table <- read_tsv("kraken_report_all_table_rarefied.txt")
meta <- read_tsv("metadata_doorknob_skin.txt") %>% select(Sample,Type,Location)

merge <- left_join(table,meta)

sub_table <- merge %>%
    filter(Type == "Indoor surface") %>%
    filter(Location == "Residence 3") %>%
    select(-Type,-Sample,-Location)

write_tsv(sub_table,"kraken_report_all_table_rarefied_doorknob_residence3.txt")

pool <- read_tsv("kraken_report_all_table_rarefied_skin_residence3.txt")
spp <- read_tsv("kraken_report_all_table_rarefied_doorknob_residence3.txt")

#If stats=FALSE the function will return a table of observed and predicted values for each otu.
sncm.fit <- function(spp, pool=TRUE, stats=FALSE, taxon=NULL){
	require(minpack.lm)
	require(Hmisc)
	require(stats4)
	
	options(warn=-1)

	#Calculate the number of individuals per community
	N <- mean(apply(spp, 1, sum))
    
	
	#Calculate the average relative abundance of each taxa across communities
	if(is.null(pool)){
		p.m <- apply(spp, 2, mean)
		p.m <- p.m[p.m != 0]
		p <- p.m/N
	} else {
		p.m <- apply(pool, 2, mean)
		p.m <- p.m[p.m != 0]
		p <- p.m/N
	}

	#Calculate the occurrence frequency of each taxa across communities
	spp.bi <- 1*(spp>0)
	freq <- apply(spp.bi, 2, mean)
	freq <- freq[freq != 0]

	#Combine
	C <- merge(p, freq, by=0)
	C <- C[order(C[,2]),]
	C <- as.data.frame(C)
	C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
	p <- C.0[,2]
	freq <- C.0[,3]
	names(p) <- C.0[,1]
	names(freq) <- C.0[,1]

	#Calculate the limit of detection
	d = 1/N

	##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
	m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
	m.ci <- confint(m.fit, 'm', level=0.95)
	
	##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
	sncm.LL <- function(m, sigma){
		R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
		R = dnorm(R, 0, sigma)
		-sum(log(R))
	}
	m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
	
	##Calculate Akaike's Information Criterion (AIC)
	aic.fit <- AIC(m.mle, k=2)
	bic.fit <- BIC(m.mle)

	##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
	freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
	Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
	
	pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Calculate AIC for binomial model
	bino.LL <- function(mu, sigma){
		R = freq - pbinom(d, N, p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
	
	aic.bino <- AIC(bino.mle, k=2)
	bic.bino <- BIC(bino.mle)
	
	##Goodness of fit for binomial model
	bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
	Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))

	bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Calculate AIC for Poisson model
	pois.LL <- function(mu, sigma){
		R = freq - ppois(d, N*p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
	
	aic.pois <- AIC(pois.mle, k=2)
	bic.pois <- BIC(pois.mle)
	
	##Goodness of fit for Poisson model
	pois.pred <- ppois(d, N*p, lower.tail=FALSE)
	Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))

	pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

	##Results
	if(stats==TRUE){
		fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
		fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
		return(fitstats)
	} else {
		A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
		A <- as.data.frame(A)
		colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
		if(is.null(taxon)){
			B <- A[order(A[,1]),]
		} else {
			B <- merge(A, taxon, by=0, all=TRUE)
			row.names(B) <- B[,1]
			B <- B[,-1]
			B <- B[order(B[,1]),]
		}
		return(B)
	}
}

stats.otu <- sncm.fit(spp,pool)
write.table(stats.otu, "skin_to_doorknob_residence3.txt", sep="\t", col.names=NA)

#If stats=TRUE the function will return fitting statistics.
sncm.fit <- function(spp, pool=TRUE, stats=TRUE, taxon=NULL){
	require(minpack.lm)
	require(Hmisc)
	require(stats4)
	
	options(warn=-1)

	#Calculate the number of individuals per community
	N <- mean(apply(spp, 1, sum))
	
	#Calculate the average relative abundance of each taxa across communities
	if(is.null(pool)){
		p.m <- apply(spp, 2, mean)
		p.m <- p.m[p.m != 0]
		p <- p.m/N
	} else {
		p.m <- apply(pool, 2, mean)
		p.m <- p.m[p.m != 0]
		p <- p.m/N
	}

	#Calculate the occurrence frequency of each taxa across communities
	spp.bi <- 1*(spp>0)
	freq <- apply(spp.bi, 2, mean)
	freq <- freq[freq != 0]

	#Combine
	C <- merge(p, freq, by=0)
	C <- C[order(C[,2]),]
	C <- as.data.frame(C)
	C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
	p <- C.0[,2]
	freq <- C.0[,3]
	names(p) <- C.0[,1]
	names(freq) <- C.0[,1]

	#Calculate the limit of detection
	d = 1/N

	##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
	m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
	m.ci <- confint(m.fit, 'm', level=0.95)
	
	##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
	sncm.LL <- function(m, sigma){
		R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
		R = dnorm(R, 0, sigma)
		-sum(log(R))
	}
	m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
	
	##Calculate Akaike's Information Criterion (AIC)
	aic.fit <- AIC(m.mle, k=2)
	bic.fit <- BIC(m.mle)

	##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
	freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
	Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
	
	pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Calculate AIC for binomial model
	bino.LL <- function(mu, sigma){
		R = freq - pbinom(d, N, p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
	
	aic.bino <- AIC(bino.mle, k=2)
	bic.bino <- BIC(bino.mle)
	
	##Goodness of fit for binomial model
	bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
	Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))

	bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Calculate AIC for Poisson model
	pois.LL <- function(mu, sigma){
		R = freq - ppois(d, N*p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
	
	aic.pois <- AIC(pois.mle, k=2)
	bic.pois <- BIC(pois.mle)
	
	##Goodness of fit for Poisson model
	pois.pred <- ppois(d, N*p, lower.tail=FALSE)
	Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))

	pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

	##Results
	if(stats==TRUE){
		fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
		fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
		return(fitstats)
	} else {
		A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
		A <- as.data.frame(A)
		colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
		if(is.null(taxon)){
			B <- A[order(A[,1]),]
		} else {
			B <- merge(A, taxon, by=0, all=TRUE)
			row.names(B) <- B[,1]
			B <- B[,-1]
			B <- B[order(B[,1]),]
		}
		return(B)
	}
}

stats.all <- sncm.fit(spp,pool)
write_tsv(stats.all, "skin_to_doorknob_residence3_stats.txt")
```

### R visulization

```jsx
library(readr)
library(dplyr)
library(ggplot2)

skin <- read_tsv("skin_to_doorknob_residence3.txt")
names(skin)[1] <- c("Species")
skin <- skin %>% mutate(Partition=ifelse(freq < pred.lwr, "Below", "Neutral"))
skin <- skin %>% mutate(Partition=ifelse(freq > pred.upr, "Above", Partition))
write_tsv(skin,"skin_to_doorknob_residence3.txt")

spp <- read_tsv("kraken_report_all_table_rarefied_skin_residence3.txt")
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p <- as.data.frame(p)
colnames(p) <- c("mean_abundance")
write.table(p, "kraken_species_mean_abundance_skin_residence3.txt", sep="\t", row.names=TRUE)
abundance <- read_tsv("kraken_species_mean_abundance_skin_residence3.txt")
names(abundance)[1] <- c("Species")
merge <- left_join(skin, abundance)
write_tsv(merge, "skin_to_doorknob_residence3_merge.txt")

#Plot the Sloan neutral model for all individuals
cbPalette <- c("#537bff", "#ff3c6d", "#b248ff")
Table <- read_tsv("skin_to_doorknob_residence3.txt") %>% filter(!Species == "unclassified")
Plot <- ggplot(Table,aes(x=log(mean_abundance)))
Plot <- Plot + geom_point(aes(y=freq, colour=Partition), size=1.5,alpha=0.8)
Plot <- Plot + geom_line(aes(y=freq.pred))
Plot <- Plot + geom_ribbon(aes(ymin=pred.lwr, ymax=pred.upr),fill="#969696",alpha=0.5)
Plot <- Plot + scale_color_manual(values = cbPalette)
Plot <- Plot + theme_classic()
Plot <- Plot + theme(legend.title = element_blank())
Plot <- Plot + theme(plot.title = element_text(hjust = 0.5, face="bold"))
Plot <- Plot + ggtitle("Residence 3")
Plot <- Plot + xlab(paste0("log (mean relative abundance of species on skin)")) +ylab(paste0("Occurrence frequency of species on doorknob"))
Plot
ggsave("sloan_neutral_model_skin_to_doorknob_residence3.png", width=12, height=7)
```

![sloan_neutral_model_skin_to_doorknob_residence3.png](Data%20Analysis%20Script%20and%20Output/sloan_neutral_model_skin_to_doorknob_residence3.png)

## Drep

[Quick Start - drep 2.0.0 documentation](https://drep.readthedocs.io/en/latest/quick_start.html)

[drep：微生物基因组快速去冗余-文章解读+帮助文档+实战_刘永鑫Adam的博客-CSDN博客](https://blog.csdn.net/woodcorpse/article/details/108557761)

```jsx
[zengyg@lnxbio zengyg]$ cat dreptestplus.slurm
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

dRep dereplicate dreptestplus/ -comp 50 -g /data/bio/zengyg/sample_genome_assemblies_genome_fasta/*.fa
```

### R visualization

[R语言ggplot2包绘制散点图详解](https://www.jianshu.com/p/49357878c555)

[R 数据可视化 —— ggplot 图例设置](https://zhuanlan.zhihu.com/p/362736163)

[R语言 绘图组合布局  一页多图_Nobeli的博客-CSDN博客_r语言多图布局](https://blog.csdn.net/ysl416/article/details/106266467)

```jsx
# Load library
library(ggplot2) 
library(dplyr)

#Import data
dRep_rMAGs <- read.csv("E:/BIO303_FYP/Data_Analysis/dRep_rMAGs.csv")
dRep_rMAGs_filter <- dRep_rMAGs %>% filter(Genotype == "Cluster A")

#Plot
plot <- ggplot(data = dRep_rMAGs_filter)
plot <- plot + geom_point(mapping=aes(x = factor(Sampling_day,levels = c("Day 1 day","Day 1 night","Day 2 day","Day 2 night","Day 3 day","Day 3 night","Day 4 day","Day 4 night","Day 5 day","Day 5 night","Day 6 day","Day 6 night","Day 7 day","Day 7 night","Day 8 day","Day 8 night","Day 9 day","Day 9 night","Day 10 day","Day 10 night")), y=Sample_type, colour = factor(Sample_type), size = 2))
plot <- plot + labs(x = "Sampling day", y = "Sample type", colour = "Sample type")
cbPalette <- c("#7998ff", "#ff6780", "#c371ff")
plot <- plot + scale_color_manual(values = cbPalette)
plot <- plot + facet_wrap(~classification) # arrange based on classification
plot <- plot + theme(axis.text.x = element_text(hjust=0,angle=90)) #change the direction of x axis
plot <- plot + guides(size = "none") #remove size legend
plot

#Export figure
ggsave("taxonomy_day.png", width=12, height=7)
```

### Results

![taxonomy_day.png](Data%20Analysis%20Script%20and%20Output/taxonomy_day.png)

![taxonomy_day.png](Data%20Analysis%20Script%20and%20Output/taxonomy_day%201.png)

## CoverM

```jsx
#压缩文件为gz
[zengyg@lnxbio read2]$ cat gzip2.slurm 
#!/bin/bash
#SBATCH -J gzip
#SBATCH -p cpunode
#SBATCH -o gzip-%j.out
#SBATCH -e gzip-%j.err
#SBATCH -n 8

export LC_ALL="en_US.UTF-8"

ml load kraken2-2.1.1-gcc-8.5.0-5vbredn
ml load parallel-20210922-gcc-8.5.0-v7ilxoy

cat initial_samples.txt | parallel -j 1 gzip ./{}_repaired_paired_2.fastq
```

[coverm genome usage](https://wwood.github.io/CoverM/coverm-genome.html)

[coverm 0.6.1 - Docs.rs](https://docs.rs/crate/coverm/latest)

```jsx
[zengyg@lnxbio zengyg]$ cat covermtest.slurm
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

cat initial_samples.txt | parallel -j 1 coverm genome --coupled /data/bio/zengyg/Coverm_Samples/read1/{}_repaired_paired_1.fastq.gz /data/bio/zengyg/Coverm_Samples/read2/{}_repaired_paired_2.fastq.gz --genome-fasta-files /data/bio/zengyg/Drep/drep_g_agrococcus/dereplicated_genomes/SRR12893903_bin5.fa -o output.tsv
```

### R visualization

```jsx
#Install and Load Library
install.packages("devtools")
install.packages("stringr")
library(devtools)
install_github('wilkox/wilkoxmisc')
install.packages("starsExtra")
install.packages("tidyverse")

library(wilkoxmisc)
library(reshape2)
library(dplyr)
library(readr)
library(ggplot2)
library(starsExtra)
library(tidyverse)
library(readxl)

#Merge Data
Meta <- metadata_doorknob_skin
Meta2 <- Meta %>% filter(Location == "Residence 3")
Match <- match
MergeData <- merge(Meta2, Match, by = "Sample", all.x = TRUE)
write_tsv(MergeData, "Information.txt")
coverm_output <- coverm_results

#Plot
plot <- ggplot(data = coverm_output) 
plot <- plot + geom_point(mapping=aes(x = factor(Sampling_day,levels = c("Day 1 day","Day 1 night","Day 2 day","Day 2 night","Day 3 day","Day 3 night","Day 4 day","Day 4 night","Day 5 day","Day 5 night","Day 6 day","Day 6 night","Day 7 day","Day 7 night","Day 8 day","Day 8 night","Day 9 day","Day 9 night","Day 10 day","Day 10 night")), y = Sample_type, colour = Sample_type, size = RelativeAbundance))
plot <- plot + labs(x = "Sampling day", y = "Sample type", colour = "Sample type", size = "RelativeAbundance(%)")
plot <- plot + facet_wrap(~classification) # arrange based on classification
plot <- plot + theme(axis.text.x = element_text(hjust=0,angle=90)) #change the direction of x axis
#plot <- plot + guides(size = "none") #remove size legend
plot

#Export figure
ggsave("coverM_transmission.png", width=12, height=7)
```

![coverM_transmission.png](Data%20Analysis%20Script%20and%20Output/coverM_transmission.png)

## Resfams

### prokka转化（Metawrap)

[Prokka：快速原核基因组、宏基因组基因注释_刘永鑫Adam的博客-CSDN博客_prokka](https://blog.csdn.net/woodcorpse/article/details/101364994)

[原核生物基因组快速注释——Prokka](https://zhuanlan.zhihu.com/p/379981540)

```jsx
[zengyg@lnxbio prokka]$ cat prokka.sh 
#!/bin/bash
#SBATCH -J Prokka
#SBATCH -p cpunode
#SBATCH -o prokka-%j.out
#SBATCH -e prokka-%j.err
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
ml load py-biopython-1.73-gcc-4.8.5-qfbrlp5
ml load py-checkm-genome-1.0.13-gcc-8.5.0-sdtawxc
ml load prodigal-2.6.3-gcc-8.5.0-ixjx62u
ml load prokka-1.14.6-gcc-8.5.0-iveilfb
ml load tbl2asn-2022-04-26-gcc-8.5.0-iqp3oyy

prokka –setupdb

metawrap annotate_bins -o FUNCT_ANNOT -b rMAG_bins/
```

### Resfams

[Resfams - Dantas Lab](http://www.dantaslab.org/resfams)

```jsx
#下载Resfam数据库
wget http://dantaslab.wustl.edu/resfams/Resfams.hmm.gz
#解压
gunzip Resfams.hmm.gz
```

[使用 HMMER 进行 PFAM 注释](http://www.chenlianfu.com/?p=2274)

### HMMER3

[2016_Boolchandani_Patel_FunctionalMetagenomicsAntibioticResistanceProtocol_MethodsMolBiol+(1).pdf](Data%20Analysis%20Script%20and%20Output/2016_Boolchandani_Patel_FunctionalMetagenomicsAntibioticResistanceProtocol_MethodsMolBiol(1).pdf)

![微信图片_20230228155219.png](Data%20Analysis%20Script%20and%20Output/%25E5%25BE%25AE%25E4%25BF%25A1%25E5%259B%25BE%25E7%2589%2587_20230228155219.png)

```jsx
ml load ml load hmmer-3.3.2-gcc-8.5.0-r4vrjwg

hmmpress Resfams.hmm

hmmscan --cut_ga --tblout hmmtest1.tsv Resfams.hmm /data/bio/zengyg/FUNCT_ANNOT/bin_translated_genes/SRR12893950_bin5.faa

awk '$5 < 0.001 {print $1"\t"$3"\t"$4}' hmm_g_B.tsv > hmm_g_B_filter.tsv
```

## Plasflow https://github.com/smaegol/PlasFlow

## mobileOG

[适合宏基因组的可移动元件数据库mobileOG的使用_刘永鑫Adam的博客-CSDN博客](https://metagenome.blog.csdn.net/article/details/128607914?spm=1001.2101.3001.6650.3&utm_medium=distribute.pc_relevant.none-task-blog-2~default~YuanLiJiHua~Position-3-128607914-blog-108257592.pc_relevant_3mothn_strategy_recovery&depth_1-utm_source=distribute.pc_relevant.none-task-blog-2~default~YuanLiJiHua~Position-3-128607914-blog-108257592.pc_relevant_3mothn_strategy_recovery&utm_relevant_index=4)

```jsx
ml load r-4.1.3-gcc-8.5.0-tdeuqoe
ml load cmake-3.23.1-gcc-8.5.0-v7krosh
ml load libxml2-2.9.13-gcc-8.5.0-5ul2v2s
ml load curl-7.83.0-gcc-8.5.0-v3scdyk
ml load r-jpeg-0.1-9-gcc-8.5.0-jjptfjn

# 使用diamond构建索引
ml load diamond-2.0.15-gcc-8.5.0-ltico7y
diamond makedb --in ./mobileOG-db_beatrix-1.6.All.faa -d ./mobileOG-db-beatrix-1.X.dmnd

# 直接使用diamond进行序列比对：
mkdir /data/bio/zengyg/mobileOG/
diamond blastp -q /data/bio/zengyg/FUNCT_ANNOT/bin_translated_genes/SRR12893831_bin1.faa --db /data/bio/zengyg/mobileOG-db/mobileOG-db-beatrix-1.X.dmnd --outfmt 6 stitle qtitle pident bitscore slen evalue qlen sstart send qstart qend -k 15 -o /data/bio/zengyg/mobileOG/SRR12893831_bin1_mobileOG.tsv -e 1e-20 --query-cover 90 --id 90
 
# 作者写的py程序进行整理
ml load python-3.9.12-gcc-8.5.0-ews4x7x
python /data/bio/zengyg/mobileOG-db/mobileOG-pl/mobileOGs-pl-kyanite.py --o /data/bio/zengyg/mobileOG/SRR12893831_bin1 --i /data/bio/zengyg/mobileOG/SRR12893950_bin14_mobileOG2.tsv -m /data/bio/zengyg/mobileOG-db/mobileOG-db-beatrix-1.6-All.csv

# 修改一下吧，作者提供的工具与分析习惯不同。
diamond blastp --db /data/bio/zengyg/mobileOG-db/mobileOG-db-beatrix-1.X.dmnd --query /data/bio/zengyg/FUNCT_ANNOT/bin_translated_genes/SRR12893950_bin4.faa --outfmt 6 --max-target-seqs 1 -e 1e-5 --sensitive --out /data/bio/zengyg/mobileOG/SRR12893950_bin14_mobileOG2.tsv
 
# 提取基因对应基因家族
cut -f 1,2 /data/bio/zengyg/mobileOG/SRR12893950_bin14_mobileOG2.tsv | uniq | \
      sed '1 i Name\tResGeneID' > data/bio/zengy/mobileOG/gene_mobileOG.list
```

## MEGfinder https://github.com/bhattlab/MGEfinder

[[文献解读#1] 可移动遗传元件助力细菌适应性进化](https://zhuanlan.zhihu.com/p/201110975)

[User Manual](https://github.com/bhattlab/MGEfinder/wiki/User-manual)

[A Bioinformatic Analysis of Integrative Mobile Genetic Elements Highlights Their Role in Bacterial Adaptation](https://www.sciencedirect.com/science/article/pii/S1931312819305463)

[Installation — Cutadapt 4.3 documentation](https://cutadapt.readthedocs.io/en/stable/installation.html)

```jsx
ml load mgefinder 3seq-1.7

# FASTQ files were then deduplicated using the SuperDeduper command from the HTStream toolset with the following command:
ml  load  htstream_v1.3.0
hts_SuperDeduper -1 SRR12893849_repaired_paired_1.fastq.gz -2 SRR12893849_repaired_paired_2.fastq.gz -f SRR12893849.nodup

# FASTQ files were trimmed using the Trime Galore package and the command:
ml load trimgalore-0.6.6-gcc-8.5.0-pctdye3
ml load py-cutadapt-2.10-gcc-8.5.0-jtye6k6
trim_galore --fastqc --paired SRR12893849.nodup_R1.fastq.gz SRR12893849.nodup_R2.fastq.gz

# The reference genome was indexed for alignment with BWA MEM and the alignment was performed using the commands:
ml load bwa-0.7.17-gcc-8.5.0-oo4toav 
ml load bwa-mem2
bwa index dietzia_cinnamea_refseq.fasta
bwa mem NC_013441.1.fasta SRR12893912.nodup_R1_val_1.fq.gz SRR12893912.nodup_R2_val_2.fq.gz > SRR12893912.NC_013441.1.sam

# There is a built-in command in mgefinder, called mgefinder formatbam which takes this SAM file and prepares it for analysis by mgefinder:
[zengyg@lnxbio zengyg]$ cat forma.sh 
#!/bin/bash
#SBATCH -J mgefinder
#SBATCH -p cpunode
#SBATCH -o mgefinder-%j.out
#SBATCH -e mgefinder-%j.err
#SBATCH -n 8

export LC_ALL="en_US.UTF-8"

export SINGULARITY_BINDPATH="/data/bio/zengyg"

ml load  mgefinder 3seq-1.7

singularity exec mgefinder.sif mgefinder formatbam /data/bio/zengyg/SRR12893912.NC_013441.1.sam SRR12893912.NC_013441.1.bam

generate the assembled contig file SRR12893849.fna, you can use the SPAdes command:
ml load spades-3.15.3-gcc-8.5.0-7lafyz5
spades.py -1 SRR12893912.nodup_R1_val_1.fq.gz -2 SRR12893912.nodup_R2_val_2.fq.gz -o SRR12893912
```

```jsx
[zengyg@lnxbio zengyg]$ cat mgefindertest.sh 
#!/bin/bash
#SBATCH -J mgefinder
#SBATCH -p cpunode
#SBATCH -o mgefinder-%j.out
#SBATCH -e mgefinder-%j.err
#SBATCH -n 8

export LC_ALL="en_US.UTF-8"

export SINGULARITY_BINDPATH="/data/bio/zengyg"

ml load  mgefinder 3seq-1.7

singularity exec /data/bio/sif/mgefinder.sif mgefinder workflow denovo /data/bio/zengyg/mgefinder_test_workdir/
```

本教程为这 10 个屎肠球菌分离株生成了候选综合移动遗传元件列表，并根据屎肠球菌参考基因组对插入进行了基因分型。接下来采取的步骤将在很大程度上取决于您要执行的分析类型。如果您想在这些候选插入中搜索抗生素抗性基因，您可以将`04.makefasta.efae_GCF_900639545.all_seqs.fna`FASTA 文件上传到[ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/)。[如果您希望识别转座酶，我们推荐使用ISEScan](https://www.ncbi.nlm.nih.gov/pubmed/29077810)提供的 HMM 。如果您想识别噬菌体元素，可以将`04.makefasta.efae_GCF_900639545.repr_seqs.fna`文件上传到[PHASTER](http://phaster.ca/)。

### Another way (use bin)

```jsx
gzip SRR12893912_bin5.fq

ml load bwa-0.7.17-gcc-8.5.0-oo4toav 
ml load bwa-mem2
bwa index NC_013441.fna

bwa mem NC_013441.fna SRR12893912_bin5.fq.gz > SRR12893912.NC_013441.sam

[zengyg@lnxbio zengyg]$ cat forma.sh 
#!/bin/bash
#SBATCH -J formatbam
#SBATCH -p cpunode
#SBATCH -o formatbam-%j.out
#SBATCH -e formatbam-%j.err
#SBATCH -n 8

export LC_ALL="en_US.UTF-8"

export SINGULARITY_BINDPATH="/data/bio/zengyg"

ml load  mgefinder 3seq-1.7

singularity exec mgefinder.sif mgefinder formatbam /data/bio/zengyg/SRR12893912.NC_013441.sam SRR12893912.NC_013441.bam

```

### ResFinderhttps://github.com/cadms/resfinder

[ResFinder – an open online resource for identification of antimicrobial resistance genes in next-generation sequencing data and prediction of phenotypes from genotypes](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000748)

[神器在手耐药基因预测分析不难！！！_http](https://www.sohu.com/a/404994663_464200)

### ISEScan

```jsx
ml   load  python-3.9.12-gcc-8.5.0-ews4x7x
ml   load  isescan-1.7.2.1-gcc-8.5.0-x6j4wpg
isescan.py rMAG_bins/SRR12893912_bin5.fa isescan_results/ isescan_hmm/

isescan.py [-h] [--version] [--removeShortIS] [--no-FragGeneScan] [--nthread NTHREAD] seqfile path2proteome path2hmm
```

[佳作推荐|ISEScan: 自动化识别原核基因组中的插入序列-行业动态-广东美格基因科技有限公司](http://www.magigene.com/Article/show/670.html)

[https://github.com/xiezhq/ISEScan](https://github.com/xiezhq/ISEScan)

[ISEScan: automated identification of insertion sequence elements in prokaryotic genomes](https://academic.oup.com/bioinformatics/article/33/21/3340/3930124?login=true)

### PHAST

[PHAST](http://phast.wishartlab.com/)

[http://phast.wishartlab.com/cgi-bin/Results.cgi?num=1679576159&multi=1](http://phast.wishartlab.com/cgi-bin/Results.cgi?num=1679576159&multi=1)

### oriTfinder

[oriTfinder: a web-based tool for the identification of origin of transfers in DNA sequences of bacterial mobile genetic elements](https://academic.oup.com/nar/article/46/W1/W229/4992657?login=true)

### mobile genetic elements

[https://www.nature.com/articles/nrmicro1235/#:~:text=Mobile genetic elements (MGEs) are segments of DNA,(intracellular mobility) or between bacterial cells (intercellular mobility)](https://www.nature.com/articles/nrmicro1235/#:~:text=Mobile%20genetic%20elements%20%28MGEs%29%20are%20segments%20of%20DNA,%28intracellular%20mobility%29%20or%20between%20bacterial%20cells%20%28intercellular%20mobility%29).

[环境科学技术论文：基于宏基因组学技术分析污水处理系统中耐药基因和可移动遗传元件分布与丰度](https://www.docin.com/p-616674147.html)

## 病原微生物数据库

[Pathogenic Bacteria – Database](https://globalrph.com/bacteria/)

[NCBI - WWW Error Blocked Diagnostic](https://www.ncbi.nlm.nih.gov/pathogens/)

感觉这个比较靠谱

[](https://www.nprc.org.cn/)

[Global Catalogue of Microorganisms Pathogenic Microorganism Database](https://nmdc.cn/gcpathogen/)

## 组间网络分析

[在线作图|绘制组间网络分析(Network Analysis)](https://www.jianshu.com/p/13864dca130e)