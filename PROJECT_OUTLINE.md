# Project Summary

This document is a single-file walkthrough of the project pipeline, code, and analysis logic behind my undergraduate microbial metagenomics thesis. It is meant to serve as a reviewer-friendly portfolio artifact: someone can read this one file and see the full technical story without needing to open separate scripts.

The project focused on microbial community structure and transmission across household indoor surfaces and skin samples using metagenomic sequencing data. The work combined HPC job orchestration, metagenomics tools, downstream ecological analysis, and figure generation in R.

# What This File Contains

- the biological and computational goal of the project
- the end-to-end pipeline in execution order
- the main shell, Slurm, Python, and R code used at each stage
- examples of the downstream statistical and visualization workflow
- notes about outputs and interpretation

# Computing Context

Most of the pipeline was executed in an HPC environment using Slurm batch jobs and environment modules. For that reason, many paths are absolute and many commands begin with `ml load ...`. Those details are intentionally preserved because they show the original execution environment and workflow orchestration.

# Project Goal

The central aim of this project was to analyze metagenomic samples from skin and indoor surfaces in order to:

1. clean and preprocess sequencing reads
2. remove host-associated contamination
3. assemble and bin genomes
4. classify taxonomic composition
5. evaluate genomic quality and redundancy
6. study abundance patterns, similarity structure, and potential transmission signals

# End-to-End Workflow

## 1. Read Quality Control

The pipeline began with `FastQC` to assess read quality across paired-end FASTQ files.

```bash
#!/bin/bash
#SBATCH -J fastqc
#SBATCH -p cpunode
#SBATCH -o fastqc1.out
#SBATCH -e fastqc1.err
#SBATCH -N 1 -n 4

ml load fastqc-0.11.9-gcc-8.5.0-couiege

FASTQDIR=/data/ncbi_data/fastq

for fastq12 in \
  SRR12893917 SRR12893916 SRR12893913 SRR12893912 SRR12893911 \
  SRR12893910 SRR12893909 SRR12893908 SRR12893907 SRR12893906 \
  SRR12893905 SRR12893904 SRR12893903 SRR12893902 SRR12893901 \
  SRR12893900 SRR12893899 SRR12893898 SRR12893897 SRR12893896 \
  SRR12893895 SRR12893894 SRR12893884 SRR12893883 SRR12893881 \
  SRR12893880 SRR12893879 SRR12893859 SRR12893854 SRR12893852 \
  SRR12893849 SRR12893845 SRR12893844 SRR12893842 SRR12893841 \
  SRR12893834 SRR12893831 SRR12893830 SRR12893829 SRR12893826 \
  SRR12893824 SRR12893823 SRR12893822 SRR12893821 SRR12893820 \
  SRR12893818 SRR12893816 SRR12893815 SRR12893814 SRR12893812 \
  SRR12893811 SRR12893810 SRR12893804 SRR12893793 SRR12893782 \
  SRR12893771 SRR12893694 SRR12893972 SRR12893950
do
  fastqc "$FASTQDIR/${fastq12}_1.fastq" "$FASTQDIR/${fastq12}_2.fastq" \
    -o /data/bio/zengyg/Sample_Fastqc/
done
```

Quality summaries were then aggregated using `MultiQC`.

```bash
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

## 2. Decompression and Adapter Trimming

Some raw files were decompressed in parallel before downstream processing.

```bash
#!/bin/bash
#SBATCH -J gunzip
#SBATCH -p cpunode
#SBATCH -o gunzip-%j.out
#SBATCH -e gunzip-%j.err
#SBATCH -N 1 -n 4

ml load parallel-20210922-gcc-8.5.0-v7ilxoy

cat Samples.txt | parallel -j 4 gunzip /data/bio/jiangmj/X101SC22124970-Z01-J001/00.CleanData/{}/*.fq.gz
```

Adapter trimming was performed with `AdapterRemoval`.

```bash
#!/bin/bash
#SBATCH -J AdapterRemoval
#SBATCH -p cpunode
#SBATCH -o adapterremoval_sample.out
#SBATCH -e adapterremoval_sample.err
#SBATCH -N 1 -n 4

ml load AdapterRemoval-2.3.3

FASTQDIR=/data/ncbi_data/fastq

for fastq12 in \
  SRR12893917 SRR12893916 SRR12893913 SRR12893912 SRR12893911 \
  SRR12893910 SRR12893909 SRR12893908 SRR12893907 SRR12893906 \
  SRR12893905 SRR12893904 SRR12893903 SRR12893902 SRR12893901 \
  SRR12893900 SRR12893899 SRR12893898 SRR12893897 SRR12893896 \
  SRR12893895 SRR12893894 SRR12893884 SRR12893883 SRR12893881 \
  SRR12893880 SRR12893879 SRR12893859 SRR12893854 SRR12893852 \
  SRR12893849 SRR12893845 SRR12893844 SRR12893842 SRR12893841 \
  SRR12893834 SRR12893831 SRR12893830 SRR12893829 SRR12893826 \
  SRR12893824 SRR12893823 SRR12893822 SRR12893821 SRR12893820 \
  SRR12893818 SRR12893816 SRR12893815 SRR12893814 SRR12893812 \
  SRR12893811 SRR12893810 SRR12893804 SRR12893793 SRR12893782 \
  SRR12893771 SRR12893694 SRR12893972 SRR12893950
do
  AdapterRemoval \
    --file1 "$FASTQDIR/${fastq12}_1.fastq" \
    --file2 "$FASTQDIR/${fastq12}_2.fastq" \
    --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --output1 "/data/bio/zengyg/Sample_AdapterRemoval/${fastq12}_adapterremove_1.fastq" \
    --output2 "/data/bio/zengyg/Sample_AdapterRemoval/${fastq12}_adapterremove_2.fastq"
done
```

## 3. Host and Contaminant Removal with KneadData

After trimming, `KneadData` was used to clean reads and remove contamination.

```bash
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
```

For samples, I also used a second-stage filtering strategy against a control-derived Bowtie2 reference.

```bash
ml load kneaddata-0.12.0
bowtie2-build final_assembly.fasta Controls
```

```bash
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

cat samples.txt | parallel -j 8 \
  kneaddata --bypass-trf \
    -i1 /data/bio/tmp1sst1/kneaddata1103/{}/{}_adapterremove_1_kneaddata_paired_1.fastq \
    -i2 /data/bio/tmp1sst1/kneaddata1103/{}/{}_adapterremove_1_kneaddata_paired_2.fastq \
    -db /data/bio/tmp1sst1/Controls_reference/Controls \
    -t 2 -p 4 \
    --remove-intermediate-output \
    -o /data/bio/tmp1sst1/kneaddata1103/{} \
    --bowtie2 "$BOWTIE2PATH" \
    --trimmomatic "$TRIMMOMATIC"
```

## 4. Assembly and Genome Binning

After generating clean reads, I used `metaWRAP` to assemble reads and construct bins.

### Control Co-assembly

```bash
cat CLEAN_READS/ERR*_1.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/ERR*_2.fastq > CLEAN_READS/ALL_READS_2.fastq
```

```bash
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

metawrap assembly \
  -1 $METAWRAPDIR/ALL_READS_1.fastq \
  -2 $METAWRAPDIR/ALL_READS_2.fastq \
  -t 4 \
  -o /data/bio/zengyg/Co_assembly_control
```

### Sample-level Assembly

```bash
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

metawrap assembly \
  -1 /data/bio/zengyg/Samples_Clean/SRR12893811_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq \
  -2 /data/bio/zengyg/Samples_Clean/SRR12893811_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq \
  -t 8 \
  -o /data/bio/zengyg/SRR12893811_assembly
```

### FASTQ Repair and Renaming

```bash
ml load bbmap-39.01
ml load parallel-20210922-gcc-8.5.0-v7ilxoy

cat samples.txt | parallel -j 30 \
  repair.sh \
    in=/data/bio/zengyg/Samples_Clean/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq \
    in2=/data/bio/zengyg/Samples_Clean/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq \
    out=/data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq \
    out2=/data/bio/zengyg/Samples_Clean_repair/{}_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq
```

```bash
rename adapterremove_1_kneaddata_paired_1_kneaddata repaired *
```

### Binning

```bash
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
```

### Bin Refinement

```bash
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

cat samples-binning.txt | parallel -j 8 \
  metawrap bin_refinement \
    -o /data/bio/zengyg/Bin_Refinement_1SRR_plot_numpy/SRR12893880 \
    -t 2 \
    -A /data/bio/zengyg/Initial_Binning/SRR12893880_backup2/metabat2_bins \
    -B /data/bio/zengyg/Initial_Binning/SRR12893880/concoct_bins \
    -c 50 \
    -x 10
```

## 5. Taxonomic Profiling with Kraken2

Taxonomic profiling was performed with `Kraken2`.

```bash
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
```

A sample-specific test run was also used during setup and debugging.

```bash
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

kraken2 \
  --db /data/bio/zengyg/minikraken2_v1_8GB_Update \
  --report /data/bio/zengyg/kraken2_report/SRR12893972_report \
  --paired \
  /data/bio/zengyg/Samples_Clean_repair/SRR12893972_adapterremove_1_kneaddata_paired_1_kneaddata_paired_1.fastq \
  /data/bio/zengyg/Sample_Clean_repair/SRR12893972_adapterremove_1_kneaddata_paired_1_kneaddata_paired_2.fastq \
  --use-names \
  --report-zero-counts \
  --output /data/bio/zengyg/kraken2_output
```

## 6. Genome Quality and Taxonomic Annotation

To interpret bins and MAGs, I used `CheckM` and `GTDB-Tk`.

### CheckM Summary Conversion Helper

One small Python utility converted structured `CheckM` output into a flatter tabular form that was easier to inspect and export.

```python
#!/usr/bin/env python3

import json
import sys

HEADER = [
    "Bin Id",
    "Marker Lineage",
    "Genomes",
    "Markers",
    "Marker Sets",
    "0",
    "1",
    "2",
    "3",
    "4",
    "5+",
    "Completeness",
    "Contamination",
    "GC",
    "GC std",
    "Genome Size",
    "Ambiguous Bases",
    "Scaffolds",
    "Contigs",
    "Translation Table",
    "Predicted Genes",
]

def load_stats(path):
    stats = {}
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.replace("'", '"').rstrip("\n")
            bin_id, payload = line.split("\t", 1)
            normalized_id = "bin_" + bin_id.replace(".bin.", "_")
            stats[normalized_id] = json.loads(payload)
    return stats

def write_summary(stats, path):
    with open(path, "w", encoding="utf-8") as output:
        output.write("\t".join(HEADER) + "\n")
        for key, values in stats.items():
            row = [
                key,
                values["marker lineage"],
                str(values["# genomes"]),
                str(values["# markers"]),
                str(values["# marker sets"]),
                str(values["0"]),
                str(values["1"]),
                str(values["2"]),
                str(values["3"]),
                str(values["4"]),
                str(values["5+"]),
                str(values["Completeness"]),
                str(values["Contamination"]),
                str(values["GC"]),
                str(values["GC std"]),
                str(values["Genome size"]),
                str(values["# ambiguous bases"]),
                str(values["# scaffolds"]),
                str(values["# contigs"]),
                str(values["Translation table"]),
                str(values["# predicted genes"]),
            ]
            output.write("\t".join(row) + "\n")

stats = load_stats(sys.argv[1])
write_summary(stats, sys.argv[2])
```

### GTDB-Tk Classification

```bash
export LC_ALL="en_US.UTF-8"

ml load python-3.9.12-gcc-8.5.0-ews4x7x
ml load pplacer-1.1.alpha19-gcc-8.5.0-a6jmrxw
ml load hmmer-3.3.2-gcc-8.5.0-r4vrjwg
ml load fasttree-2.1.10-gcc-8.5.0-jzumr6k
export GTDBTK_DATA_PATH=/data/public_data/gtdbtk-data/release207_v2

gtdbtk classify_wf \
  --genome_dir /data/bio/zengyg/SRRBIN \
  --extension fa \
  --out_dir /data/bio/zengyg/GTDB_classify_wf \
  --cpus 1
```

## 7. Species Abundance Analysis and Visualization in R

Once taxonomic counts were generated, I used R for data wrangling and visualization.

### Relative Abundance Workflow

This script restructures species-level Kraken output, merges it with metadata, computes relative abundance, and generates a stacked abundance plot.

```r
library(devtools)
library(wilkoxmisc)
library(reshape2)
library(dplyr)
library(readr)
library(ggplot2)
library(starsExtra)

OTU <- read.csv("kraken_report_all_species.csv")

OTU <- OTU %>%
  gather(Sample, Count, 2:ncol(OTU)) %>%
  filter(Count > 0)

names(OTU)[1] <- "Species"
OTU$Sample <- gsub("kraken.report.", "", OTU$Sample)
OTU$Sample <- gsub("_kraken_report", "", OTU$Sample)

Meta <- read_tsv("metadata_doorknob_skin.txt")
OTU <- OTU %>% filter(OTU$Sample %in% Meta$Sample)

OTU <- ddply(OTU, .(Sample, Species), summarise, count = sum(Count))
OTU <- ddply(OTU, .(Sample), mutate, RelativeAbundance = (count * 100) / sum(count))

write_tsv(OTU, "kraken_species_abundance.tidy.txt")

OTU <- read_tsv("kraken_species_abundance.tidy.txt")
Meta <- read_tsv("metadata_doorknob_skin.txt")
Meta2 <- Meta %>% filter(Location == "Residence 3")
OTUTable <- merge(OTU, Meta2, by = "Sample", all.x = TRUE)
OTU3 <- OTUTable %>% filter(Location == "Residence 3")

Top12 <- collapse_taxon_table(OTU3, n = 12, Rank = "Species")
Top12species <- merge(Top12, Meta2, by = "Sample", all.x = TRUE)
write_tsv(Top12species, "Top12species.txt")

SpeciesAbundance <- read_tsv("Top12Species.txt")

Plot <- ggplot(
  SpeciesAbundance,
  aes(
    x = factor(
      Sampling_day,
      levels = c(
        "Day 1 day", "Day 1 night", "Day 2 day", "Day 2 night",
        "Day 3 day", "Day 3 night", "Day 4 day", "Day 4 night",
        "Day 5 day", "Day 5 night", "Day 6 day", "Day 6 night",
        "Day 7 day", "Day 7 night", "Day 8 day", "Day 8 night",
        "Day 9 day", "Day 9 night", "Day 10 day", "Day 10 night"
      )
    ),
    y = RelativeAbundance,
    fill = factor(Species)
  )
)

Plot <- Plot + geom_bar(stat = "identity")
Plot <- Plot + scale_fill_brewer(palette = "Paired")
Plot <- Plot + theme_classic()
Plot <- Plot + facet_grid(Sample_type ~ Location, scales = "free", space = "free")
Plot <- Plot + scale_y_continuous(expand = c(0, 0))
Plot <- Plot + ylab("Relative Abundance (%)") + xlab("Sampling day") + labs(fill = "Species")
Plot <- Plot + theme(axis.text.x = element_text(hjust = 0, angle = 90))

ggsave("taxonomy_general.png", width = 12, height = 7)
ggsave("top12species_by_sample_type.pdf", dpi = 1086)
```

### Bray-Curtis Dissimilarity

This stage transformed abundance tables into ecological distance matrices and then compared doorknob and palm communities across time and residence.

```r
library(reshape2)
library(wilkoxmisc)
library(ggplot2)
library(dplyr)
library(cultevo)
library(readr)
library(vegan)
library(MASS)
library(cowplot)

tax_cast <- read.table("kraken_species_abundance.cast.txt", header = TRUE, row.names = 1)
distance <- vegdist(tax_cast, method = "bray")
Matrix <- as.matrix(distance)
write.matrix(Matrix, sep = " ", "bray_curtis_distance_matrix.txt")

UniFrac <- read_tsv("bray_curtis_distance_matrix.txt")
UniFrac <- melt(UniFrac, value.name = "Distance")
names(UniFrac)[1:2] <- c("Sample1", "Sample2")

Samples <- read_tsv("metadata_doorknob_skin.txt") %>%
  select(Sample, Sample_type, Location, Sampling_day)

names(Samples)[1] <- "Sample2"
UniFrac <- merge(UniFrac, Samples, by = "Sample2", all.x = FALSE)
names(UniFrac)[4:6] <- c("Sample_type2", "Location2", "Sampling_day2")

names(Samples)[1] <- "Sample1"
UniFrac <- merge(UniFrac, Samples, by = "Sample1", all.x = FALSE)
names(UniFrac)[7:9] <- c("Sample_type1", "Location1", "Sampling_day1")

time <- read_tsv("sampling_time_hour_conversion.txt")
names(time)[1] <- "Sampling_day1"
merged <- left_join(UniFrac, time)
names(merged)[10] <- "Time_hour1"

names(time)[1] <- "Sampling_day2"
merged <- left_join(merged, time)
names(merged)[11] <- "Time_hour2"

UniFrac <- merged
UniFrac <- UniFrac[which(!UniFrac$Sample1 == UniFrac$Sample2), ]
UniFrac$HouseholdType <- ifelse(UniFrac$Location1 == UniFrac$Location2, "Same household", "Different household")

table <- UniFrac %>%
  filter(HouseholdType == "Same household") %>%
  filter(Sample_type2 == "Door knob") %>%
  filter(Sample_type1 %in% c("Right palm", "Left palm"))

table$Group <- paste0(table$Sample_type2, " vs. ", table$Sample_type1)
table <- table %>% mutate(time_decay = Time_hour1 - Time_hour2)

Plot1 <- ggplot(table, aes(x = time_decay, y = Distance))
Plot1 <- Plot1 + geom_point(color = "#455def", alpha = 0.4)
Plot1 <- Plot1 + geom_smooth(color = "#fd017d", fill = "#fd017d")
Plot1 <- Plot1 + theme_bw()
Plot1 <- Plot1 + geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5)
Plot1 <- Plot1 + ylab("Bray-Curtis dissimilarity between doorknob and palm")
Plot1 <- Plot1 + xlab("Time interval (hours)")
Plot1 <- Plot1 + theme(legend.position = "none")

ggsave(plot = Plot1, "Bray-Curtis_dissimilarity_between_doorknob_and_palm.png", width = 5, height = 5)
write_tsv(table, "doorkob_palm_comparison_Bray_Curtis.txt")
```

### Sloan Neutral Model

This portion modeled the relationship between mean species abundance on skin and observed occurrence frequency on the doorknob.

```r
library(readr)
library(Hmisc)
library(bbmle)
library(wilkoxmisc)
library(reshape2)

table <- read_tsv("kraken_report_all_table_rarefied.txt")
meta <- read_tsv("metadata_doorknob_skin.txt") %>% select(Sample, Type, Location)

merge <- left_join(table, meta)

sub_table <- merge %>%
  filter(Type == "Indoor surface") %>%
  filter(Location == "Residence 3") %>%
  select(-Type, -Sample, -Location)

write_tsv(sub_table, "kraken_report_all_table_rarefied_doorknob_residence3.txt")

pool <- read_tsv("kraken_report_all_table_rarefied_skin_residence3.txt")
spp <- read_tsv("kraken_report_all_table_rarefied_doorknob_residence3.txt")

sncm.fit <- function(spp, pool = TRUE, stats = FALSE, taxon = NULL) {
  require(minpack.lm)
  require(Hmisc)
  require(stats4)

  options(warn = -1)
  N <- mean(apply(spp, 1, sum))

  if (is.null(pool)) {
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m / N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m / N
  }

  spp.bi <- 1 * (spp > 0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]

  C <- merge(p, freq, by = 0)
  C <- C[order(C[, 2]), ]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]
  p <- C.0[, 2]
  freq <- C.0[, 3]
  names(p) <- C.0[, 1]
  names(freq) <- C.0[, 1]

  d <- 1 / N
  m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))
  freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
  pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = TRUE)

  A <- cbind(p, freq, freq.pred, pred.ci[, 2:3])
  A <- as.data.frame(A)
  colnames(A) <- c("p", "freq", "freq.pred", "pred.lwr", "pred.upr")
  A[order(A[, 1]), ]
}

stats.otu <- sncm.fit(spp, pool)
write.table(stats.otu, "skin_to_doorknob_residence3.txt", sep = "\t", col.names = NA)
```

The visualization step highlighted taxa above, below, or within neutral expectations.

```r
library(readr)
library(dplyr)
library(ggplot2)

skin <- read_tsv("skin_to_doorknob_residence3.txt")
names(skin)[1] <- "Species"
skin <- skin %>% mutate(Partition = ifelse(freq < pred.lwr, "Below", "Neutral"))
skin <- skin %>% mutate(Partition = ifelse(freq > pred.upr, "Above", Partition))

spp <- read_tsv("kraken_report_all_table_rarefied_skin_residence3.txt")
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m / N
p <- as.data.frame(p)
colnames(p) <- "mean_abundance"
write.table(p, "kraken_species_mean_abundance_skin_residence3.txt", sep = "\t", row.names = TRUE)

abundance <- read_tsv("kraken_species_mean_abundance_skin_residence3.txt")
names(abundance)[1] <- "Species"
merged <- left_join(skin, abundance)
write_tsv(merged, "skin_to_doorknob_residence3_merge.txt")

cbPalette <- c("#537bff", "#ff3c6d", "#b248ff")
Table <- read_tsv("skin_to_doorknob_residence3.txt") %>% filter(!Species == "unclassified")

Plot <- ggplot(Table, aes(x = log(mean_abundance)))
Plot <- Plot + geom_point(aes(y = freq, colour = Partition), size = 1.5, alpha = 0.8)
Plot <- Plot + geom_line(aes(y = freq.pred))
Plot <- Plot + geom_ribbon(aes(ymin = pred.lwr, ymax = pred.upr), fill = "#969696", alpha = 0.5)
Plot <- Plot + scale_color_manual(values = cbPalette)
Plot <- Plot + theme_classic()
Plot <- Plot + theme(legend.title = element_blank())
Plot <- Plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
Plot <- Plot + ggtitle("Residence 3")
Plot <- Plot + xlab("log (mean relative abundance of species on skin)")
Plot <- Plot + ylab("Occurrence frequency of species on doorknob")

ggsave("sloan_neutral_model_skin_to_doorknob_residence3.png", width = 12, height = 7)
```

## 8. Dereplication and Genome Tracking

### dRep

To reduce redundancy among genomes, I used `dRep`.

```bash
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
```

I then used R to visualize temporal/sample-type patterns for dereplicated genomes.

```r
library(ggplot2)
library(dplyr)

dRep_rMAGs <- read.csv("dRep_rMAGs.csv")
dRep_rMAGs_filter <- dRep_rMAGs %>% filter(Genotype == "Cluster A")

plot <- ggplot(data = dRep_rMAGs_filter)
plot <- plot + geom_point(
  mapping = aes(
    x = factor(
      Sampling_day,
      levels = c(
        "Day 1 day", "Day 1 night", "Day 2 day", "Day 2 night",
        "Day 3 day", "Day 3 night", "Day 4 day", "Day 4 night",
        "Day 5 day", "Day 5 night", "Day 6 day", "Day 6 night",
        "Day 7 day", "Day 7 night", "Day 8 day", "Day 8 night",
        "Day 9 day", "Day 9 night", "Day 10 day", "Day 10 night"
      )
    ),
    y = Sample_type,
    colour = factor(Sample_type),
    size = 2
  )
)

plot <- plot + labs(x = "Sampling day", y = "Sample type", colour = "Sample type")
cbPalette <- c("#7998ff", "#ff6780", "#c371ff")
plot <- plot + scale_color_manual(values = cbPalette)
plot <- plot + facet_wrap(~classification)
plot <- plot + theme(axis.text.x = element_text(hjust = 0, angle = 90))
plot <- plot + guides(size = "none")

ggsave("taxonomy_day.png", width = 12, height = 7)
```

### CoverM

`CoverM` was used to estimate abundance and track a selected genome across samples.

```bash
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
```

The result was visualized in R as follows:

```r
library(wilkoxmisc)
library(reshape2)
library(dplyr)
library(readr)
library(ggplot2)
library(starsExtra)
library(tidyverse)
library(readxl)

Meta <- metadata_doorknob_skin
Meta2 <- Meta %>% filter(Location == "Residence 3")
Match <- match
MergeData <- merge(Meta2, Match, by = "Sample", all.x = TRUE)
write_tsv(MergeData, "Information.txt")

coverm_output <- coverm_results

plot <- ggplot(data = coverm_output)
plot <- plot + geom_point(
  mapping = aes(
    x = factor(
      Sampling_day,
      levels = c(
        "Day 1 day", "Day 1 night", "Day 2 day", "Day 2 night",
        "Day 3 day", "Day 3 night", "Day 4 day", "Day 4 night",
        "Day 5 day", "Day 5 night", "Day 6 day", "Day 6 night",
        "Day 7 day", "Day 7 night", "Day 8 day", "Day 8 night",
        "Day 9 day", "Day 9 night", "Day 10 day", "Day 10 night"
      )
    ),
    y = Sample_type,
    colour = Sample_type,
    size = RelativeAbundance
  )
)

plot <- plot + labs(
  x = "Sampling day",
  y = "Sample type",
  colour = "Sample type",
  size = "RelativeAbundance(%)"
)
plot <- plot + facet_wrap(~classification)
plot <- plot + theme(axis.text.x = element_text(hjust = 0, angle = 90))

ggsave("coverM_transmission.png", width = 12, height = 7)
```

# Outputs and Deliverables

The main outputs from this project include:

- quality control summaries
- cleaned and filtered reads
- assemblies and bins
- taxonomic abundance tables
- distance matrices and ecological comparisons
- neutral model outputs
- dereplication summaries
- abundance-tracking plots
- final thesis writing and presentation materials

# Portfolio Value of This Project

This project reflects several types of technical work:

- pipeline orchestration on HPC systems
- large-scale file handling and tool chaining
- reproducible batch job scripting
- data wrangling across shell, Python, and R
- applied ecological and microbiome-style statistical analysis
- scientific figure generation for presentation and thesis writing

# Notes

- Paths and module names reflect the original compute environment.
- Raw sequencing data are not included in this portfolio repository.
- The goal of this document is to tell the full technical story in one place, so some code has been cleaned for readability while preserving the original logic and workflow.

