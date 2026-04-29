# Pipeline Run Notes

Detailed notes on parameter choices, runtime observations, and gotchas
encountered during the 2023 analysis.

---

## Step 1 — Quality Control (FastQC + Trimmomatic)

**Tool versions**: FastQC 0.11.9, Trimmomatic 0.39  
**Cluster**: 4 CPUs, 8 GB RAM, ~20 min per sample

### Adapter choice
All libraries were prepared with TruSeq Nano DNA HT (PE) kits.  
Adapter file: `TruSeq3-PE-2.fa` (ships with Trimmomatic).

### Trimming parameters
```
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true
LEADING:3  TRAILING:3  SLIDINGWINDOW:4:20  MINLEN:50
```
- `SLIDINGWINDOW:4:20` removes runs of quality < Q20.
- `MINLEN:50` was chosen over the default 36 because very short reads
  map ambiguously to the Kraken2 database.

### Pre/post-trim read counts
Average read retention after trimming: **94.8%** (range 91–97%).
All samples passed with > 3 M read pairs remaining.

---

## Step 2 — Host Decontamination (Bowtie2)

**Tool version**: Bowtie2 2.5.1  
**Host genome**: *Zea mays* B73 RefGen v5 (GCF_902167145)  
**Index location on HPC**: `/scratch/$USER/databases/host_genome/`

### Why Bowtie2 `--sensitive` and not `--very-sensitive`?
Pilot tests on 2 samples showed that `--very-sensitive` added ~45 min
per sample with < 0.01% difference in removed reads.  `--sensitive` was
used for all 24 samples.

### Host mapping rates
Average: **2.3%** of read pairs mapped to the host genome and were removed.
One agricultural sample (AG-07) had a higher host rate (8.1%), consistent
with higher plant root material in that core.

---

## Step 3 — Taxonomic Profiling (Kraken2 + Bracken)

**Tool versions**: Kraken2 2.1.3, Bracken 2.8  
**Database**: Standard Kraken2 DB built 2023-03-14 (bacteria, archaea, viral,
human, plasmid, UniVec_Core)  
**DB size on disk**: ~66 GB (requires highmem partition, 128 GB RAM allocated)

### Confidence threshold
`--confidence 0.1` was used (10% k-mer agreement required for classification).
Lower values (0.0) increased classification rate slightly but introduced
more spurious assignments at species level; 0.1 is the literature-recommended
default for environmental samples.

### Bracken settings
- Read length: 150 bp (matches actual read length post-trimming)
- Threshold: 10 (minimum reads at rank to re-estimate)
- Species-level re-estimation used for ecological analysis.
- Genus-level re-estimation used for DESeq2 differential abundance (more
  statistical power due to aggregation).

### Classification rates
Average: **68%** of reads classified (range 52–79%).
Lower rates in forest samples are consistent with higher novelty / diversity.

---

## Step 4 — De Novo Assembly (MEGAHIT)

**Tool version**: MEGAHIT 1.2.9  
**Cluster**: 24 CPUs, 128 GB RAM, 4–18 h per sample

### Per-sample vs. co-assembly
Per-sample assembly was chosen because:
1. Samples span a strong environmental gradient — co-assembly can chimera
   contigs from ecologically distinct communities.
2. Per-sample assembly preserves strain-level variation needed for accurate
   binning.

### k-mer list
`--k-list 21,29,39,59,79,99,119,141`  
Odd k-mers avoid palindromic sequences; large range handles both
low-coverage novel taxa and high-coverage dominant taxa.

### N50 summary (representative)
| Sample  | Contigs (≥1 kb) | N50 (bp) | Largest (bp) |
|---------|----------------|----------|-------------|
| FOR-01  | 42 318         | 2 841    | 187 209     |
| FOR-06  | 38 710         | 3 102    | 214 551     |
| AG-01   | 28 504         | 2 218    | 98 430      |
| AG-12   | 31 847         | 2 409    | 112 007     |

---

## Step 5 — Coverage Mapping (Bowtie2)

**Tool version**: Bowtie2 2.5.1, SAMtools 1.17  
**Cluster**: 16 CPUs, 32 GB RAM, ~2 h per sample

Mapped trimmed reads back to per-sample contigs to generate depth profiles
for MetaBAT2. Average mapping rate back to own assembly: **78%**.

---

## Step 6 — Genome Binning (MetaBAT2)

**Tool version**: MetaBAT2 2.15  
**Min contig length**: 2 500 bp  
**Cluster**: 16 CPUs, 64 GB RAM, ~1 h per sample

MetaBAT2 uses tetranucleotide frequency + coverage depth jointly.
`--unbinned` flag retains unbinned contigs as a separate FASTA for
subsequent CONCOCT comparison (not shown in final analysis).

---

## Step 7 — Bin Quality (CheckM2)

**Tool version**: CheckM2 1.0.1  
**DB**: `uniref100.KO.1.dmnd` (DIAMOND database, ~3 GB)  
**Cluster**: 16 CPUs, 32 GB RAM, ~30 min per sample

CheckM2 uses a trained neural network on DIAMOND alignments to reference
genomes. More accurate than CheckM1 marker-gene approach for novel lineages.

### Final MAG tally
| Tier            | Criteria                   | Count |
|-----------------|----------------------------|-------|
| High quality    | ≥ 90% comp, ≤ 5% cont      | 12    |
| Medium quality  | ≥ 70% comp, ≤ 10% cont     | 35    |
| Low quality     | passing bins below MQ      | 61    |
| **Total MAGs**  |                            | **108** |

---

## Step 8 — Phylogenomic Classification (GTDB-Tk)

**Tool version**: GTDB-Tk 2.3.2, **Reference DB**: GTDB r214  
**Cluster**: 48 CPUs, 256 GB RAM, ~36 h (single run on all passing MAGs)

`--skip_ani_screen` was used because the pplacer tree is more informative
for the novel Acidobacteria bins; ANI screening can prematurely classify them
at genus level to the nearest reference.

### Notable taxonomic findings
- 3 MAGs placed in unnamed families (f__) within *Acidobacteriota*,
  phyla p__Acidobacteriota; average ANI to nearest named genus < 85%.
- 8 MAGs in *Verrucomicrobiota*, consistent with soil-typical taxa.
- 2 MAGs in candidate phylum *FCPU426* — entirely novel to this soil type.

---

## Software Modules on Cluster

All tools were loaded via the university's `module` system (Lmod).
Examples:
```bash
module load fastqc/0.11.9
module load trimmomatic/0.39
module load bowtie2/2.5.1
module load samtools/1.17
module load kraken2/2.1.3
module load bracken/2.8
module load megahit/1.2.9
module load metabat2/2.15
module load checkm2/1.0.1
module load gtdbtk/2.3.2
```

Conda was used for R packages and Python utilities (see `docs/environment.yml`).
