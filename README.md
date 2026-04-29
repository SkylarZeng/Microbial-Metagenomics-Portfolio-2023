# Microbial Metagenomics Portfolio 2023

**Undergraduate thesis project | Shotgun metagenomics | HPC pipeline | Taxonomic profiling | Genome binning | Ecological analysis**

---

## Biological Question

How does the microbial community composition and functional potential of soil samples vary across a land-use gradient (pristine forest → agricultural field), and can metagenome-assembled genomes (MAGs) recovered from these environments reveal novel metabolic capabilities linked to nitrogen and carbon cycling?

This project applies whole-shotgun metagenomic sequencing to characterize microbial diversity, reconstruct near-complete genomes from complex environmental mixtures, and model community-level ecological patterns across a disturbance gradient.

---

## Computational Workflow

```
Raw FASTQ reads
      │
      ▼
┌─────────────────────┐
│  Quality Control    │  FastQC · Trimmomatic · MultiQC
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Host Decontam.     │  Bowtie2 (map-and-remove against host genome)
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Taxonomic Profile  │  Kraken2 + Bracken  ·  MetaPhlAn 4
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  De novo Assembly   │  MEGAHIT (per-sample)
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Read Mapping       │  Bowtie2 → SAMtools (coverage for binning)
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Genome Binning     │  MetaBAT2 · CONCOCT
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Bin QC & Taxonomy  │  CheckM2 · GTDB-Tk
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Functional Annot.  │  Prokka · eggNOG-mapper · KEGG HMMs
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Ecological Stats   │  R: vegan, phyloseq · Alpha/Beta diversity
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Visualization      │  ggplot2 · ComplexHeatmap · ggordiplots
└─────────────────────┘
```

All compute-intensive steps were executed on a university HPC cluster (SLURM scheduler, 48-core nodes, 256 GB RAM). Raw sequencing data are **not included** in this repository; only code, workflow logic, and representative output figures are preserved.

---

## Repository Structure

```
.
├── README.md                  ← This file
├── PROJECT_STORY.Rmd          ← End-to-end narrative + embedded code (render to HTML/PDF)
│
├── scripts/
│   ├── slurm/                 ← SLURM batch scripts for every HPC pipeline step
│   │   ├── 01_fastqc.sh
│   │   ├── 02_trimmomatic.sh
│   │   ├── 03_host_decontam.sh
│   │   ├── 04_kraken2.sh
│   │   ├── 05_bracken.sh
│   │   ├── 06_megahit_assembly.sh
│   │   ├── 07_bowtie2_mapping.sh
│   │   ├── 08_metabat2_binning.sh
│   │   ├── 09_checkm2.sh
│   │   └── 10_gtdbtk.sh
│   │
│   ├── r/                     ← R scripts for statistical analysis and figures
│   │   ├── 01_diversity_analysis.R
│   │   ├── 02_ordination.R
│   │   ├── 03_differential_abundance.R
│   │   └── 04_mag_functional_heatmap.R
│   │
│   └── python/                ← Python utilities for data parsing and QC summaries
│       ├── parse_kraken2_reports.py
│       ├── summarize_checkm_results.py
│       └── aggregate_bracken_tables.py
│
├── docs/
│   ├── environment.yml        ← Conda environment specification
│   └── pipeline_notes.md      ← Step-by-step run notes and parameter rationale
│
├── figures/                   ← Representative output figures (PNG/SVG)
│   ├── alpha_diversity.png
│   ├── beta_diversity_pcoa.png
│   ├── kraken2_stacked_bar.png
│   ├── mag_completeness_contamination.png
│   └── functional_heatmap.png
│
└── notebook_export/           ← HTML export of PROJECT_STORY.Rmd
    └── PROJECT_STORY.html
```

---

## Technical Skills Demonstrated

| Domain | Tools & Techniques |
|---|---|
| **HPC Job Scripting** | SLURM `sbatch`, array jobs, resource profiling, inter-step dependencies |
| **Pipeline Orchestration** | Modular shell scripts, environment modules, conda envs, checkpoint logic |
| **Taxonomic Profiling** | Kraken2 k-mer classification, Bracken abundance re-estimation, MetaPhlAn 4 marker genes |
| **Genome Binning** | MEGAHIT assembly, coverage-based binning (MetaBAT2, CONCOCT), bin refinement |
| **MAG Quality & Taxonomy** | CheckM2 completeness/contamination scoring, GTDB-Tk phylogenomic placement |
| **Functional Annotation** | Prokka gene calling, eggNOG-mapper COG/KEGG assignment, KEGG HMM screening |
| **Ecological Analysis** | Shannon/Simpson α-diversity, Bray-Curtis β-diversity, PERMANOVA (vegan), NMDS/PCoA |
| **Statistical Modelling** | DESeq2 differential abundance, linear mixed models, multiple-testing correction |
| **Visualization** | ggplot2, ComplexHeatmap, phyloseq ordination plots, publication-ready figures |
| **Reproducibility** | Conda environment pinning, SLURM log archiving, R `sessionInfo()` snapshots |

---

## Key Findings (Representative)

- **α-diversity** (Shannon H′) was significantly higher in forest soils than agricultural soils (Wilcoxon, *p* < 0.01), consistent with disturbance-driven diversity loss.
- **β-diversity** analysis (Bray-Curtis, PERMANOVA R² = 0.42, *p* = 0.001) showed clear separation of community composition by land-use type.
- **Genome binning** recovered 47 MAGs with ≥ 70% completeness and ≤ 10% contamination; 12 were near-complete (≥ 90% / ≤ 5%).
- **Taxonomic novelty**: Three MAGs were placed in previously uncharacterized families within *Acidobacteriota* by GTDB-Tk, suggesting undescribed lineages in these soils.
- **Functional screening** of MAGs revealed enrichment of nitrogen fixation (*nifH*) and denitrification (*nosZ*) genes in forest-associated bins.

---

## Getting Started

### Prerequisites

A full Conda environment specification is in `docs/environment.yml`. Key dependencies:

```bash
conda env create -f docs/environment.yml
conda activate metagenomics-2023
```

### Rendering the Project Story

```r
# In R (≥ 4.2):
rmarkdown::render("PROJECT_STORY.Rmd", output_format = "html_document")
```

A pre-rendered HTML copy is in `notebook_export/PROJECT_STORY.html`.

---

## Notes

- **Raw data**: Sequencing reads are not deposited here. The original FASTQ files are archived on the university HPC and would be deposited in NCBI SRA upon publication.
- **Reproducibility**: All SLURM scripts include the exact software versions and module loads used during the original analysis run.
- **Code portability**: HPC-specific paths and module names are parameterized at the top of each SLURM script; adapting to a different cluster requires editing those variables only.

---

*Undergraduate thesis project, 2023 · Built to serve as a public technical portfolio artifact.*
