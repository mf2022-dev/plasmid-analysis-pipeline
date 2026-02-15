# PlasmidScope

**Comprehensive Bacterial Plasmid Genomic Analysis Pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.xxxx%2Fxxxxx-blue.svg)]()

---

## Overview

**PlasmidScope** is a state-of-the-art Nextflow pipeline for comprehensive bacterial plasmid characterization, mobile genetic element (MGE) analysis, and genomic environment visualization. It is designed for researchers who need publication-quality plasmid analysis suitable for high-impact journals such as *The Lancet Microbe*, *Nature Microbiology*, and *Antimicrobial Agents and Chemotherapy*.

The pipeline integrates over 20 bioinformatics tools into a single, reproducible workflow that takes raw sequencing data (or pre-assembled genomes) and produces:

- Complete plasmid reconstruction and enumeration
- Replicon typing (Inc groups) and mobility classification
- AMR and virulence gene detection with plasmid-level mapping
- Insertion sequence, transposon, and integron identification
- Circular (BRIG) and linear (Clinker) comparative visualizations
- Plasmid clustering and transmission network analysis
- Publication-quality HTML report with figures and summary tables

---

## Pipeline Architecture

```
                    ┌─────────────────────────────────┐
                    │       INPUT SAMPLESHEET          │
                    │  (FASTQ / FASTA assemblies)      │
                    └──────────────┬──────────────────┘
                                   │
                    ┌──────────────▼──────────────────┐
                    │     PHASE 1: QUALITY CONTROL     │
                    │  FastQC │ NanoPlot │ Fastp │      │
                    │  Filtlong                         │
                    └──────────────┬──────────────────┘
                                   │
                    ┌──────────────▼──────────────────┐
                    │     PHASE 2: ASSEMBLY            │
                    │  Unicycler (hybrid)               │
                    │  Flye (long-read)                 │
                    │  Hybracter (automated)            │
                    │  QUAST + Bandage QC               │
                    └──────────────┬──────────────────┘
                                   │
                    ┌──────────────▼──────────────────┐
                    │     PHASE 3: ANNOTATION          │
                    │  Bakta (standardized)             │
                    └──────────────┬──────────────────┘
                                   │
              ┌────────────────────┼────────────────────┐
              │                    │                     │
   ┌──────────▼─────────┐ ┌───────▼──────────┐ ┌───────▼──────────┐
   │  PHASE 4: AMR/VIR  │ │  PHASE 5: TYPING │ │  PHASE 6: MGE    │
   │  ABRicate (CARD,   │ │  PlasmidFinder   │ │  ISfinder BLAST  │
   │  VFDB)             │ │  MOB-suite       │ │  IntegronFinder  │
   │  AMRFinderPlus     │ │  (replicon, MOB, │ │  MobileElement   │
   │                    │ │   mobility)      │ │  Finder          │
   └──────────┬─────────┘ └───────┬──────────┘ └───────┬──────────┘
              │                    │                     │
              └────────────────────┼─────────────────────┘
                                   │
                    ┌──────────────▼──────────────────┐
                    │  PHASE 7: COMPARATIVE ANALYSIS   │
                    │  BRIG (circular comparison)       │
                    │  Clinker (linear synteny)         │
                    │  Mash (distance + clustering)     │
                    └──────────────┬──────────────────┘
                                   │
                    ┌──────────────▼──────────────────┐
                    │     PHASE 8: REPORTING           │
                    │  HTML Report + Summary Tables     │
                    │  Publication-quality Figures      │
                    │  MultiQC                          │
                    └──────────────────────────────────┘
```

---

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) >= 23.04
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) or [Conda](https://docs.conda.io/)

### Installation

```bash
# Clone the repository
git clone https://github.com/alghoribima/plasmid-analysis-pipeline.git
cd plasmid-analysis-pipeline

# Verify Nextflow installation
nextflow -version
```

### Run with test data

```bash
nextflow run main.nf -profile test,docker
```

### Run with your data

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --assembly_mode hybrid \
    --bakta_db /path/to/bakta_db/ \
    -profile docker
```

---

## Input Samplesheet

Create a CSV file with the following columns:

```csv
sample,short_reads_1,short_reads_2,long_reads,assembly
isolate_001,/data/isolate_001_R1.fastq.gz,/data/isolate_001_R2.fastq.gz,/data/isolate_001_nanopore.fastq.gz,
isolate_002,/data/isolate_002_R1.fastq.gz,/data/isolate_002_R2.fastq.gz,,
isolate_003,,,,/data/isolate_003_assembly.fasta
```

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `short_reads_1` | No | Path to Illumina R1 FASTQ |
| `short_reads_2` | No | Path to Illumina R2 FASTQ |
| `long_reads` | No | Path to ONT/PacBio FASTQ |
| `assembly` | No | Path to pre-assembled FASTA |

Each sample must have at least one of: short reads, long reads, or a pre-assembled genome.

---

## Parameters

### Core Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | *required* | Path to input samplesheet CSV |
| `--outdir` | `./results` | Output directory |
| `--assembly_mode` | `hybrid` | Assembly mode: `hybrid`, `long`, `short`, `auto` |
| `--skip_assembly` | `false` | Skip assembly (use pre-assembled genomes) |

### Annotation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bakta_db` | *required* | Path to Bakta annotation database |
| `--genus` | `Klebsiella` | Genus for annotation |
| `--species` | `pneumoniae` | Species for annotation |

### Plasmid Typing Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--plasmidfinder_db` | `enterobacteriaceae` | PlasmidFinder database |
| `--plasmidfinder_cov` | `0.60` | Minimum coverage threshold |
| `--plasmidfinder_id` | `0.95` | Minimum identity threshold |

### MGE Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_mge` | `false` | Skip MGE analysis |
| `--isfinder_identity` | `90` | Minimum IS element identity (%) |
| `--isfinder_evalue` | `1e-10` | E-value threshold for IS detection |

### Visualization Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_viz` | `false` | Skip visualization steps |
| `--reference` | `null` | Reference plasmid for BRIG comparison |
| `--mash_threshold` | `0.04` | Mash distance threshold for clustering |

### Resource Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | `32` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `72.h` | Maximum time per process |

---

## Output Structure

```
results/
├── 01_qc/
│   ├── fastqc/              # Short-read quality reports
│   ├── nanoplot/             # Long-read quality reports
│   ├── fastp/                # Trimming reports
│   └── filtlong/             # Filtered long reads
├── 02_assembly/
│   ├── unicycler/            # Hybrid assemblies
│   └── quast/                # Assembly quality metrics
├── 03_annotation/
│   └── bakta/                # Genome annotations (GBK, GFF3, FAA)
├── 04_amr_virulence/
│   ├── abricate/             # AMR/virulence gene screening
│   └── amrfinderplus/        # NCBI AMR detection
├── 05_plasmid_typing/
│   ├── plasmidfinder/        # Replicon typing results
│   └── mob_suite/            # MOB typing, mobility prediction
├── 06_mge_analysis/
│   ├── isfinder/             # Insertion sequence detection
│   ├── integron_finder/      # Integron identification
│   └── mobile_element_finder/ # Comprehensive MGE detection
├── 07_comparative/
│   ├── brig/                 # Circular comparison figures
│   ├── clinker/              # Linear synteny figures
│   └── mash/                 # Distance matrix and clustering
├── 08_report/
│   ├── plasmid_analysis_report.html  # Main HTML report
│   ├── multiqc_report.html           # MultiQC report
│   ├── summary_tables/               # TSV summary tables
│   └── figures/                      # Publication-quality figures
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    └── pipeline_dag.svg
```

---

## Execution Profiles

| Profile | Description | Command |
|---------|-------------|---------|
| `docker` | Run with Docker containers | `-profile docker` |
| `singularity` | Run with Singularity containers | `-profile singularity` |
| `conda` | Run with Conda environments | `-profile conda` |
| `slurm` | Submit to SLURM HPC cluster | `-profile slurm,docker` |
| `pbs` | Submit to PBS HPC cluster | `-profile pbs,singularity` |
| `sge` | Submit to SGE HPC cluster | `-profile sge,docker` |
| `test` | Run with minimal test data | `-profile test,docker` |

---

## Tools and Citations

| Tool | Version | Purpose | Citation |
|------|---------|---------|----------|
| FastQC | 0.12.1 | Short-read QC | Andrews (2010) |
| NanoPlot | 1.42.0 | Long-read QC | De Coster et al. (2018) |
| Fastp | 0.23.4 | Read trimming | Chen et al. (2018) |
| Filtlong | 0.2.1 | Long-read filtering | Wick (2021) |
| Unicycler | 0.5.0 | Hybrid assembly | Wick et al. (2017) |
| Flye | 2.9.3 | Long-read assembly | Kolmogorov et al. (2019) |
| Hybracter | 0.7.0 | Automated hybrid assembly | Bouras et al. (2024) |
| QUAST | 5.2.0 | Assembly QC | Gurevich et al. (2013) |
| Bakta | 1.9.3 | Genome annotation | Schwengers et al. (2021) |
| ABRicate | 1.0.1 | AMR/virulence screening | Seemann (2024) |
| AMRFinderPlus | 3.12 | NCBI AMR detection | Feldgarden et al. (2021) |
| PlasmidFinder | 2.1 | Replicon typing | Carattoli et al. (2014) |
| MOB-suite | 3.1.0 | Plasmid typing/reconstruction | Robertson & Nash (2018) |
| ISfinder | DB | IS element detection | Siguier et al. (2006) |
| IntegronFinder | 2.0 | Integron detection | Neron et al. (2022) |
| MobileElementFinder | 1.1 | MGE detection | Johansson et al. (2021) |
| BRIG | 0.95 | Circular comparison | Alikhan et al. (2011) |
| Clinker | 0.0.28 | Linear synteny | Gilchrist & Chooi (2021) |
| Mash | 2.3 | Distance estimation | Ondov et al. (2016) |
| MultiQC | 1.21 | Report aggregation | Ewels et al. (2016) |

---

## Companion SOP

This pipeline is accompanied by a comprehensive **Standard Operating Procedure (SOP)** document:

> **SOP_Plasmid_Genomic_Analysis.md** - A state-of-the-art protocol for plasmid characterization, mobile genetic element analysis, and genomic environment visualization. This document provides detailed step-by-step instructions, decision trees, quality control checkpoints, and publication standards.

The SOP is located in the `docs/` directory and serves as the definitive reference for interpreting pipeline outputs and designing downstream analyses.

---

## Contributing

Contributions are welcome. Please submit issues and pull requests on GitHub.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Commit changes (`git commit -am 'Add new analysis module'`)
4. Push to branch (`git push origin feature/new-analysis`)
5. Open a Pull Request

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Acknowledgments

This pipeline was developed as part of the ST11 *Klebsiella pneumoniae* genomic epidemiology project. We acknowledge the developers of all integrated tools and databases.

---

## Contact

- **Author:** Majed Alghoribi
- **Email:** alghoribima@gmail.com
- **Institution:** King Abdulaziz Medical City, Riyadh, Saudi Arabia
