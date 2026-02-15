# Installation Guide

## Prerequisites

### 1. Nextflow

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verify
nextflow -version
```

### 2. Container Engine (choose one)

**Docker (recommended for local use):**
```bash
sudo apt-get update && sudo apt-get install -y docker.io
sudo usermod -aG docker $USER
```

**Singularity (recommended for HPC):**
```bash
# Follow instructions at https://sylabs.io/guides/latest/user-guide/quick_start.html
```

**Conda (alternative):**
```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create environment
conda env create -f envs/plasmidscope.yml
conda activate plasmidscope
```

### 3. Databases

**Bakta database (required for annotation):**
```bash
bakta_db download --output /path/to/bakta_db/ --type full
```

**ISfinder database (required for IS detection):**
```bash
# Download from https://isfinder.biotoul.fr/
# Set environment variable:
export ISFINDER_DB=/path/to/isfinder_db
```

**AMRFinderPlus database:**
```bash
amrfinder_update --force_update
```

## Quick Installation

```bash
git clone https://github.com/alghoribima/plasmid-analysis-pipeline.git
cd plasmid-analysis-pipeline
nextflow run main.nf -profile test,docker
```

## HPC Installation

For SLURM clusters:
```bash
git clone https://github.com/alghoribima/plasmid-analysis-pipeline.git
cd plasmid-analysis-pipeline
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --bakta_db /path/to/bakta_db/ \
    -profile slurm,singularity
```

For PBS clusters:
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --bakta_db /path/to/bakta_db/ \
    -profile pbs,singularity
```
