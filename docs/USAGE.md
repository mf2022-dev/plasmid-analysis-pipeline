# Usage Guide

## Basic Usage

### 1. Prepare your samplesheet

Create a CSV file listing your samples:

```csv
sample,short_reads_1,short_reads_2,long_reads,assembly
KP001,/data/KP001_R1.fastq.gz,/data/KP001_R2.fastq.gz,/data/KP001_ont.fastq.gz,
KP002,/data/KP002_R1.fastq.gz,/data/KP002_R2.fastq.gz,,
KP003,,,,/data/KP003_assembly.fasta
```

### 2. Run the pipeline

**Hybrid assembly mode (recommended):**
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --assembly_mode hybrid \
    --bakta_db /path/to/bakta_db/ \
    -profile docker
```

**Pre-assembled genomes (skip assembly):**
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --skip_assembly \
    --bakta_db /path/to/bakta_db/ \
    -profile docker
```

**With reference plasmid for BRIG comparison:**
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --reference /path/to/reference_plasmid.fasta \
    --bakta_db /path/to/bakta_db/ \
    -profile docker
```

### 3. Skip specific analyses

```bash
# Skip MGE analysis (faster)
nextflow run main.nf --input samplesheet.csv --skip_mge -profile docker

# Skip visualization (faster)
nextflow run main.nf --input samplesheet.csv --skip_viz -profile docker
```

### 4. Resume failed runs

```bash
nextflow run main.nf --input samplesheet.csv -profile docker -resume
```

## Advanced Usage

### Custom resource allocation

Edit `nextflow.config` or use command-line overrides:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --max_cpus 64 \
    --max_memory '256.GB' \
    --max_time '96.h' \
    -profile slurm,singularity
```

### Custom PlasmidFinder thresholds

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --plasmidfinder_cov 0.80 \
    --plasmidfinder_id 0.98 \
    -profile docker
```

## Interpreting Results

See the companion SOP document (`docs/SOP_Plasmid_Genomic_Analysis.md`) for detailed guidance on interpreting all pipeline outputs.
