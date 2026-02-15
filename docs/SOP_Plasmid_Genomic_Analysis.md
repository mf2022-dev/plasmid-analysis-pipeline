# Standard Operating Procedure: Comprehensive Plasmid Genomic Analysis

## A State-of-the-Art Protocol for Plasmid Characterization, Mobile Genetic Element Analysis, and Genomic Environment Visualization

**Document ID:** MANUS-SOP-PLASMID-2026-01
**Version:** 2.0
**Effective Date:** 2026-02-15
**Classification:** Open Access Reference Protocol
**Author:** Manus AI

---

## Table of Contents

1. [Purpose and Rationale](#10-purpose-and-rationale)
2. [Scope](#20-scope)
3. [Definitions and Terminology](#30-definitions-and-terminology)
4. [Materials: Software, Databases, and Hardware Requirements](#40-materials-software-databases-and-hardware-requirements)
5. [Phase 1: Sequencing Strategy and Data Quality Control](#50-phase-1-sequencing-strategy-and-data-quality-control)
6. [Phase 2: Gold-Standard Genome Assembly](#60-phase-2-gold-standard-genome-assembly)
7. [Phase 3: Complete Plasmid Enumeration and Validation](#70-phase-3-complete-plasmid-enumeration-and-validation)
8. [Phase 4: Annotation and Foundational Characterization](#80-phase-4-annotation-and-foundational-characterization)
9. [Phase 5: Replicon Typing and Mobility Classification](#90-phase-5-replicon-typing-and-mobility-classification)
10. [Phase 6: Comparative Plasmid Genomics](#100-phase-6-comparative-plasmid-genomics)
11. [Phase 7: Mobile Genetic Element Deep-Dive](#110-phase-7-mobile-genetic-element-deep-dive)
12. [Phase 8: Genomic Environment and Synteny Analysis](#120-phase-8-genomic-environment-and-synteny-analysis)
13. [Phase 9: Plasmid Epidemiology and Transmission Analysis](#130-phase-9-plasmid-epidemiology-and-transmission-analysis)
14. [Phase 10: Interpretation, Reporting, and Publication Standards](#140-phase-10-interpretation-reporting-and-publication-standards)
15. [Quality Assurance and Control](#150-quality-assurance-and-control)
16. [Decision Trees](#160-decision-trees)
17. [Data Management and Reproducibility](#170-data-management-and-reproducibility)
18. [References](#180-references)

---

## 1.0 Purpose and Rationale

Bacterial plasmids are among the most consequential drivers of antimicrobial resistance (AMR) dissemination and virulence evolution in clinical pathogens [1]. The horizontal transfer of plasmids carrying carbapenemase genes (e.g., *bla*<sub>KPC</sub>, *bla*<sub>NDM</sub>, *bla*<sub>OXA-48</sub>) and virulence determinants (e.g., *iucABCD*, *rmpA*) has been directly linked to the emergence of high-risk clones, including *Klebsiella pneumoniae* ST11 and ST258 [2] [3]. Despite their clinical significance, plasmid analysis remains one of the most challenging aspects of microbial genomics, owing to the structural complexity of these elements, the abundance of repetitive mobile genetic elements (MGEs), and the limitations of short-read sequencing in resolving complete plasmid architectures [4].

This Standard Operating Procedure (SOP) provides a definitive, step-by-step protocol for the comprehensive genomic analysis of bacterial plasmids. It is designed to serve as a **gold-standard reference** for any researcher or laboratory seeking to characterize plasmids at a level suitable for publication in high-impact journals such as *The Lancet Microbe*, *Nature Microbiology*, and *Antimicrobial Agents and Chemotherapy*. The protocol integrates the latest tools, databases, and methodological advances as of 2026, and is grounded in the peer-reviewed literature from the field's leading research groups.

## 2.0 Scope

This SOP is applicable to the genomic analysis of plasmids from bacterial isolates, with a primary focus on clinically relevant Gram-negative pathogens (Enterobacteriaceae, *Acinetobacter*, *Pseudomonas*). The protocol is designed for isolates sequenced using long-read (Oxford Nanopore Technologies or PacBio) and/or short-read (Illumina) platforms. It covers the complete analytical pipeline from raw sequencing data to publication-quality figures, tables, and manuscript-ready text. The protocol is not limited to a single species or sequence type and can be adapted for any bacterial pathogen.

## 3.0 Definitions and Terminology

| Term | Definition |
| :--- | :--- |
| **Plasmid** | An extrachromosomal, self-replicating DNA molecule, typically circular, that can carry accessory genes including those for AMR and virulence. |
| **Mobile Genetic Element (MGE)** | A segment of DNA that can move within or between genomes. This encompasses plasmids, insertion sequences (IS), transposons (Tn), integrons, integrative conjugative elements (ICEs), and genomic islands. |
| **Insertion Sequence (IS)** | The simplest autonomous transposable element, typically 700-2500 bp, encoding only a transposase gene flanked by inverted repeats (IRs). |
| **Transposon (Tn)** | A mobile DNA element larger than an IS, carrying accessory genes (e.g., AMR genes) in addition to transposition machinery. Composite transposons are flanked by IS elements; unit transposons (e.g., Tn*3* family) carry their own *tnpA* and *tnpR* genes. |
| **Integron** | A genetic element that captures and expresses gene cassettes via site-specific recombination mediated by an integrase (*intI*). Class 1 integrons are the most clinically significant. |
| **Replicon** | The minimal genetic unit required for autonomous replication. In plasmid biology, the replicon defines the incompatibility (Inc) group. |
| **Inc Group** | Incompatibility Group: A classification system for plasmids based on their replication and partitioning systems. Two plasmids of the same Inc group cannot coexist stably in the same cell. |
| **Hybrid Assembly** | The computational process of combining long-read and short-read sequencing data to produce a genome assembly that is both complete (from long reads) and accurate (from short reads). |
| **Synteny** | The conserved order and orientation of genes or genomic blocks between different sequences. |
| **Convergent Plasmid** | A plasmid that carries both AMR and virulence determinants, representing a particularly dangerous evolutionary event [5]. |
| **Plasmid Taxonomic Unit (PTU)** | A classification unit for grouping related plasmids based on shared genomic content, analogous to species concepts for organisms [6]. |

## 4.0 Materials: Software, Databases, and Hardware Requirements

### 4.1 Hardware Recommendations

| Component | Minimum | Recommended |
| :--- | :--- | :--- |
| **CPU** | 8 cores | 32+ cores |
| **RAM** | 16 GB | 64-128 GB |
| **Storage** | 500 GB SSD | 2+ TB NVMe SSD |
| **OS** | Ubuntu 20.04+ / CentOS 7+ | Ubuntu 22.04 LTS |

### 4.2 Core Bioinformatics Software

The following table lists all software tools required for this protocol, organized by their role in the analytical pipeline. All tools are open-source unless otherwise noted.

| Phase | Tool | Version | Purpose | Citation |
| :--- | :--- | :--- | :--- | :--- |
| **QC** | FastQC | v0.12.1+ | Short-read quality assessment | [7] |
| **QC** | NanoPlot | v1.42+ | Long-read quality assessment and statistics | [8] |
| **QC** | Filtlong | v0.2.1+ | Long-read quality filtering | [9] |
| **Assembly** | Unicycler | v0.5.0+ | Hybrid assembly (gold standard for plasmid circularization) | [10] |
| **Assembly** | Trycycler | v0.5.5+ | Consensus long-read assembly (highest accuracy) | [11] |
| **Assembly** | Flye | v2.9.3+ | Long-read de novo assembly | [12] |
| **Assembly** | Hybracter | v0.7.0+ | Automated hybrid assembly pipeline | [13] |
| **Assembly** | Plassembler | v1.6.0+ | Dedicated plasmid assembly from hybrid data | [14] |
| **Assembly QC** | QUAST | v5.2.0+ | Assembly quality assessment | [15] |
| **Assembly QC** | Bandage | v0.9.0+ | Assembly graph visualization | [16] |
| **Annotation** | Bakta | v1.9.3+ | Rapid, standardized genome annotation | [17] |
| **Annotation** | Prokka | v1.14.6+ | Genome annotation (alternative to Bakta) | [18] |
| **AMR/Virulence** | ABRicate | v1.0.1+ | Mass screening for AMR, virulence, and plasmid genes | [19] |
| **AMR/Virulence** | AMRFinderPlus | v3.12+ | NCBI AMR gene detection | [20] |
| **AMR/Virulence** | Kleborate | v2.4+ | *Klebsiella*-specific genotyping and virulence scoring | [21] |
| **Replicon Typing** | PlasmidFinder | v2.1+ | Replicon-based plasmid Inc group identification | [22] |
| **Plasmid Typing** | MOB-suite | v3.1.0+ | Plasmid reconstruction, replicon typing, MOB typing, mobility prediction | [23] |
| **Plasmid Typing** | COPLA | v1.0+ | Plasmid Taxonomic Unit (PTU) classification | [6] |
| **MGE Detection** | ISfinder | Database | Reference database for insertion sequence identification | [24] |
| **MGE Detection** | MobileElementFinder | v1.1+ | Automated IS and composite transposon detection | [25] |
| **MGE Detection** | IntegronFinder | v2.0+ | Integron and gene cassette identification | [26] |
| **MGE Detection** | mobileOG-db | v1.6+ | Comprehensive MGE protein family database | [27] |
| **MGE Detection** | Beav | v1.0+ | Automated bacterial genome and MGE annotation pipeline | [28] |
| **Visualization** | BRIG | v0.95+ | Circular comparative plasmid visualization | [29] |
| **Visualization** | Clinker | v0.0.28+ | Linear gene cluster comparison figures | [30] |
| **Visualization** | pyGenomeViz | v0.4.2+ | Programmatic genome visualization (Python) | [31] |
| **Visualization** | SnapGene Viewer | Latest | Interactive plasmid map viewer | [32] |
| **Visualization** | DNAPlotter | v18.2+ | Circular genome/plasmid maps | [33] |
| **Clustering** | Mash | v2.3+ | Fast genome/plasmid distance estimation | [34] |
| **Clustering** | pling | v1.0+ | Structural-aware plasmid clustering | [35] |
| **Network** | Cytoscape | v3.10+ | Plasmid network visualization | [36] |

### 4.3 Key Databases

| Database | URL | Description | Last Updated |
| :--- | :--- | :--- | :--- |
| **PLSDB** | https://ccb-microbe.cs.uni-saarland.de/plsdb/ | Comprehensive plasmid sequence database (72,360+ sequences) | May 2024 |
| **NCBI Nucleotide** | https://www.ncbi.nlm.nih.gov/nucleotide/ | Primary sequence repository | Continuous |
| **ISfinder** | https://isfinder.biotoul.fr/ | Reference centre for bacterial insertion sequences | Continuous |
| **CARD** | https://card.mcmaster.ca/ | Comprehensive Antibiotic Resistance Database | Continuous |
| **VFDB** | http://www.mgc.ac.cn/VFs/ | Virulence Factor Database | Continuous |
| **pMLST** | https://pubmlst.org/organisms/plasmid-mlst | Plasmid MLST typing schemes | Continuous |
| **mobileOG-db** | https://github.com/clb21565/mobileOG-db | MGE protein family database | 2022 |

---

## 5.0 Phase 1: Sequencing Strategy and Data Quality Control

### 5.1 Rationale for Sequencing Strategy

The choice of sequencing platform is the single most critical decision in plasmid analysis. Short-read sequencing alone (Illumina) produces fragmented assemblies that cannot resolve the repetitive IS elements and transposons that are hallmarks of AMR plasmids [4]. Long-read sequencing (Oxford Nanopore Technologies [ONT] or PacBio) can span these repetitive regions, enabling the reconstruction of complete, circular plasmid sequences. Hybrid assembly, which combines the completeness of long reads with the base-level accuracy of short reads, represents the **gold standard** for plasmid genomics [37].

> **Critical Decision:** If only short-read data is available, the analysis is fundamentally limited to plasmid *prediction* (identifying which contigs are likely plasmid-derived) rather than plasmid *reconstruction* (obtaining complete, circular sequences). This distinction must be clearly stated in any resulting publication.

### 5.2 Sequencing Depth Recommendations

| Platform | Minimum Depth | Recommended Depth | Notes |
| :--- | :--- | :--- | :--- |
| **Illumina (short-read)** | 50x | 100x | Paired-end 150 bp reads preferred |
| **ONT (long-read)** | 30x | 50-100x | N50 read length >10 kb preferred |
| **PacBio HiFi** | 25x | 50x | Highest per-read accuracy |

### 5.3 Quality Control Steps

**Step 5.3.1: Short-Read QC**

Assess the quality of Illumina reads using FastQC. Key metrics to evaluate include per-base quality scores (Phred >Q30 for >80% of bases), adapter contamination, and GC content distribution. Trim adapters and low-quality bases using Fastp or Trimmomatic.

```bash
# Quality assessment
fastqc -t 8 -o fastqc_results/ sample_R1.fastq.gz sample_R2.fastq.gz

# Adapter trimming and quality filtering
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
      -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz \
      --detect_adapter_pe --qualified_quality_phred 20 \
      --length_required 50 --thread 8
```

**Step 5.3.2: Long-Read QC**

Assess the quality and length distribution of ONT reads using NanoPlot. Filter reads by quality and length using Filtlong to remove short, low-quality reads that can confound assembly.

```bash
# Quality assessment
NanoPlot --fastq sample_nanopore.fastq.gz -o nanoplot_results/ --threads 8

# Quality and length filtering
filtlong --min_length 1000 --keep_percent 95 \
         --target_bases 500000000 sample_nanopore.fastq.gz \
         | gzip > filtered_nanopore.fastq.gz
```

**Step 5.3.3: QC Pass/Fail Criteria**

| Metric | Pass Criteria | Action if Fail |
| :--- | :--- | :--- |
| Illumina Q30 | >80% bases | Re-sequence or aggressive trimming |
| ONT N50 read length | >5 kb | Consider re-extraction or re-sequencing |
| Estimated genome coverage (short) | >50x | Pool additional sequencing runs |
| Estimated genome coverage (long) | >30x | Pool additional sequencing runs |
| Contamination (Kraken2/Mash screen) | <5% non-target reads | Remove contaminant reads |

---

## 6.0 Phase 2: Gold-Standard Genome Assembly

### 6.1 Assembly Strategy Decision Tree

The choice of assembly strategy depends on the available data and the desired level of accuracy. The following decision framework should be applied:

**If both long-read and short-read data are available (Preferred):**

Use **Unicycler** in hybrid mode for a streamlined, automated approach, or **Hybracter** for a more modern, scalable pipeline. For the highest possible accuracy (e.g., for a reference genome), use the **Trycycler** consensus approach followed by short-read polishing.

**If only long-read data is available:**

Use **Flye** for initial assembly, followed by **Medaka** polishing. Be aware that small plasmids (<10 kb) may be missed or duplicated [38]. Use **Plassembler** as a dedicated plasmid assembly tool if hybrid data becomes available later.

**If only short-read data is available (Limited Analysis):**

Use **SPAdes** for assembly. Plasmid contigs must be identified computationally using tools like **MOB-recon**, **mlplasmids**, or **geNomad**. Complete plasmid reconstruction is generally not possible. This limitation must be explicitly acknowledged.

### 6.2 Recommended Assembly Commands

**Option A: Unicycler Hybrid Assembly (Recommended for most users)**

```bash
unicycler -1 trimmed_R1.fastq.gz -2 trimmed_R2.fastq.gz \
          -l filtered_nanopore.fastq.gz \
          -o unicycler_output/ \
          --threads 32 --mode normal
```

**Option B: Hybracter Automated Pipeline**

```bash
hybracter hybrid --input sample_info.csv \
                 --output hybracter_output/ \
                 --threads 32
```

**Option C: Trycycler Consensus (Highest Accuracy)**

```bash
# Step 1: Generate multiple independent assemblies from subsampled long reads
trycycler subsample --reads filtered_nanopore.fastq.gz --out_dir subsampled/ --count 12

# Step 2: Assemble each subsample with different assemblers
flye --nano-raw subsampled/sample_01.fastq --out-dir assemblies/flye_01/
miniasm_and_minipolish.sh subsampled/sample_02.fastq 8 > assemblies/miniasm_02.fasta
raven --threads 8 subsampled/sample_03.fastq > assemblies/raven_03.fasta
# ... repeat for all subsamples

# Step 3: Cluster contigs into replicon groups
trycycler cluster --assemblies assemblies/*.fasta --reads filtered_nanopore.fastq.gz --out_dir trycycler_clusters/

# Step 4: For each cluster, reconcile, MSA, partition, and consensus
trycycler reconcile --reads filtered_nanopore.fastq.gz --cluster_dir trycycler_clusters/cluster_001/
trycycler msa --cluster_dir trycycler_clusters/cluster_001/
trycycler partition --reads filtered_nanopore.fastq.gz --cluster_dirs trycycler_clusters/cluster_*/
trycycler consensus --cluster_dir trycycler_clusters/cluster_001/

# Step 5: Polish with Medaka (long-read) then Polypolish (short-read)
medaka_polish consensus.fasta filtered_nanopore.fastq.gz medaka_output/
polypolish polish medaka_output/consensus.fasta aligned_R1.sam aligned_R2.sam > final_assembly.fasta
```

### 6.3 Assembly Quality Assessment

**Step 6.3.1: QUAST Assessment**

```bash
quast -o quast_results/ final_assembly.fasta --threads 8
```

**Step 6.3.2: Assembly Graph Visualization with Bandage**

Open the assembly graph file (`.gfa`) in Bandage to visually inspect the assembly. Circular contigs (closed loops) represent complete replicons. Linear contigs may represent incomplete assemblies or linear plasmids.

**Step 6.3.3: Assembly QC Metrics**

| Metric | Expected Value | Interpretation |
| :--- | :--- | :--- |
| Number of contigs | 1 (chromosome) + N (plasmids) | Each complete replicon should be a single contig |
| Largest contig | ~5.0-5.5 Mb (for *K. pneumoniae*) | Should correspond to the chromosome |
| Circular contigs | All major replicons | Confirms completeness |
| Total assembly size | ~5.2-5.8 Mb (for *K. pneumoniae*) | Consistent with species |
| GC content | ~57% (for *K. pneumoniae*) | Consistent with species |

---

## 7.0 Phase 3: Complete Plasmid Enumeration and Validation

### 7.1 Rationale

This is a **critical phase** that is often overlooked. Before any downstream analysis, the researcher must definitively determine: (a) how many plasmids are present, (b) which are complete (circular), and (c) which are potentially incomplete or chimeric. This step provides the foundation for all subsequent characterization and is essential for publication credibility.

### 7.2 Procedure

**Step 7.2.1: Identify Circular Contigs**

From the Unicycler or Trycycler output, identify all contigs annotated as circular. In Unicycler, these are marked with `circular=true` in the FASTA header. In Trycycler, each cluster that successfully completes the consensus step represents a circular replicon.

```bash
# Extract circular contigs from Unicycler output
grep "circular=true" unicycler_output/assembly.fasta | sed 's/>//' > circular_contigs.txt
```

**Step 7.2.2: Distinguish Plasmids from Chromosome**

The largest circular contig is typically the chromosome. All other circular contigs are candidate plasmids. Verify by checking size (plasmids are typically 1 kb to 500 kb), GC content (plasmids may differ from the chromosome), and gene content (presence of replication and transfer genes).

**Step 7.2.3: Validate Completeness**

For each candidate plasmid, confirm completeness by:

1. **Checking for the presence of a replication initiation gene** (*repA* or equivalent) using BLAST against the PlasmidFinder database.
2. **Mapping reads back to the plasmid** using minimap2 and checking for even coverage across the entire sequence. Gaps or regions of zero coverage indicate assembly errors.
3. **Checking for overlapping ends** in the assembly, which is a hallmark of circularization.

```bash
# Map long reads back to plasmid
minimap2 -a plasmid_01.fasta filtered_nanopore.fastq.gz | samtools sort -o plasmid_01_mapped.bam
samtools index plasmid_01_mapped.bam
samtools depth plasmid_01_mapped.bam > plasmid_01_depth.txt
```

**Step 7.2.4: Create the Definitive Plasmid Inventory**

Compile a master table documenting every replicon identified in the isolate. This table forms the backbone of all subsequent analyses and should be included in the manuscript supplementary materials.

| Replicon | Size (bp) | Topology | GC (%) | Coverage (x) | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Chromosome | 5,333,942 | Circular | 57.2 | 85x | Complete |
| pKPN-A | 185,521 | Circular | 52.1 | 120x | Complete |
| pKPN-B | 50,234 | Circular | 46.8 | 95x | Complete |
| pKPN-C | 8,712 | Circular | 50.3 | 210x | Complete |
| pKPN-D | 4,156 | Linear | 48.9 | 45x | Incomplete (fragment) |

---

## 8.0 Phase 4: Annotation and Foundational Characterization

### 8.1 Genome Annotation

**Step 8.1.1: Annotate with Bakta**

Bakta is the recommended annotation tool as it provides standardized, high-quality annotations with consistent gene nomenclature. It outputs GenBank (.gbk) files that are directly compatible with downstream visualization tools (BRIG, Clinker, pyGenomeViz).

```bash
# Annotate the complete genome (chromosome + all plasmids)
bakta --db /path/to/bakta_db/ \
      --output bakta_results/ \
      --prefix sample_name \
      --genus Klebsiella --species "pneumoniae" \
      --strain "sample_ID" \
      --complete \
      --threads 16 \
      final_assembly.fasta
```

**Step 8.1.2: Extract Individual Plasmid Annotations**

Separate the GenBank file into individual files for each replicon. This is necessary for comparative analysis.

```bash
# Use a script or tool to split multi-record GenBank files
# Each plasmid should have its own .gbk file
python3 split_genbank.py bakta_results/sample_name.gbk --output plasmid_annotations/
```

### 8.2 AMR and Virulence Gene Screening

**Step 8.2.1: Screen with ABRicate**

Run ABRicate against multiple databases to identify AMR genes, virulence factors, and plasmid markers. Running against multiple databases provides cross-validation.

```bash
# Screen against CARD (AMR genes)
abricate --db card final_assembly.fasta > abricate_card.tsv

# Screen against VFDB (virulence factors)
abricate --db vfdb final_assembly.fasta > abricate_vfdb.tsv

# Screen against PlasmidFinder (replicon markers)
abricate --db plasmidfinder final_assembly.fasta > abricate_plasmidfinder.tsv

# Generate a summary across all databases
abricate --summary abricate_card.tsv abricate_vfdb.tsv > abricate_summary.tsv
```

**Step 8.2.2: Validate with AMRFinderPlus**

AMRFinderPlus uses a curated NCBI database and provides more detailed annotations, including point mutations conferring resistance.

```bash
amrfinder --nucleotide final_assembly.fasta \
          --protein bakta_results/sample_name.faa \
          --gff bakta_results/sample_name.gff3 \
          --organism "Klebsiella_pneumoniae" \
          --plus --threads 8 \
          --output amrfinder_results.tsv
```

**Step 8.2.3: Map Genes to Specific Replicons**

This is a critical step. For each AMR and virulence gene identified, determine which replicon (chromosome or specific plasmid) it resides on. This is achieved by cross-referencing the contig name in the ABRicate/AMRFinderPlus output with the plasmid inventory from Phase 3.

| Gene | Product | Replicon | Contig Position | Database | % Identity | % Coverage |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| *bla*<sub>KPC-2</sub> | Carbapenemase | pKPN-A | 45,231-46,092 | CARD | 100 | 100 |
| *bla*<sub>NDM-1</sub> | Metallo-beta-lactamase | pKPN-B | 12,456-13,269 | CARD | 100 | 100 |
| *iucA* | Aerobactin biosynthesis | pKPN-A | 98,123-99,876 | VFDB | 99.8 | 100 |
| *rmpA* | Regulator of mucoid phenotype | pKPN-A | 102,345-103,012 | VFDB | 100 | 100 |

---

## 9.0 Phase 5: Replicon Typing and Mobility Classification

### 9.1 Replicon (Inc Group) Typing with PlasmidFinder

PlasmidFinder identifies plasmid replicon types by detecting conserved replication initiation genes. This is the standard method for classifying plasmids into incompatibility groups [22].

```bash
# Run PlasmidFinder on each plasmid individually
plasmidfinder.py -i plasmid_01.fasta -o pf_results/ -p enterobacteriaceae \
                 -l 0.60 -t 0.95
```

### 9.2 Comprehensive Typing with MOB-suite

MOB-suite provides a more comprehensive characterization than PlasmidFinder alone, including relaxase (MOB) typing, mate-pair formation (MPF) typing, and mobility prediction [23].

```bash
# Run MOB-recon on the complete assembly
mob_recon --infile final_assembly.fasta --outdir mob_suite_results/ --num_threads 16

# Run MOB-typer on individual plasmids for detailed typing
mob_typer --infile plasmid_01.fasta --outdir mob_typer_results/
```

### 9.3 Plasmid Taxonomic Unit (PTU) Classification

For deeper evolutionary context, classify plasmids into PTUs using COPLA, which groups plasmids based on pairwise Average Nucleotide Identity (ANI) [6].

```bash
copla.py plasmid_01.fasta -o copla_results/ --threads 8
```

### 9.4 Comprehensive Plasmid Profile Table

Compile all typing results into a single, comprehensive profile for each plasmid. This table is a key deliverable for any publication.

| Plasmid | Size (kb) | Replicon(s) | MOB Type | MPF Type | Mobility | PTU | Key Cargo |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| pKPN-A | 185.5 | IncFII(K), IncFIB(pB171) | MOBF12 | MPFF | Conjugative | PTU-FE | *bla*<sub>KPC-2</sub>, *iucABCD*, *rmpA* |
| pKPN-B | 50.2 | IncX3 | MOBP3 | MPFT | Conjugative | PTU-X3 | *bla*<sub>NDM-1</sub> |
| pKPN-C | 8.7 | ColRNAI | - | - | Non-mobilizable | PTU-C10 | None |

---

## 10.0 Phase 6: Comparative Plasmid Genomics

### 10.1 Objective

The goal of comparative plasmid genomics is to determine the relationships between plasmids from different isolates, identify conserved backbone regions, locate variable accessory regions (often carrying AMR/virulence genes), and identify the "successful" plasmid types that are driving the dissemination of key genes across the population.

### 10.2 Circular Comparison with BRIG

BRIG (BLAST Ring Image Generator) is the standard tool for generating circular comparison figures of complete plasmids [29]. It produces publication-quality images showing regions of homology and divergence between a reference plasmid and multiple query plasmids.

**Procedure:**

1. Select a **reference plasmid**. This should be either a well-characterized reference from NCBI (e.g., pKPC-LK30 for IncFII/KPC plasmids) or the most representative plasmid from your dataset.
2. Prepare query plasmid FASTA files from your isolates.
3. Configure BRIG with the reference and query files.
4. Annotate the outermost ring with key genes (AMR, virulence, IS elements) from the Bakta annotation.
5. Generate the image.

**Interpretation:** Regions with consistent color across all rings represent the conserved plasmid backbone (replication, transfer, maintenance genes). Gaps or regions of different color represent variable accessory regions, which often correspond to MGE-mediated insertions carrying AMR or virulence genes.

### 10.3 Linear Synteny Comparison with Clinker

For comparing specific regions of interest (e.g., an AMR island or a virulence locus) across multiple plasmids, Clinker provides linear synteny diagrams that clearly show gene order, orientation, and homology [30].

```bash
# Compare specific regions from multiple plasmids
clinker plasmid_A.gbk plasmid_B.gbk plasmid_C.gbk \
        -p clinker_output.html \
        --identity 0.3
```

### 10.4 Plasmid Clustering and Network Analysis

For large datasets (>10 plasmids of the same type), cluster plasmids using Mash distance and visualize the relationships using Cytoscape.

```bash
# Calculate pairwise Mash distances
mash sketch -o plasmid_sketches *.fasta
mash dist plasmid_sketches.msh plasmid_sketches.msh > mash_distances.tsv

# Import into Cytoscape for network visualization
# Nodes = plasmids, Edges = Mash distance < threshold (e.g., 0.04)
```

---

## 11.0 Phase 7: Mobile Genetic Element (MGE) Deep-Dive

### 11.1 Rationale

This phase is the **heart of the analysis** for understanding why certain plasmids are successful and how they evolve. MGEs, particularly IS elements, transposons, and integrons, are the molecular machinery that drives the acquisition, rearrangement, and dissemination of AMR and virulence genes on plasmids. A thorough MGE analysis is what distinguishes a descriptive study from a mechanistic, high-impact publication.

### 11.2 Insertion Sequence (IS) Identification

**Step 11.2.1: Automated IS Detection with ISfinder/MobileElementFinder**

Use MobileElementFinder for automated, high-throughput IS detection against the curated ISfinder database [24] [25].

```bash
# Run MobileElementFinder
mobileelement_finder.py -i plasmid_01.fasta -o mef_results/
```

Alternatively, perform a BLASTn search directly against the ISfinder database:

```bash
# Download ISfinder database sequences
# BLASTn search
blastn -query plasmid_01.fasta -db ISfinder_db \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
       -evalue 1e-10 -perc_identity 90 \
       -out is_blast_results.tsv
```

**Step 11.2.2: Manual Curation of IS Elements**

Automated results must be manually curated. For each IS element identified:

1. Verify the presence of **inverted repeats (IRs)** flanking the IS element.
2. Confirm the presence of a complete **transposase gene** (*tnpA*).
3. Determine the **orientation** of the IS element relative to surrounding genes.
4. Identify **direct repeats (DRs)** flanking the IS element, which are generated upon insertion and serve as evidence of a transposition event.
5. Note any **truncated or partial IS elements**, which may indicate ancient insertion events.

**Step 11.2.3: IS Element Inventory**

| IS Element | Family | Size (bp) | Copies on Plasmid | Locations | Flanking Genes | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| IS*26* | IS*6* | 820 | 6 | Multiple | *bla*<sub>KPC-2</sub>, *aph(3')-Ia* | Major driver of rearrangement |
| IS*Kpn26* | IS*5* | 1,195 | 2 | 45,231; 98,765 | *iucA*, *rmpA* | Flanks virulence region |
| IS*Kpn28* | IS*1* | 1,443 | 1 | 123,456 | *bla*<sub>SHV-12</sub> | Upstream of SHV |

### 11.3 Transposon Characterization

**Step 11.3.1: Identify Composite Transposons**

Composite transposons are bounded by two copies of the same IS element (in direct or inverted orientation) and carry cargo genes between them. Identify these by looking for pairs of IS elements of the same type that flank AMR or virulence genes.

**Key Transposons in *K. pneumoniae* ST11:**

| Transposon | Flanking IS | Cargo Genes | Significance |
| :--- | :--- | :--- | :--- |
| Tn*4401* | IS*Kpn7* | *bla*<sub>KPC-2</sub> | Primary KPC carrier |
| Tn*125* | IS*Aba125* | *bla*<sub>NDM-1</sub> | Primary NDM carrier |
| Tn*1999* | IS*1999* | *bla*<sub>OXA-48</sub> | Primary OXA-48 carrier |
| IS*26*-composite | IS*26* | Variable (AMR cassettes) | Major rearrangement driver |

**Step 11.3.2: Identify Unit Transposons**

Unit transposons (e.g., Tn*3* family) carry their own transposase (*tnpA*) and resolvase (*tnpR*) genes. Identify these by searching for *tnpA* and *tnpR* genes in the Bakta annotation and examining the surrounding genetic context.

### 11.4 Integron and Gene Cassette Analysis

**Step 11.4.1: Run IntegronFinder 2.0**

IntegronFinder identifies integrons based on the presence of an integrase gene (*intI*), an *attI* recombination site, and associated gene cassettes [26].

```bash
integron_finder plasmid_01.fasta \
                --local-max \
                --func-annot \
                --outdir integron_results/ \
                --cpu 8
```

**Step 11.4.2: Characterize Gene Cassettes**

For each integron identified, document the gene cassettes and their order. The cassette array is a key indicator of the integron's evolutionary history and resistance profile.

| Integron Class | Integrase | Cassette Array | Resistance Profile |
| :--- | :--- | :--- | :--- |
| Class 1 | *intI1* | *dfrA12*-*aadA2*-*qacE*Δ1-*sul1* | Trimethoprim, Streptomycin, Sulfonamide |
| Class 1 | *intI1* | *bla*<sub>VIM-1</sub>-*aacA4*-*aadA1* | Carbapenem, Aminoglycoside |

### 11.5 Comprehensive MGE Annotation with Beav

For an automated, all-in-one approach, use the Beav pipeline, which integrates multiple annotation tools (Bakta, ISfinder, IntegronFinder, ICEfinder, and others) into a single workflow [28].

```bash
beav --input final_assembly.fasta --output beav_results/ --threads 16
```

---

## 12.0 Phase 8: Genomic Environment and Synteny Analysis

### 12.1 Rationale

This phase produces the **key figures** for any high-impact publication on plasmid biology. Genomic environment analysis involves extracting a region of interest (e.g., the flanking region of *bla*<sub>KPC-2</sub>) from multiple isolates and comparing them to reveal the MGE-mediated rearrangements that have shaped the evolution of these regions. This is the analysis that tells the **evolutionary story** of how a gene arrived on a plasmid and how it has been mobilized.

### 12.2 Extracting Regions of Interest

**Step 12.2.1: Identify the Target Region**

Select the gene or region of interest (e.g., *bla*<sub>KPC-2</sub>, *iucABCD* operon, or a large AMR island). Use the Bakta annotation to determine its exact coordinates on the plasmid.

**Step 12.2.2: Extract Flanking Regions**

Extract the target gene plus a defined flanking region (e.g., 10-20 kb upstream and downstream) from each plasmid. This can be done using Samtools, Bedtools, or custom scripts.

```bash
# Extract a 40 kb region centered on bla_KPC-2 (e.g., positions 25,000-65,000)
samtools faidx plasmid_01.fasta contig_name:25000-65000 > kpc_region_isolate1.fasta
```

**Step 12.2.3: Re-annotate Extracted Regions**

Re-annotate the extracted regions with Bakta to generate GenBank files suitable for visualization.

```bash
bakta --db /path/to/bakta_db/ --output region_annotation/ kpc_region_isolate1.fasta
```

### 12.3 Linear Comparison with Clinker

Clinker is the recommended tool for generating publication-quality linear comparison figures. It automatically identifies homologous genes between input GenBank files and draws links between them [30].

```bash
clinker kpc_region_isolate1.gbk kpc_region_isolate2.gbk \
        kpc_region_isolate3.gbk kpc_region_reference.gbk \
        -p kpc_comparison.html \
        --identity 0.3
```

The resulting HTML file is interactive and can be customized. Export as SVG for publication.

### 12.4 Advanced Visualization with pyGenomeViz

For highly customized figures with precise control over colors, labels, and layout, use pyGenomeViz [31]. This Python library allows programmatic generation of genome comparison figures.

```python
from pygenomeviz import GenomeViz

gv = GenomeViz(fig_width=15, fig_track_height=0.7, feature_track_ratio=0.3)

# Add tracks for each isolate's region
track1 = gv.add_feature_track("Isolate_1", size=40000)
track1.add_feature(start=15000, end=15862, strand=1, label="bla_KPC-2", 
                   facecolor="red", linewidth=1)
track1.add_feature(start=14180, end=15000, strand=1, label="ISKpn7", 
                   facecolor="orange", linewidth=1)
# ... add more features

# Add similarity links between tracks
gv.add_link(("Isolate_1", 1000, 20000), ("Isolate_2", 1000, 20000), 
            color="grey", inverted=False, curve=True)

gv.savefig("kpc_genomic_environment.png", dpi=300)
```

### 12.5 Circular Plasmid Maps

For complete plasmid visualization, use DNAPlotter or SnapGene Viewer to create circular maps showing all annotated features (CDS, IS elements, AMR genes, virulence genes, transfer genes) [33].

---

## 13.0 Phase 9: Plasmid Epidemiology and Transmission Analysis

### 13.1 Identifying Plasmid Transmission Events

When analyzing a collection of isolates (e.g., from a hospital outbreak or a population study), it is essential to determine whether the same plasmid is shared across multiple isolates, which would indicate horizontal transfer.

**Step 13.1.1: Pairwise Plasmid Comparison**

Use Mash to calculate pairwise distances between all plasmids of the same Inc type across your dataset [34].

```bash
mash dist -t plasmid_sketches.msh plasmid_sketches.msh > pairwise_distances.tsv
```

**Step 13.1.2: Define Similarity Threshold**

A commonly used threshold for inferring recent plasmid transmission is **>95% gene content at >99% nucleotide identity** [4]. For Mash-based analysis, a distance of **<0.01** (corresponding to ~99% ANI) is a reasonable starting point for identifying highly similar plasmids.

**Step 13.1.3: Structural-Aware Clustering with pling**

For a more sophisticated analysis that accounts for structural rearrangements (insertions, deletions, inversions) common in plasmids, use pling [35].

```bash
pling cluster --input plasmid_list.txt --output pling_results/ --threads 16
```

### 13.2 Plasmid Network Visualization

Visualize plasmid relationships as a network using Cytoscape, where nodes represent plasmids (colored by isolate, host, or location) and edges represent similarity above the defined threshold [36].

---

## 14.0 Phase 10: Interpretation, Reporting, and Publication Standards

### 14.1 Key Questions to Address

A high-impact plasmid analysis should address the following questions:

1. **How many distinct plasmids are present?** What are their Inc types, sizes, and mobility?
2. **Which plasmids carry the key AMR and virulence genes?** Are these genes chromosomal or plasmid-borne?
3. **Are there convergent plasmids?** Do any plasmids carry both AMR and virulence genes?
4. **What MGEs are driving gene acquisition?** Which IS elements, transposons, and integrons are present?
5. **Is there evidence of plasmid transmission?** Are the same plasmids shared across isolates from different patients, wards, or hospitals?
6. **What makes these plasmids "successful"?** What structural features (e.g., conjugation machinery, addiction systems, IS-mediated flexibility) contribute to their persistence and spread?

### 14.2 Required Figures for Publication

| Figure | Content | Tool |
| :--- | :--- | :--- |
| **Figure 1** | Circular plasmid map with all annotated features | DNAPlotter / SnapGene |
| **Figure 2** | BRIG circular comparison of plasmids across isolates | BRIG |
| **Figure 3** | Genomic environment of key AMR gene (e.g., *bla*<sub>KPC-2</sub>) | Clinker / pyGenomeViz |
| **Figure 4** | Genomic environment of key virulence locus (e.g., *iuc* operon) | Clinker / pyGenomeViz |
| **Figure 5** | Plasmid transmission network | Cytoscape |
| **Supplementary** | Complete plasmid inventory table | Manual compilation |
| **Supplementary** | Full IS element and MGE inventory | ISfinder / IntegronFinder |

### 14.3 Manuscript Methods Section Template

The following is a template for the Methods section of a manuscript:

> **Plasmid Analysis.** Complete bacterial genomes were assembled using Unicycler v0.5.0 [10] in hybrid mode, combining Illumina short reads and Oxford Nanopore long reads. Assembly quality was assessed using QUAST v5.2.0 [15] and Bandage v0.9.0 [16]. Genome annotation was performed using Bakta v1.9.3 [17]. Antimicrobial resistance genes were identified using ABRicate v1.0.1 [19] against the CARD database [39] and validated with AMRFinderPlus v3.12 [20]. Virulence factors were identified using the VFDB [40]. Plasmid replicon typing was performed using PlasmidFinder v2.1 [22], and comprehensive plasmid characterization, including relaxase typing and mobility prediction, was performed using MOB-suite v3.1.0 [23]. Insertion sequences were identified by BLASTn search against the ISfinder database [24] and validated using MobileElementFinder v1.1 [25]. Integrons were detected using IntegronFinder v2.0 [26]. Comparative plasmid analysis was performed using BRIG v0.95 [29] for circular comparisons and Clinker v0.0.28 [30] for linear synteny analysis. Plasmid clustering was performed using Mash v2.3 [34] and visualized using Cytoscape v3.10 [36].

---

## 15.0 Quality Assurance and Control

### 15.1 QC Checkpoints

| Phase | Checkpoint | Pass Criteria | Action if Fail |
| :--- | :--- | :--- | :--- |
| Assembly | Circular contigs identified | All major replicons circular | Re-sequence with long reads |
| Assembly | Read mapping coverage | Even coverage, no gaps >100 bp | Investigate assembly errors |
| Annotation | Key genes detected | All expected AMR/virulence genes found | Check assembly completeness |
| Replicon Typing | Inc type assigned | At least one replicon type per plasmid | Check against PLSDB |
| IS Detection | IS elements validated | IRs and transposase confirmed | Manual BLAST verification |
| Integrons | Cassette arrays complete | Full array from *intI* to *qacE*Δ1 | Check assembly for breaks |

### 15.2 Common Pitfalls and Solutions

| Pitfall | Cause | Solution |
| :--- | :--- | :--- |
| Missing small plasmids | Long-read assemblers struggle with <10 kb plasmids | Use Plassembler or check raw reads |
| Chimeric contigs | IS elements cause mis-assembly | Verify with read mapping; use Trycycler |
| Duplicated plasmids | Flye may duplicate circular sequences | Check for redundant contigs with Mash |
| Incomplete IS annotation | Truncated IS at contig boundaries | Use long-read assembly to resolve |
| False-positive plasmid prediction | Chromosomal regions misclassified | Verify with read depth and gene content |

---

## 16.0 Decision Trees

### 16.1 Sequencing Strategy Decision

```
START: What sequencing data is available?
│
├── Long-read + Short-read (Hybrid) → BEST OPTION
│   ├── Need highest accuracy? → Trycycler + Polypolish
│   ├── Need automated pipeline? → Hybracter
│   └── Standard analysis? → Unicycler (hybrid mode)
│
├── Long-read only
│   ├── ONT reads → Flye + Medaka polishing
│   ├── PacBio HiFi → Flye (HiFi mode)
│   └── WARNING: Small plasmids may be missed
│
└── Short-read only → LIMITED ANALYSIS
    ├── Assembly: SPAdes
    ├── Plasmid prediction: MOB-recon / geNomad / mlplasmids
    └── WARNING: Complete plasmid reconstruction NOT possible
```

### 16.2 Plasmid Characterization Decision

```
START: Is the plasmid complete (circular)?
│
├── YES (Complete circular plasmid)
│   ├── Replicon typing → PlasmidFinder + MOB-typer
│   ├── Full annotation → Bakta
│   ├── AMR/Virulence screening → ABRicate + AMRFinderPlus
│   ├── IS element detection → ISfinder + MobileElementFinder
│   ├── Integron detection → IntegronFinder
│   ├── Comparative analysis → BRIG + Clinker
│   └── Genomic environment → pyGenomeViz / Clinker
│
└── NO (Draft/fragmented plasmid)
    ├── Plasmid prediction → MOB-recon
    ├── Replicon typing → PlasmidFinder (may be incomplete)
    ├── AMR gene detection → ABRicate (gene presence only)
    └── WARNING: Genomic context analysis is LIMITED
        └── Consider long-read sequencing for key isolates
```

### 16.3 MGE Analysis Decision

```
START: What type of MGE analysis is needed?
│
├── IS Element Analysis
│   ├── Automated detection → MobileElementFinder
│   ├── Database search → BLASTn vs ISfinder
│   ├── Manual curation → Check IRs, DRs, transposase
│   └── Output: IS inventory table with locations
│
├── Transposon Analysis
│   ├── Composite Tn → Identify paired IS elements flanking cargo
│   ├── Unit Tn → Search for tnpA/tnpR genes (Tn3 family)
│   ├── Known Tn → Tn4401 (KPC), Tn125 (NDM), Tn1999 (OXA-48)
│   └── Output: Transposon structure diagrams
│
├── Integron Analysis
│   ├── Detection → IntegronFinder 2.0
│   ├── Cassette identification → Gene cassette array
│   └── Output: Integron class, cassette order, resistance profile
│
└── Comprehensive MGE Pipeline
    └── Beav → Automated all-in-one annotation
```

---

## 17.0 Data Management and Reproducibility

### 17.1 Directory Structure

A standardized directory structure ensures reproducibility and facilitates collaboration:

```
project_name/
├── 00_raw_data/
│   ├── illumina/
│   └── nanopore/
├── 01_qc/
│   ├── fastqc/
│   ├── nanoplot/
│   └── filtered_reads/
├── 02_assembly/
│   ├── unicycler/
│   └── assembly_qc/
├── 03_annotation/
│   ├── bakta/
│   └── individual_replicons/
├── 04_plasmid_typing/
│   ├── plasmidfinder/
│   ├── mob_suite/
│   └── copla/
├── 05_amr_virulence/
│   ├── abricate/
│   └── amrfinderplus/
├── 06_mge_analysis/
│   ├── isfinder/
│   ├── integron_finder/
│   └── mobile_element_finder/
├── 07_comparative/
│   ├── brig/
│   ├── clinker/
│   └── mash/
├── 08_visualization/
│   ├── plasmid_maps/
│   ├── genomic_environment/
│   └── networks/
├── 09_results/
│   ├── tables/
│   ├── figures/
│   └── supplementary/
└── scripts/
    └── analysis_pipeline.sh
```

### 17.2 Data Submission Requirements

All data generated under this protocol should be submitted to public repositories:

| Data Type | Repository | Format |
| :--- | :--- | :--- |
| Raw sequencing reads | NCBI SRA | FASTQ |
| Assembled genomes | NCBI GenBank | FASTA + GenBank |
| Complete plasmid sequences | NCBI GenBank | GenBank |
| AMR gene profiles | NCBI BioSample | Antibiogram metadata |

---

## 18.0 References

[1] Murray, C.J.L., et al. (2022). Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. *The Lancet*, 399:629-655.

[2] Wyres, K.L., et al. (2020). Genomic surveillance for hypervirulence and multi-drug resistance in invasive *Klebsiella pneumoniae* from South and Southeast Asia. *Genome Medicine*, 12:11.

[3] Gibbon, M.J., et al. (2026). Convergence and global molecular epidemiology of *Klebsiella pneumoniae* plasmids harbouring the *iuc3* virulence locus. *The Lancet Microbe*, 7(2):101236.

[4] Beh, J.Q., et al. (2025). Challenges and considerations for whole-genome-based antimicrobial resistance plasmid investigations. *Antimicrobial Agents and Chemotherapy*, 69(12).

[5] Lam, M.M.C., et al. (2019). A genomic surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its related species complex. *Nature Communications*, 12:4188.

[6] Redondo-Salvo, S., et al. (2021). COPLA, a taxonomic classifier of plasmids. *BMC Bioinformatics*, 22:390.

[7] Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data. Babraham Bioinformatics.

[8] De Coster, W., et al. (2018). NanoPlot: long-read sequencing data visualization and quality control. *Bioinformatics*, 34(15):2666-2669.

[9] Wick, R.R. (2021). Filtlong: quality filtering tool for long reads. GitHub. https://github.com/rrwick/Filtlong

[10] Wick, R.R., et al. (2017). Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. *PLoS Computational Biology*, 13(6):e1005595.

[11] Wick, R.R., et al. (2021). Trycycler: consensus long-read assemblies for bacterial genomes. *Genome Biology*, 22(1):266.

[12] Kolmogorov, M., et al. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, 37:540-546.

[13] Bouras, G., et al. (2024). Hybracter: enabling scalable, automated, complete and accurate bacterial genome assemblies. *Microbial Genomics*, 10(1).

[14] Bouras, G., et al. (2023). Plassembler: an automated bacterial plasmid assembly tool. *Bioinformatics*, 39(7):btad409.

[15] Gurevich, A., et al. (2013). QUAST: quality assessment tool for genome assemblies. *Bioinformatics*, 29(8):1072-1075.

[16] Wick, R.R., et al. (2015). Bandage: interactive visualization of de novo genome assemblies. *Bioinformatics*, 31(20):3350-3352.

[17] Schwengers, O., et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics*, 7(11).

[18] Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. *Bioinformatics*, 30(14):2068-2069.

[19] Seemann, T. (2024). ABRicate: mass screening of contigs for antimicrobial resistance or virulence genes. GitHub. https://github.com/tseemann/abricate

[20] Feldgarden, M., et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*, 11:12728.

[21] Lam, M.M.C., et al. (2021). A genomic surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its related species complex. *Nature Communications*, 12:4188.

[22] Carattoli, A., et al. (2014). In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrobial Agents and Chemotherapy*, 58(7):3895-3903.

[23] Robertson, J., & Nash, J.H.E. (2018). MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. *Microbial Genomics*, 4(8):e000206.

[24] Siguier, P., et al. (2006). ISfinder: the reference centre for bacterial insertion sequences. *Nucleic Acids Research*, 34(Database issue):D32-D36.

[25] Johansson, M.H.K., et al. (2021). Detection of mobile genetic elements associated with antibiotic resistance in *Salmonella enterica* using a newly developed web tool: MobileElementFinder. *Journal of Antimicrobial Chemotherapy*, 76(1):101-109.

[26] Néron, B., et al. (2022). IntegronFinder 2.0: identification and analysis of integrons across bacteria. *Nucleic Acids Research*, 50(W1):W530-W536.

[27] Brown, C.L., et al. (2022). mobileOG-db: a manually curated database of protein families mediating the life cycle of bacterial mobile genetic elements. *Applied and Environmental Microbiology*, 88(18):e00991-22.

[28] Jung, J.M., et al. (2024). Beav: a bacterial genome and mobile element annotation pipeline. *mSphere*, 9(7):e00209-24.

[29] Alikhan, N.F., et al. (2011). BLAST Ring Image Generator (BRIG): simple prokaryote genome comparisons. *BMC Genomics*, 12:402.

[30] Gilchrist, C.L.M., & Chooi, Y.H. (2021). clinker & clustermap.js: automatic generation of gene cluster comparison figures. *Bioinformatics*, 37(16):2473-2475.

[31] Shimoyama, Y. (2022). pyGenomeViz: a genome visualization python package for comparative genomics. GitHub. https://github.com/moshi4/pyGenomeViz

[32] SnapGene. SnapGene Viewer. https://www.snapgene.com/snapgene-viewer

[33] Carver, T., et al. (2009). DNAPlotter: circular and linear interactive genome visualization. *Bioinformatics*, 25(1):119-120.

[34] Ondov, B.D., et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*, 17:132.

[35] Barilar, I., et al. (2024). pling: epidemiological clustering of plasmids accounting for structural events. *bioRxiv*.

[36] Shannon, P., et al. (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks. *Genome Research*, 13(11):2498-2504.

[37] van Almsick, V., et al. (2026). Long-read sequencing for bacterial plasmid analysis: a brief overview. *FEMS Microbiology Letters*, fnag014.

[38] Wick, R.R., et al. (2023). Long read genome assemblers struggle with small plasmids. *Microbial Genomics*, 9(5).

[39] Alcock, B.P., et al. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Research*, 51(D1):D419-D427.

[40] Liu, B., et al. (2022). VFDB 2022: a general classification scheme for bacterial virulence factors. *Nucleic Acids Research*, 50(D1):D912-D917.

[41] Liu, X., et al. (2025). ST11 carbapenem-resistant *Klebsiella pneumoniae* integrates virulence plasmid fragments into the chromosome via insertion sequence. *BMC Microbiology*, 25:493.

---

**End of Document**

*This SOP is intended as a living document and should be updated as new tools, databases, and methodological advances become available. Last reviewed: February 2026.*
