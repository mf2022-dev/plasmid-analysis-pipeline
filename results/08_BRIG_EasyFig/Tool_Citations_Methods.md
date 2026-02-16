# Software Tools and Citations for Genomic Visualization

## Tools Used in This Analysis

The following table summarizes all bioinformatics tools used to generate the resistance gene environment figures, along with their versions, purposes, and formal citations for inclusion in the manuscript methods section.

| Tool | Version | Purpose | Citation |
|------|---------|---------|----------|
| **Bakta** | v1.9.4 | Rapid and standardized annotation of bacterial genomes, including CDS, IS elements, and AMR gene identification | Schwengers O, Jelonek L, Diber MA, Zurbr&uuml;gg S,3 Goesmann A. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom*. 2021;7(11):000685. doi:10.1099/mgen.0.000685 |
| **AMRFinderPlus** | v3.12.8 | Identification of antimicrobial resistance genes, stress response genes, and virulence factors using the NCBI Reference Gene Database | Feldgarden M, Brover V, Gonzalez-Escalona N, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep*. 2021;11:12728. doi:10.1038/s41598-021-91456-0 |
| **ISfinder** | Database (2006) | Reference database for bacterial insertion sequences used by Bakta for IS element classification | Siguier P, Perochon J, Lestrade L, Mahillon J, Chandler M. ISfinder: the reference centre for bacterial insertion sequences. *Nucleic Acids Res*. 2006;34(Database issue):D32-D36. doi:10.1093/nar/gkj014 |
| **IntegronFinder** | v2.0.6 | Detection and annotation of integrons in bacterial genomes | N&eacute;ron B, Littner E, Haudiquet M, Perrin A, Cury J, Rocha EPC. IntegronFinder 2.0: identification and analysis of integrons across bacteria, with a focus on antibiotic resistance in Klebsiella. *Microorganisms*. 2022;10(4):700. doi:10.3390/microorganisms10040700 |
| **HMMER** | v3.3.2 | Profile hidden Markov model searches for integron attC sites and protein domain identification | Eddy SR. Accelerated profile HMM searches. *PLoS Comput Biol*. 2011;7(10):e1002195. doi:10.1371/journal.pcbi.1002195 |
| **INFERNAL** | v1.1.5 | RNA secondary structure-based covariance model searches for non-coding RNA identification | Nawrocki EP, Eddy SR. Infernal 1.1: 100-fold faster RNA homology searches. *Bioinformatics*. 2013;29(22):2933-2935. doi:10.1093/bioinformatics/btt509 |
| **Prodigal** | v2.6.3 | Ab initio gene prediction in prokaryotic genomes | Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics*. 2010;11:119. doi:10.1186/1471-2105-11-119 |
| **DIAMOND** | v2.1.8 | Accelerated protein alignment for functional annotation of predicted CDS | Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. *Nat Methods*. 2015;12(1):59-60. doi:10.1038/nmeth.3176 |
| **BLASTn** | v2.15.0+ | Nucleotide sequence comparison for BRIG-style circular maps and EasyFig-style linear comparisons | Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421. doi:10.1186/1471-2105-10-421 |
| **PlasmidFinder** | v2.1 (web) | Replicon typing and plasmid incompatibility group identification | Carattoli A, Zankari E, Garcia-Fernandez A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob Agents Chemother*. 2014;58(7):3895-3903. doi:10.1128/AAC.02412-14 |
| **Python matplotlib** | v3.x | Generation of publication-quality scientific figures | Hunter JD. Matplotlib: a 2D graphics environment. *Comput Sci Eng*. 2007;9(3):90-95. doi:10.1109/MCSE.2007.55 |
| **pyCirclize** | v1.10.1 | Circular genome visualization (BRIG-style maps) | Shimoyama Y. pyCirclize: Circular visualization in Python. 2022. Available at: https://github.com/moshi4/pyCirclize |
| **pyGenomeViz** | v1.6.1 | Linear genome comparison visualization (EasyFig-style maps) | Shimoyama Y. pyGenomeViz: A genome visualization python package for comparative genomics. 2022. Available at: https://github.com/moshi4/pyGenomeViz |

## Suggested Methods Text

> **Genome annotation and resistance gene identification.** Long-read genome assemblies were annotated using Bakta v1.9.4 [1] with the light database (v5.1), which integrates AMRFinderPlus v3.12.8 [2] for antimicrobial resistance gene identification against the NCBI Reference Gene Catalog, and the ISfinder database [3] for insertion sequence classification. Integrons were detected using IntegronFinder v2.0.6 [4] with HMMER v3.3.2 [5] and INFERNAL v1.1.5 [6] for attC site identification. Plasmid replicon types were determined using PlasmidFinder v2.1 [7].
>
> **Comparative genomic visualization.** Circular plasmid comparison maps were generated using pyCirclize v1.10.1 [8] with BLASTn v2.15.0+ [9] for pairwise nucleotide identity calculations across isolates. Linear genetic environment maps of resistance gene cassettes were created using matplotlib v3.x [10] with gene coordinates extracted from Bakta annotations. Genes were color-coded by functional category: antimicrobial resistance genes (red), IS elements and transposases (blue), integrases and integrons (orange), mercury resistance genes (grey), and other coding sequences (light grey). Transposon cassettes were delineated based on flanking IS elements and structural boundaries identified from the annotation data.

## References

1. Schwengers O, et al. *Microb Genom*. 2021;7(11):000685.
2. Feldgarden M, et al. *Sci Rep*. 2021;11:12728.
3. Siguier P, et al. *Nucleic Acids Res*. 2006;34(Database issue):D32-D36.
4. Neron B, et al. *Microorganisms*. 2022;10(4):700.
5. Eddy SR. *PLoS Comput Biol*. 2011;7(10):e1002195.
6. Nawrocki EP, Eddy SR. *Bioinformatics*. 2013;29(22):2933-2935.
7. Carattoli A, et al. *Antimicrob Agents Chemother*. 2014;58(7):3895-3903.
8. Shimoyama Y. pyCirclize. 2022. https://github.com/moshi4/pyCirclize
9. Camacho C, et al. *BMC Bioinformatics*. 2009;10:421.
10. Hunter JD. *Comput Sci Eng*. 2007;9(3):90-95.
