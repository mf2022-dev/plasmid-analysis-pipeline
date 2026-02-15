# Changelog

All notable changes to PlasmidScope will be documented in this file.

## [1.0.0] - 2026-02-15

### Added
- Initial release of PlasmidScope pipeline
- 8-phase analysis workflow (QC -> Assembly -> Annotation -> AMR/Virulence -> Typing -> MGE -> Comparative -> Reporting)
- Support for hybrid, long-read, and short-read assembly modes
- 22 Nextflow process modules
- Docker, Singularity, and Conda profiles
- SLURM, PBS, and SGE HPC support
- Comprehensive HTML report generation
- Publication-quality figure generation
- Interactive plasmid network visualization
- Companion SOP document
- Test data and CI/CD pipeline

### Tools Integrated
- FastQC, NanoPlot, Fastp, Filtlong (QC)
- Unicycler, Flye, Hybracter (Assembly)
- Bakta (Annotation)
- ABRicate, AMRFinderPlus (AMR/Virulence)
- PlasmidFinder, MOB-suite (Plasmid Typing)
- ISfinder, IntegronFinder, MobileElementFinder (MGE)
- BRIG, Clinker, Mash (Comparative)
- MultiQC (Reporting)
