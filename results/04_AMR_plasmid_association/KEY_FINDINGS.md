# Key Findings: AMR Gene-Plasmid Associations in ST11 K. pneumoniae (Saudi Arabia)

## Analysis Date: 2026-02-15
## Samples Processed: 5/8 (JKPB244, JKPB284, JKPB381, JKPR282, KP21802)
## Pending: KP21915, MDRKP121, MDRKP122

## Critical Finding: Dual Carbapenemase Producers

### Carbapenemase Distribution

| Isolate | blaKPC-2 | blaNDM-1 | blaOXA-232 | Dual CP | KPC Plasmid | NDM/OXA Plasmid |
|---------|----------|----------|------------|---------|-------------|-----------------|
| JKPB244 | Yes | Yes | No | KPC+NDM | ~130 kb IncFII+IncR | ~363 kb IncFIB(pNDM-Mar)+IncHI1B |
| JKPB284 | Yes | Yes | No | KPC+NDM | ~130 kb IncFII+IncR | ~362 kb IncFIB(pNDM-Mar)+IncHI1B |
| JKPB381 | Yes | No | Yes | KPC+OXA | ~130 kb IncFII+IncR | ~353 kb plasmid |
| JKPR282 | Yes | No | No | No | ~130 kb IncFII+IncR | ~353 kb plasmid |
| KP21802 | Yes | Yes | No | KPC+NDM | ~130 kb IncFII+IncR | ~360 kb IncFIB(pNDM-Mar)+IncHI1B |

### Key Statistics
- **blaKPC-2**: 5/5 (100%) — universal in this ST11 lineage
- **blaNDM-1**: 3/5 (60%) — co-carried on the large virulence plasmid
- **blaOXA-232**: 1/5 (20%) — JKPB381 carries OXA-232 instead of NDM-1
- **Dual carbapenemase**: 4/5 (80%) carry two carbapenemases

### blaKPC-2 (Carbapenemase)
- **Location**: ~130 kb IncFII(pHN7A8)+IncR plasmid (contig_3)
- **Conserved position**: pos ~42,346 in JKPB244/JKPB284, pos ~40,251 in JKPB381
- **Prevalence**: 5/5 (100%) of processed samples
- **Co-located with**: blaCTX-M-65, blaSHV-12, blaTEM-1, rmtB1, fosA3, catA2

### blaNDM-1 (Metallo-beta-lactamase)
- **Location**: ~360 kb IncFIB(pNDM-Mar)+IncHI1B(pNDM-MAR) plasmid (contig_2)
- **Position**: pos ~218,912 in JKPB244, pos ~218,016 in JKPB284
- **Prevalence**: 3/5 (60%)
- **Co-located with**: armA, ble, msr(E), mph(E), mph(A), qnrS1, sul1, sul2, dfrA5

## Plasmid-Borne AMR Genes per Isolate

| Isolate | Plasmid AMR | Chromosome AMR | Total AMR | MGE Elements |
|---------|-------------|----------------|-----------|--------------|
| JKPB244 | 67 | 281 | 348 | 457 |
| JKPB284 | 68 | 284 | 352 | 468 |
| JKPB381 | 38 | 276 | 314 | 440 |
| JKPR282 | 44 | 280 | 324 | 458 |
| KP21802 | 65 | 280 | 345 | 466 |

## Plasmid Architecture Summary

### Contig_2 (Large plasmid, ~330-363 kb)
- **Replicon**: IncFIB(pNDM-Mar) + IncHI1B(pNDM-MAR)
- **AMR genes**: blaNDM-1, armA, ble, mph(A), mrx(A), sul1, qacEdelta1, dfrA5, aph(3')-Ia, msr(E), mph(E), qnrS1, aph(3')-VI, sul2
- **Heavy metals**: terABCDEWXZ (tellurium), merCPTR (mercury)
- **Virulence**: iucABCD, iutA (aerobactin siderophore), rmpA, rmpC (mucoid phenotype)

### Contig_3 (Medium plasmid, ~105-131 kb)
- **Replicon**: IncFII(pHN7A8) + IncR
- **AMR genes**: blaKPC-2, blaCTX-M-65, blaSHV-12, blaTEM-1, rmtB1, fosA3, catA2
- **Critical**: This is the KPC-carrying plasmid

### Contig_4 (~10 kb)
- **Replicon**: ColRNAI
- **No AMR genes detected**

### Contig_5 (~5.6 kb)
- **Replicon**: Col-type
- **No AMR genes detected**

## IntegronFinder Results

| Isolate | Complete Integrons | CALIN | In0 | Location |
|---------|-------------------|-------|-----|----------|
| JKPB244 | 1 | 0 | 0 | contig_2 (~363kb plasmid) |
| JKPB284 | 1 | 0 | 0 | contig_2 (~362kb plasmid) |
| JKPB381 | 0 | 0 | 0 | None |
| JKPR282 | 1 | 0 | 0 | contig_2 (~353kb plasmid) |

## Clinical Significance
- These ST11 K. pneumoniae isolates represent a **high-risk convergent lineage** combining:
  - Multiple carbapenemases (KPC-2 + NDM-1/OXA-232)
  - Extended-spectrum beta-lactamases (CTX-M-65, SHV-12, TEM-1)
  - Hypervirulence markers (aerobactin, rmpA/rmpC)
  - Heavy metal resistance (tellurium, mercury)
- The conserved plasmid architecture suggests **clonal spread** of this lineage in Saudi Arabia
