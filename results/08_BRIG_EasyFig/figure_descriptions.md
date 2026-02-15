# BRIG and EasyFig-style Visualization Descriptions

## Figure: EasyFig_KPC2_Environment
- **Type**: EasyFig-style linear comparison (pyGenomeViz + BLAST)
- **Content**: 30kb genomic environment of blaKPC-2 across all 8 isolates
- **Key observations**:
  - Highly conserved KPC-2 environment across all 8 isolates
  - blaKPC-2 (red) flanked by transposases (purple) and IS elements
  - blaSHV-12 and blaCTX-M-65 also present on same plasmid
  - BLAST similarity links (blue shading) show >99% identity between isolates
  - JKPB381 shows slight structural variation in the KPC region

## Figure: EasyFig_NDM1_Environment
- **Type**: EasyFig-style linear comparison (pyGenomeViz + BLAST)
- **Content**: 30kb genomic environment of blaNDM-1 across 6 NDM+ isolates
- **Key observations**:
  - Conserved NDM-1 cassette: sul1-dfrA5-tnp-aph(3')-Ia-tnp-ble-blaNDM-1-tnp-armA
  - armA (16S rRNA methylase) consistently downstream of NDM-1
  - mph(E) and msr(E) (macrolide resistance) further downstream
  - BLAST links show high conservation but some structural rearrangements between isolates
  - MDRKP122 shows truncation at right end of the region

## Figure: BRIG_KPC_Plasmid
- **Type**: BRIG-style circular comparison map (pyCirclize + BLAST)
- **Content**: Full ~131kb KPC plasmid (IncFII+IncR) with 7 comparison rings
- **Reference**: KP21915
- **Key observations**:
  - All 7 other isolates show near-complete coverage of the reference plasmid
  - Gaps in some rings indicate structural variations or deletions
  - AMR gene cluster visible at ~5kb (blaCTX-M-65) and ~42kb (blaKPC-2)
  - Transposase/IS elements distributed throughout the plasmid

## Figure: BRIG_NDM_Plasmid
- **Type**: BRIG-style circular comparison map (pyCirclize + BLAST)
- **Content**: Full ~354kb NDM plasmid (IncFIB(pNDM-Mar)+IncHI1B) with 5 comparison rings
- **Reference**: KP21915
- **Key observations**:
  - Large plasmid with extensive conservation across all 6 NDM+ isolates
  - AMR gene cluster concentrated around ~210-220kb region
  - Some isolates show gaps indicating structural variations
  - Conjugation/transfer genes visible in multiple regions
