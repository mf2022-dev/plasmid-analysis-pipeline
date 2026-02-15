#!/usr/bin/env python3
"""
Extract AMR genes and MGE elements from Bakta TSV annotation files.
Works with the TSV format (columns: Sequence_Id, Type, Start, Stop, Strand, Locus_Tag, Gene, Product, DbXrefs).
"""
import os
import csv
import json
import re

BAKTA_DIR = "/home/ubuntu/plasmid_analysis/04_AMR_Plasmid_Association/Bakta_Annotations"
ENUM_FILE = "/home/ubuntu/plasmid_analysis/01_plasmid_enumeration/plasmid_enumeration_full.json"
OUTPUT_DIR = "/home/ubuntu/plasmid_analysis/04_AMR_Plasmid_Association"

# Known AMR gene patterns
AMR_PATTERNS = {
    'Beta-lactamase': [
        r'bla[A-Z]', r'beta-lactamase', r'carbapenemase', r'penicillin',
        r'cephalosporin', r'metallo-beta'
    ],
    'Aminoglycoside': [
        r'aac\(', r'aph\(', r'ant\(', r'armA', r'rmtB', r'rmtC', r'rmtF',
        r'aminoglycoside', r'16S rRNA methyltransferase'
    ],
    'Fluoroquinolone': [r'qnr[A-Z]', r'quinolone', r'oqx[AB]'],
    'Tetracycline': [r'tet\(', r'tet[A-Z]\b', r'tetracycline'],
    'Sulfonamide': [r'sul[12]', r'sulfonamide', r'dihydropteroate'],
    'Trimethoprim': [r'dfr[A-Z]', r'trimethoprim', r'dihydrofolate'],
    'Macrolide': [r'mph\(', r'msr\(', r'ere[A-Z]', r'macrolide', r'erm\('],
    'Fosfomycin': [r'fosA', r'fosfomycin'],
    'Chloramphenicol': [r'cat[A-Z]', r'chloramphenicol acetyltransferase', r'floR'],
    'Rifamycin': [r'arr-', r'rifampin'],
    'Colistin': [r'mcr-', r'colistin'],
    'Bleomycin': [r'\bble\b', r'bleomycin'],
    'Efflux': [r'emrD', r'eefA', r'eefB', r'eefC', r'acrA', r'acrB', r'mdtK', r'kdeA'],
    'Heavy_metal': [r'terC', r'terD', r'terE', r'terA', r'terB', r'merA', r'merR', r'arsA', r'arsB', r'arsC', r'pcoA', r'silA'],
    'Virulence': [r'rmpA', r'rmpC', r'iuc[A-D]', r'iutA', r'iro[B-E]', r'ybt[A-Z]', r'peg-344'],
}

# IS element / MGE patterns
MGE_PATTERNS = [
    r'IS\d+', r'IS[A-Z][a-z]+\d*', r'transpos', r'insertion sequence',
    r'integrase', r'recombinase', r'resolvase', r'conjugal transfer',
    r'mobilization', r'relaxase', r'type IV', r'tra[A-Z]\b', r'trb[A-Z]\b',
]

def classify_amr(gene_name, product):
    """Classify an AMR gene by resistance class."""
    text = f"{gene_name} {product}"
    for amr_class, patterns in AMR_PATTERNS.items():
        for pattern in patterns:
            if re.search(pattern, text, re.IGNORECASE):
                return amr_class
    return None

def is_mge(gene_name, product):
    """Check if a gene is an MGE element."""
    text = f"{gene_name} {product}"
    for pattern in MGE_PATTERNS:
        if re.search(pattern, text, re.IGNORECASE):
            return True
    return False

def get_contig_info(sample_name, enum_data):
    """Get contig sizes and classify as chromosome/plasmid."""
    contig_map = {}
    # Try to find this sample in enumeration data
    for isolate in enum_data:
        if sample_name in isolate.get('isolate', ''):
            for i, contig in enumerate(isolate.get('contigs', [])):
                contig_id = f"contig_{i+1}"
                size = contig.get('size', 0)
                contig_map[contig_id] = {
                    'size': size,
                    'type': 'chromosome' if size > 1000000 else 'plasmid',
                    'replicon': contig.get('replicon', 'unknown'),
                }
            break
    return contig_map

def parse_bakta_tsv(tsv_file, sample_name, contig_map):
    """Parse a Bakta TSV file and extract AMR genes and MGE elements."""
    amr_genes = []
    mge_elements = []
    
    with open(tsv_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            
            seq_id, feat_type, start, stop, strand, locus_tag, gene, product = parts[:8]
            dbxrefs = parts[8] if len(parts) > 8 else ''
            
            # Get contig info
            contig_info = contig_map.get(seq_id, {'size': 0, 'type': 'unknown', 'replicon': 'unknown'})
            
            # Check for AMR
            amr_class = classify_amr(gene, product)
            if amr_class:
                amr_genes.append({
                    'sample': sample_name,
                    'contig': seq_id,
                    'contig_type': contig_info['type'],
                    'contig_size': contig_info['size'],
                    'start': int(start),
                    'stop': int(stop),
                    'strand': strand,
                    'gene': gene,
                    'product': product,
                    'amr_class': amr_class,
                    'dbxrefs': dbxrefs,
                })
            
            # Check for MGE
            if is_mge(gene, product):
                mge_elements.append({
                    'sample': sample_name,
                    'contig': seq_id,
                    'contig_type': contig_info['type'],
                    'contig_size': contig_info['size'],
                    'start': int(start),
                    'stop': int(stop),
                    'strand': strand,
                    'gene': gene,
                    'product': product,
                    'dbxrefs': dbxrefs,
                })
    
    return amr_genes, mge_elements

def main():
    # Load plasmid enumeration data for contig classification
    with open(ENUM_FILE) as f:
        enum_raw = json.load(f)
    enum_data = enum_raw.get('isolates', []) if isinstance(enum_raw, dict) else enum_raw
    
    all_amr = []
    all_mge = []
    summary = {}
    
    # Process each sample with a TSV file
    samples_with_tsv = []
    for sample_dir in sorted(os.listdir(BAKTA_DIR)):
        tsv_file = os.path.join(BAKTA_DIR, sample_dir, f"{sample_dir}.tsv")
        if os.path.exists(tsv_file):
            samples_with_tsv.append(sample_dir)
    
    print(f"Found {len(samples_with_tsv)} samples with Bakta TSV files: {', '.join(samples_with_tsv)}")
    
    for sample in samples_with_tsv:
        tsv_file = os.path.join(BAKTA_DIR, sample, f"{sample}.tsv")
        
        # Build contig map from enumeration data
        # Map sample names to enumeration isolate names
        name_map = {
            'JKPB244': '1_JKPB244', 'JKPB284': '2_JKPB284',
            'JKPB381': '4_JKPB381', 'JKPR282': '5_JKPR282',
            'KP21802': '6_KP21802', 'KP21915': '7_KP21915',
            'MDRKP121': '8_MDRKP121', 'MDRKP122': '9_MDRKP122',
        }
        
        # Build contig map from the FASTA file
        fasta_name = name_map.get(sample, sample)
        contig_map = {}
        
        # Get contig sizes from the TSV itself
        contig_sizes = {}
        with open(tsv_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    seq_id = parts[0]
                    stop = int(parts[3])
                    if seq_id not in contig_sizes or stop > contig_sizes[seq_id]:
                        contig_sizes[seq_id] = stop
        
        for seq_id, max_pos in contig_sizes.items():
            contig_map[seq_id] = {
                'size': max_pos,
                'type': 'chromosome' if max_pos > 1000000 else 'plasmid',
                'replicon': 'unknown',
            }
        
        # Also try to get actual sizes from enumeration
        for isolate in enum_data:
            iso_name = isolate.get('isolate', '') if isinstance(isolate, dict) else ''
            if fasta_name in iso_name:
                # Use plasmid_details to identify plasmid contigs
                plasmid_contigs = set()
                for pd in isolate.get('plasmid_details', []):
                    ch = pd.get('contig_header', '')
                    plasmid_contigs.add(f"contig_{ch}")
                    cid = f"contig_{ch}"
                    if cid in contig_map:
                        contig_map[cid]['size'] = pd.get('size_bp', contig_map[cid]['size'])
                        contig_map[cid]['type'] = 'plasmid'
                # Mark chromosome
                for cid in contig_map:
                    if cid not in plasmid_contigs:
                        if contig_map[cid]['size'] > 1000000 or contig_map[cid].get('type') == 'unknown':
                            contig_map[cid]['type'] = 'chromosome' if contig_map[cid]['size'] > 500000 else 'plasmid'
                break
        
        amr_genes, mge_elements = parse_bakta_tsv(tsv_file, sample, contig_map)
        all_amr.extend(amr_genes)
        all_mge.extend(mge_elements)
        
        # Summary
        plasmid_amr = [g for g in amr_genes if g['contig_type'] == 'plasmid']
        chromo_amr = [g for g in amr_genes if g['contig_type'] == 'chromosome']
        
        summary[sample] = {
            'total_amr': len(amr_genes),
            'plasmid_amr': len(plasmid_amr),
            'chromosome_amr': len(chromo_amr),
            'total_mge': len(mge_elements),
            'contigs': {k: v for k, v in contig_map.items()},
            'key_amr_genes': list(set(g['gene'] for g in amr_genes if g['gene'])),
        }
        
        print(f"\n  {sample}:")
        print(f"    AMR genes: {len(amr_genes)} (plasmid: {len(plasmid_amr)}, chromosome: {len(chromo_amr)})")
        print(f"    MGE elements: {len(mge_elements)}")
        print(f"    Key genes: {', '.join(sorted(set(g['gene'] for g in amr_genes if g['gene']))[:15])}")
    
    # Save results
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # AMR genes CSV
    amr_csv = os.path.join(OUTPUT_DIR, "amr_genes_all_samples.csv")
    with open(amr_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'sample', 'contig', 'contig_type', 'contig_size', 'start', 'stop',
            'strand', 'gene', 'product', 'amr_class', 'dbxrefs'
        ])
        writer.writeheader()
        writer.writerows(all_amr)
    
    # MGE elements CSV
    mge_csv = os.path.join(OUTPUT_DIR, "mge_elements_all_samples.csv")
    with open(mge_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'sample', 'contig', 'contig_type', 'contig_size', 'start', 'stop',
            'strand', 'gene', 'product', 'dbxrefs'
        ])
        writer.writeheader()
        writer.writerows(all_mge)
    
    # Summary JSON
    summary_json = os.path.join(OUTPUT_DIR, "annotation_summary.json")
    with open(summary_json, 'w') as f:
        json.dump({
            'total_samples': len(samples_with_tsv),
            'samples_processed': samples_with_tsv,
            'total_amr_genes': len(all_amr),
            'total_mge_elements': len(all_mge),
            'per_sample': summary,
        }, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"TOTAL: {len(all_amr)} AMR genes, {len(all_mge)} MGE elements across {len(samples_with_tsv)} samples")
    print(f"Files saved to: {OUTPUT_DIR}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
