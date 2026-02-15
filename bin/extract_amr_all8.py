#!/usr/bin/env python3
"""
Extract AMR genes, IS elements, and MGE from all 8 Bakta-annotated genomes.
Associates each gene with its contig (chromosome vs plasmid) using enumeration data.
"""

import json
import csv
import os
from pathlib import Path
from collections import defaultdict

BASE = Path("/home/ubuntu/plasmid_analysis")
BAKTA_DIR = BASE / "04_AMR_Plasmid_Association" / "Bakta_Annotations"
ENUM_FILE = BASE / "01_plasmid_enumeration" / "plasmid_enumeration_summary.csv"
OUTPUT_DIR = BASE / "04_AMR_Plasmid_Association"

# Known AMR gene patterns
AMR_KEYWORDS = [
    'bla', 'aac', 'aph', 'ant', 'arm', 'rmt', 'mph', 'msr', 'ere', 'erm',
    'tet', 'sul', 'dfr', 'qnr', 'oqx', 'fos', 'cat', 'cml', 'flo',
    'mcr', 'van', 'ble', 'str', 'aad', 'arr', 'lin', 'lnu', 'vat',
    'mef', 'nim', 'cfr', 'optrA', 'poxtA', 'ndm', 'kpc', 'oxa', 'ctx',
    'shv', 'tem', 'vim', 'imp', 'ges', 'per', 'veb', 'pse', 'dha',
    'cmy', 'fox', 'mox', 'act', 'mir', 'acc', 'qep', 'emr', 'ade',
    'mdt', 'acr', 'tol', 'mac', 'mex', 'opr', 'nar', 'mar', 'sox',
    'rob', 'ram', 'aar', 'bcr', 'ept', 'pmr', 'arn', 'pho', 'bas',
    'mgt', 'ire', 'iuc', 'iut', 'ybt', 'irp', 'fyu', 'sit', 'ent',
    'fep', 'fec', 'feo', 'hmu', 'chu', 'has', 'hem', 'hbp'
]

AMR_PRODUCTS = [
    'beta-lactamase', 'aminoglycoside', 'carbapenemase', 'resistance',
    'efflux', 'tetracycline', 'chloramphenicol', 'macrolide', 'quinolone',
    'sulfonamide', 'trimethoprim', 'fosfomycin', 'colistin', 'vancomycin',
    'bleomycin', 'streptomycin', 'rifampicin', 'lincosamide', 'oxacillinase',
    'metallo-beta-lactamase', 'extended-spectrum', 'ESBL', 'AmpC',
    'methyltransferase', '16S rRNA methyltransferase', 'phosphotransferase',
    'acetyltransferase', 'nucleotidyltransferase'
]

MGE_KEYWORDS = [
    'IS', 'transpos', 'integrase', 'recombinase', 'insertion sequence',
    'mobile element', 'conjugat', 'relaxase', 'mobilization', 'resolvase',
    'invertase', 'Tn', 'phage', 'plasmid'
]

# Load enumeration data to map contigs to chromosome/plasmid
def load_enumeration():
    """Load plasmid enumeration data."""
    enum_data = {}
    if ENUM_FILE.exists():
        with open(ENUM_FILE) as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample = row.get('Sample', row.get('sample', ''))
                # Clean sample name
                for prefix in ['1_', '2_', '3_', '4_', '5_', '6_', '7_', '8_', '9_']:
                    if sample.startswith(prefix):
                        sample = sample[2:]
                sample = sample.replace('_LR', '')
                enum_data[sample] = row
    return enum_data

def classify_contig(contig_name, contig_size, sample_name, enum_data):
    """Classify a contig as chromosome or plasmid based on size."""
    # Chromosome is typically >4 Mb for K. pneumoniae
    if contig_size > 4000000:
        return "Chromosome"
    elif contig_size > 300000:
        return f"Large plasmid (~{contig_size//1000}kb)"
    elif contig_size > 100000:
        return f"Medium plasmid (~{contig_size//1000}kb)"
    elif contig_size > 50000:
        return f"Small plasmid (~{contig_size//1000}kb)"
    else:
        return f"Small replicon (~{contig_size//1000}kb)"

def is_amr_gene(gene_name, product):
    """Check if a gene is an AMR gene."""
    gene_lower = (gene_name or '').lower()
    product_lower = (product or '').lower()
    
    for kw in AMR_KEYWORDS:
        if kw.lower() in gene_lower:
            return True
    for kw in AMR_PRODUCTS:
        if kw.lower() in product_lower:
            return True
    return False

def is_mge(gene_name, product):
    """Check if a gene is a mobile genetic element."""
    gene_lower = (gene_name or '').lower()
    product_lower = (product or '').lower()
    
    for kw in MGE_KEYWORDS:
        if kw.lower() in gene_lower or kw.lower() in product_lower:
            return True
    return False

def process_json(json_path, sample_name, enum_data):
    """Process a Bakta JSON file to extract AMR and MGE data."""
    amr_genes = []
    mge_elements = []
    
    with open(json_path) as f:
        data = json.load(f)
    
    # Build contig size map
    contig_sizes = {}
    for seq in data.get('sequences', []):
        contig_sizes[seq['id']] = seq.get('length', 0)
    
    # Process features
    for feat in data.get('features', []):
        contig = feat.get('contig', feat.get('sequence', ''))
        gene = feat.get('gene', '')
        product = feat.get('product', '')
        start = feat.get('start', 0)
        stop = feat.get('stop', 0)
        strand = feat.get('strand', '')
        feat_type = feat.get('type', '')
        db_xrefs = feat.get('db_xrefs', [])
        
        contig_size = contig_sizes.get(contig, 0)
        location = classify_contig(contig, contig_size, sample_name, enum_data)
        
        if is_amr_gene(gene, product):
            amr_genes.append({
                'Sample': sample_name,
                'Gene': gene,
                'Product': product,
                'Contig': contig,
                'Contig_Size': contig_size,
                'Location': location,
                'Start': start,
                'Stop': stop,
                'Strand': strand,
                'Type': feat_type,
                'DB_Xrefs': ';'.join(db_xrefs) if db_xrefs else ''
            })
        
        if is_mge(gene, product):
            mge_elements.append({
                'Sample': sample_name,
                'Gene': gene,
                'Product': product,
                'Contig': contig,
                'Contig_Size': contig_size,
                'Location': location,
                'Start': start,
                'Stop': stop,
                'Strand': strand,
                'Type': feat_type
            })
    
    return amr_genes, mge_elements

def process_tsv(tsv_path, sample_name, enum_data):
    """Process a Bakta TSV file to extract AMR and MGE data."""
    amr_genes = []
    mge_elements = []
    
    # First, try to get contig sizes from the FNA file
    fna_path = tsv_path.parent / f"{sample_name}.fna"
    contig_sizes = {}
    if fna_path.exists():
        current_contig = None
        current_size = 0
        with open(fna_path) as f:
            for line in f:
                if line.startswith('>'):
                    if current_contig:
                        contig_sizes[current_contig] = current_size
                    current_contig = line.strip().split()[0][1:]
                    current_size = 0
                else:
                    current_size += len(line.strip())
            if current_contig:
                contig_sizes[current_contig] = current_size
    
    with open(tsv_path) as f:
        # Skip comment lines (start with '# ') but keep header (starts with '#Sequence')
        lines = [l for l in f if not l.startswith('# ')]
    
    import io
    reader = csv.DictReader(io.StringIO(''.join(lines)), delimiter='\t')
    for row in reader:
        contig = row.get('#Sequence Id', row.get('Sequence Id', ''))
        gene = row.get('Gene', '') or ''
        product = row.get('Product', '') or ''
        start = row.get('Start', 0)
        stop = row.get('Stop', 0)
        strand = row.get('Strand', '')
        feat_type = row.get('Type', '')
        db_xrefs = row.get('DbXrefs', '')
        
        contig_size = contig_sizes.get(contig, 0)
        location = classify_contig(contig, contig_size, sample_name, enum_data)
        
        if is_amr_gene(gene, product):
            amr_genes.append({
                'Sample': sample_name,
                'Gene': gene,
                'Product': product,
                'Contig': contig,
                'Contig_Size': contig_size,
                'Location': location,
                'Start': start,
                'Stop': stop,
                'Strand': strand,
                'Type': feat_type,
                'DB_Xrefs': db_xrefs
            })
        
        if is_mge(gene, product):
            mge_elements.append({
                'Sample': sample_name,
                'Gene': gene,
                'Product': product,
                'Contig': contig,
                'Contig_Size': contig_size,
                'Location': location,
                'Start': start,
                'Stop': stop,
                'Strand': strand,
                'Type': feat_type
            })
    
    return amr_genes, mge_elements

def main():
    enum_data = load_enumeration()
    
    all_amr = []
    all_mge = []
    
    samples = ['JKPB244', 'JKPB284', 'JKPB381', 'JKPR282', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']
    
    for sample in samples:
        sample_dir = BAKTA_DIR / sample
        json_path = sample_dir / f"{sample}.json"
        tsv_path = sample_dir / f"{sample}.tsv"
        
        if json_path.exists():
            print(f"Processing {sample} (JSON)...")
            amr, mge = process_json(json_path, sample, enum_data)
        elif tsv_path.exists():
            print(f"Processing {sample} (TSV)...")
            amr, mge = process_tsv(tsv_path, sample, enum_data)
        else:
            print(f"WARNING: No annotation found for {sample}")
            continue
        
        all_amr.extend(amr)
        all_mge.extend(mge)
        print(f"  {sample}: {len(amr)} AMR genes, {len(mge)} MGE elements")
    
    # Write AMR CSV
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    amr_csv = OUTPUT_DIR / "amr_genes_all8_samples.csv"
    if all_amr:
        with open(amr_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=all_amr[0].keys())
            writer.writeheader()
            writer.writerows(all_amr)
        print(f"\nAMR genes saved: {amr_csv} ({len(all_amr)} records)")
    
    # Write MGE CSV
    mge_csv = OUTPUT_DIR / "mge_elements_all8_samples.csv"
    if all_mge:
        with open(mge_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=all_mge[0].keys())
            writer.writeheader()
            writer.writerows(all_mge)
        print(f"MGE elements saved: {mge_csv} ({len(all_mge)} records)")
    
    # Summary statistics
    print("\n" + "="*60)
    print("SUMMARY - ALL 8 GENOMES")
    print("="*60)
    
    # AMR by sample
    amr_by_sample = defaultdict(list)
    for a in all_amr:
        amr_by_sample[a['Sample']].append(a)
    
    print("\nAMR genes per sample:")
    for s in samples:
        genes = amr_by_sample.get(s, [])
        plasmid_genes = [g for g in genes if 'plasmid' in g['Location'].lower()]
        chrom_genes = [g for g in genes if 'Chromosome' in g['Location']]
        print(f"  {s}: {len(genes)} total ({len(plasmid_genes)} plasmid, {len(chrom_genes)} chromosome)")
    
    # Key carbapenemase distribution
    print("\nCarbapenemase distribution:")
    for gene_prefix in ['blaKPC', 'blaNDM', 'blaOXA']:
        carriers = set()
        for a in all_amr:
            if (a['Gene'] or '').startswith(gene_prefix):
                carriers.add(a['Sample'])
        if carriers:
            print(f"  {gene_prefix}: {len(carriers)}/8 isolates ({', '.join(sorted(carriers))})")
    
    # MGE by sample
    mge_by_sample = defaultdict(list)
    for m in all_mge:
        mge_by_sample[m['Sample']].append(m)
    
    print("\nMGE elements per sample:")
    for s in samples:
        elements = mge_by_sample.get(s, [])
        is_elements = [e for e in elements if 'IS' in (e.get('Gene', '') or '') or 'insertion' in (e.get('Product', '') or '').lower()]
        transposases = [e for e in elements if 'transpos' in (e.get('Product', '') or '').lower()]
        print(f"  {s}: {len(elements)} total ({len(is_elements)} IS, {len(transposases)} transposases)")

if __name__ == '__main__':
    main()
