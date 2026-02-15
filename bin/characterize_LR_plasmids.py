#!/usr/bin/env python3
"""
Comprehensive Plasmid Characterization from Long-Read Complete Genomes
Step 1: Enumerate all plasmids, determine sizes, circularity, and basic features
"""

import os
import json
import csv

def parse_fasta(filepath):
    """Parse FASTA file and return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))
    return sequences

def check_circularity(seq, window=500):
    """Check if sequence is likely circular by comparing start and end."""
    if len(seq) < window * 2:
        return False, 0
    start = seq[:window].upper()
    end = seq[-window:].upper()
    # Simple overlap check
    matches = sum(1 for a, b in zip(start, end) if a == b)
    identity = matches / window * 100
    return identity > 95, round(identity, 1)

def gc_content(seq):
    """Calculate GC content."""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    total = len(seq) - seq.count('N')
    return round(gc / total * 100, 1) if total > 0 else 0

def classify_contig(length, header):
    """Classify contig as chromosome, plasmid, or small element."""
    header_lower = header.lower()
    if length > 4000000:
        return 'Chromosome'
    elif length >= 100000:
        return 'Large Plasmid'
    elif length >= 10000:
        return 'Medium Plasmid'
    elif length >= 2000:
        return 'Small Plasmid'
    else:
        return 'Small Fragment'

# Main analysis
LR_DIR = '/home/ubuntu/plasmid_analysis/LR_assemblies'
OUTPUT_DIR = '/home/ubuntu/plasmid_analysis/01_plasmid_enumeration'
os.makedirs(OUTPUT_DIR, exist_ok=True)

all_isolates = []
all_plasmids = []
plasmid_id = 0

print("=" * 100)
print("STEP 1: COMPLETE PLASMID ENUMERATION FROM LONG-READ ASSEMBLIES")
print("=" * 100)

for fname in sorted(os.listdir(LR_DIR)):
    if not (fname.endswith('.fasta') or fname.endswith('.fa')):
        continue
    
    filepath = os.path.join(LR_DIR, fname)
    sequences = parse_fasta(filepath)
    
    # Extract isolate name
    isolate_name = fname.replace('.fasta', '').replace('_LR', '').replace('_LR*', '').rstrip('*')
    
    print(f"\n{'━' * 100}")
    print(f"ISOLATE: {isolate_name} ({fname})")
    print(f"{'━' * 100}")
    
    isolate_data = {
        'isolate': isolate_name,
        'filename': fname,
        'total_contigs': len(sequences),
        'total_size_bp': sum(len(s) for _, s in sequences),
        'chromosomes': 0,
        'plasmids': 0,
        'plasmid_details': []
    }
    
    # Sort by size (largest first)
    sequences_sorted = sorted(sequences, key=lambda x: len(x[1]), reverse=True)
    
    for i, (header, seq) in enumerate(sequences_sorted):
        seqlen = len(seq)
        gc = gc_content(seq)
        classification = classify_contig(seqlen, header)
        
        # Check circularity from header
        header_lower = header.lower()
        is_circular_header = 'circular' in header_lower
        
        # Also check sequence-based circularity
        is_circular_seq, circ_identity = check_circularity(seq)
        is_circular = is_circular_header or is_circular_seq
        
        topology = 'Circular' if is_circular else 'Linear'
        
        if classification == 'Chromosome':
            isolate_data['chromosomes'] += 1
            print(f"  CHROMOSOME: {seqlen:,} bp ({seqlen/1e6:.2f} Mb) | GC: {gc}% | {topology}")
        else:
            isolate_data['plasmids'] += 1
            plasmid_id += 1
            
            plasmid_info = {
                'plasmid_id': f'p{plasmid_id:03d}',
                'isolate': isolate_name,
                'contig_header': header.split()[0] if header else f'contig_{i+1}',
                'size_bp': seqlen,
                'size_kb': round(seqlen / 1000, 1),
                'gc_content': gc,
                'topology': topology,
                'classification': classification,
                'circular_identity': circ_identity if is_circular_seq else 'N/A'
            }
            
            isolate_data['plasmid_details'].append(plasmid_info)
            all_plasmids.append(plasmid_info)
            
            print(f"  PLASMID {plasmid_info['plasmid_id']}: {seqlen:,} bp ({plasmid_info['size_kb']} kb) | "
                  f"GC: {gc}% | {topology} | {classification}")
    
    # Also extract plasmid sequences to individual files
    plasmid_seq_dir = os.path.join(OUTPUT_DIR, 'plasmid_sequences', isolate_name)
    os.makedirs(plasmid_seq_dir, exist_ok=True)
    
    plasmid_idx = 0
    for header, seq in sequences_sorted:
        seqlen = len(seq)
        if seqlen < 4000000:  # Not chromosome
            plasmid_idx += 1
            out_file = os.path.join(plasmid_seq_dir, f'{isolate_name}_plasmid_{plasmid_idx}_{seqlen}bp.fasta')
            with open(out_file, 'w') as f:
                f.write(f'>{isolate_name}_plasmid_{plasmid_idx} length={seqlen}\n')
                # Write sequence in 80-char lines
                for j in range(0, len(seq), 80):
                    f.write(seq[j:j+80] + '\n')
    
    all_isolates.append(isolate_data)
    
    print(f"  SUMMARY: {isolate_data['chromosomes']} chromosome(s), {isolate_data['plasmids']} plasmid(s)")
    print(f"  Plasmid sequences extracted to: {plasmid_seq_dir}")

# Overall summary
print(f"\n{'=' * 100}")
print("OVERALL SUMMARY")
print(f"{'=' * 100}")
print(f"  Total isolates analyzed:     {len(all_isolates)}")
print(f"  Total plasmids identified:   {len(all_plasmids)}")
print(f"  Average plasmids/isolate:    {len(all_plasmids)/len(all_isolates):.1f}")

# Size distribution
sizes = [p['size_kb'] for p in all_plasmids]
print(f"\n  Plasmid size range:          {min(sizes):.1f} kb - {max(sizes):.1f} kb")
print(f"  Mean plasmid size:           {sum(sizes)/len(sizes):.1f} kb")

# Topology
circular = sum(1 for p in all_plasmids if p['topology'] == 'Circular')
print(f"\n  Circular plasmids:           {circular}/{len(all_plasmids)} ({circular/len(all_plasmids)*100:.1f}%)")

# Classification distribution
from collections import Counter
class_counts = Counter(p['classification'] for p in all_plasmids)
print(f"\n  Classification:")
for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
    print(f"    {cls}: {count}")

# Save results
# CSV summary
csv_file = os.path.join(OUTPUT_DIR, 'plasmid_enumeration_summary.csv')
with open(csv_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['plasmid_id', 'isolate', 'contig_header', 'size_bp', 'size_kb', 
                                            'gc_content', 'topology', 'classification', 'circular_identity'])
    writer.writeheader()
    writer.writerows(all_plasmids)

# JSON with full details
json_file = os.path.join(OUTPUT_DIR, 'plasmid_enumeration_full.json')
with open(json_file, 'w') as f:
    json.dump({
        'isolates': all_isolates,
        'plasmids': all_plasmids,
        'summary': {
            'total_isolates': len(all_isolates),
            'total_plasmids': len(all_plasmids),
            'avg_plasmids_per_isolate': round(len(all_plasmids)/len(all_isolates), 1),
            'size_range_kb': [min(sizes), max(sizes)],
            'circular_count': circular,
            'classification_counts': dict(class_counts)
        }
    }, f, indent=2)

# Isolate-level summary
isolate_csv = os.path.join(OUTPUT_DIR, 'isolate_plasmid_summary.csv')
with open(isolate_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Isolate', 'Total_Contigs', 'Chromosomes', 'Plasmids', 'Plasmid_Sizes_kb', 'Total_Plasmid_Size_kb'])
    for iso in all_isolates:
        plasmid_sizes = [str(p['size_kb']) for p in iso['plasmid_details']]
        total_plasmid_size = sum(p['size_kb'] for p in iso['plasmid_details'])
        writer.writerow([
            iso['isolate'], iso['total_contigs'], iso['chromosomes'], iso['plasmids'],
            '; '.join(plasmid_sizes), round(total_plasmid_size, 1)
        ])

print(f"\n  Results saved to:")
print(f"    {csv_file}")
print(f"    {json_file}")
print(f"    {isolate_csv}")
print(f"    Plasmid sequences: {os.path.join(OUTPUT_DIR, 'plasmid_sequences')}/")
