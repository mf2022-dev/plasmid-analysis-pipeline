#!/usr/bin/env python3
"""
Assembly Quality Assessment for Plasmid Analysis
Determines sequencing type (LR vs SR), counts contigs, identifies potential plasmids,
and assesses whether complete circular plasmids are present.
"""

import os
import sys
import json
from collections import defaultdict

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

def assess_assembly(filepath):
    """Comprehensive assembly quality assessment."""
    filename = os.path.basename(filepath)
    sequences = parse_fasta(filepath)
    
    # Basic stats
    num_contigs = len(sequences)
    lengths = sorted([len(seq) for _, seq in sequences], reverse=True)
    total_size = sum(lengths)
    
    # N50 calculation
    cumsum = 0
    n50 = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total_size / 2:
            n50 = l
            break
    
    # Largest and smallest contigs
    largest = lengths[0] if lengths else 0
    smallest = lengths[-1] if lengths else 0
    
    # Determine sequencing type from filename
    if '_LR' in filename:
        seq_type = 'Long-Read'
    elif '_SR' in filename:
        seq_type = 'Short-Read'
    else:
        # Infer from assembly characteristics
        if num_contigs <= 10 and n50 > 1000000:
            seq_type = 'Long-Read (inferred)'
        elif num_contigs > 100:
            seq_type = 'Short-Read (inferred)'
        else:
            seq_type = 'Unknown'
    
    # Identify potential plasmid contigs
    # Chromosome: typically >4 Mb for K. pneumoniae
    # Plasmids: typically 2 kb - 400 kb
    chromosome_contigs = []
    plasmid_candidates = []
    small_fragments = []
    
    for header, seq in sequences:
        seqlen = len(seq)
        if seqlen > 4000000:  # >4 Mb = chromosome
            chromosome_contigs.append((header, seqlen))
        elif seqlen >= 2000:  # 2 kb - 4 Mb = potential plasmid
            plasmid_candidates.append((header, seqlen))
        else:  # <2 kb = small fragment
            small_fragments.append((header, seqlen))
    
    # Check for circular indicators in headers
    circular_contigs = []
    for header, seq in sequences:
        header_lower = header.lower()
        if 'circular' in header_lower or 'plasmid' in header_lower or 'complete' in header_lower:
            circular_contigs.append((header, len(seq)))
    
    # GC content
    total_gc = 0
    total_bases = 0
    for _, seq in sequences:
        seq_upper = seq.upper()
        total_gc += seq_upper.count('G') + seq_upper.count('C')
        total_bases += len(seq_upper) - seq_upper.count('N')
    gc_content = (total_gc / total_bases * 100) if total_bases > 0 else 0
    
    # Assembly quality verdict
    if seq_type.startswith('Long-Read'):
        if num_contigs <= 10 and largest > 5000000:
            quality = 'COMPLETE GENOME - Suitable for full plasmid analysis'
        elif num_contigs <= 20:
            quality = 'NEAR-COMPLETE - Good for plasmid analysis'
        else:
            quality = 'FRAGMENTED LR - May need reassembly'
    else:
        if num_contigs <= 50:
            quality = 'HIGH-QUALITY DRAFT - Replicon typing possible, complete plasmids unlikely'
        elif num_contigs <= 200:
            quality = 'DRAFT - Replicon typing possible, plasmid structure fragmented'
        else:
            quality = 'FRAGMENTED - Limited plasmid analysis possible'
    
    return {
        'filename': filename,
        'sequencing_type': seq_type,
        'num_contigs': num_contigs,
        'total_size_bp': total_size,
        'total_size_mb': round(total_size / 1e6, 2),
        'n50': n50,
        'n50_kb': round(n50 / 1000, 1),
        'largest_contig': largest,
        'largest_contig_mb': round(largest / 1e6, 2),
        'smallest_contig': smallest,
        'gc_content': round(gc_content, 1),
        'chromosome_contigs': len(chromosome_contigs),
        'plasmid_candidates': len(plasmid_candidates),
        'small_fragments': len(small_fragments),
        'circular_contigs': len(circular_contigs),
        'quality_verdict': quality,
        'plasmid_candidate_sizes': [(h.split()[0], round(s/1000, 1)) for h, s in sorted(plasmid_candidates, key=lambda x: x[1], reverse=True)[:10]],
        'chromosome_sizes': [(h.split()[0], round(s/1e6, 2)) for h, s in chromosome_contigs]
    }

# Main execution
fasta_dir = '/home/ubuntu/plasmid_analysis/sample_fastas'
results = []

print("=" * 80)
print("ASSEMBLY QUALITY ASSESSMENT FOR PLASMID ANALYSIS")
print("=" * 80)

for fname in sorted(os.listdir(fasta_dir)):
    if fname.endswith('.fasta') or fname.endswith('.fa'):
        filepath = os.path.join(fasta_dir, fname)
        result = assess_assembly(filepath)
        results.append(result)
        
        print(f"\n{'─' * 80}")
        print(f"FILE: {result['filename']}")
        print(f"{'─' * 80}")
        print(f"  Sequencing Type:    {result['sequencing_type']}")
        print(f"  Total Contigs:      {result['num_contigs']}")
        print(f"  Total Size:         {result['total_size_mb']} Mb")
        print(f"  N50:                {result['n50_kb']} kb")
        print(f"  Largest Contig:     {result['largest_contig_mb']} Mb")
        print(f"  GC Content:         {result['gc_content']}%")
        print(f"  Chromosome(s):      {result['chromosome_contigs']}")
        print(f"  Plasmid Candidates: {result['plasmid_candidates']}")
        print(f"  Circular Contigs:   {result['circular_contigs']}")
        print(f"  QUALITY:            {result['quality_verdict']}")
        
        if result['plasmid_candidate_sizes']:
            print(f"  Plasmid Sizes (kb): {', '.join([f'{name}={size}kb' for name, size in result['plasmid_candidate_sizes']])}")
        if result['chromosome_sizes']:
            print(f"  Chromosome(s):      {', '.join([f'{name}={size}Mb' for name, size in result['chromosome_sizes']])}")

# Summary
print(f"\n{'=' * 80}")
print("SUMMARY")
print(f"{'=' * 80}")
lr_count = sum(1 for r in results if 'Long' in r['sequencing_type'])
sr_count = sum(1 for r in results if 'Short' in r['sequencing_type'])
complete_count = sum(1 for r in results if 'COMPLETE' in r['quality_verdict'])

print(f"  Long-Read assemblies:  {lr_count}")
print(f"  Short-Read assemblies: {sr_count}")
print(f"  Complete genomes:      {complete_count}")
print(f"  Total plasmid candidates across all: {sum(r['plasmid_candidates'] for r in results)}")

# Save results
with open('/home/ubuntu/plasmid_analysis/assembly_assessment.json', 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to: /home/ubuntu/plasmid_analysis/assembly_assessment.json")
