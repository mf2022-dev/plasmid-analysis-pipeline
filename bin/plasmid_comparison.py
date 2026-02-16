#!/usr/bin/env python3
"""
Comprehensive whole-plasmid comparison:
1. Pairwise alignment using nucmer/dnadiff between our 8 plasmids and 10 global references
2. Calculate ANI, SNP counts, structural differences
3. Generate pairwise identity matrix
4. Publication-quality visualizations
"""

import os
import subprocess
import re
import json
import numpy as np
import pandas as pd
from itertools import combinations

WORKDIR = "/home/ubuntu/plasmid_analysis/12_Plasmid_Comparison"
OUR_DIR = os.path.join(WORKDIR, "our_plasmids")
REF_DIR = os.path.join(WORKDIR, "reference_plasmids")
DNADIFF_DIR = os.path.join(WORKDIR, "dnadiff_results")
os.makedirs(DNADIFF_DIR, exist_ok=True)

def get_plasmid_files():
    """Get all plasmid FASTA files with short labels."""
    plasmids = {}
    
    # Our plasmids
    for f in sorted(os.listdir(OUR_DIR)):
        if f.endswith('.fasta'):
            # Extract isolate name
            parts = f.split('_')
            label = f"SA_{parts[1]}"  # SA = Saudi Arabia
            plasmids[label] = {
                'file': os.path.join(OUR_DIR, f),
                'group': 'This study',
                'origin': 'Saudi Arabia'
            }
    
    # Reference plasmids
    for f in sorted(os.listdir(REF_DIR)):
        if f.endswith('.fasta'):
            acc = f.split('_')[0]
            # Get origin from filename
            if 'Egypt' in f:
                origin = 'Egypt'
            elif 'UK' in f:
                origin = 'UK'
            else:
                origin = 'China'
            
            # Short label
            label = f"Ref_{acc}"
            plasmids[label] = {
                'file': os.path.join(REF_DIR, f),
                'group': 'Reference',
                'origin': origin
            }
    
    return plasmids

def run_dnadiff(ref_file, query_file, prefix):
    """Run dnadiff and parse results."""
    try:
        # Run nucmer first
        nucmer_cmd = f"nucmer --prefix={prefix} {ref_file} {query_file}"
        subprocess.run(nucmer_cmd, shell=True, capture_output=True, timeout=120)
        
        # Run dnadiff
        dnadiff_cmd = f"dnadiff -p {prefix} -d {prefix}.delta"
        subprocess.run(dnadiff_cmd, shell=True, capture_output=True, timeout=120)
        
        # Parse report
        report_file = f"{prefix}.report"
        if os.path.exists(report_file):
            return parse_dnadiff_report(report_file)
    except Exception as e:
        print(f"    Error: {e}")
    
    return None

def parse_dnadiff_report(report_file):
    """Parse dnadiff report for key metrics."""
    result = {}
    with open(report_file) as f:
        for line in f:
            line = line.strip()
            
            # Average identity
            if line.startswith('AvgIdentity'):
                parts = line.split()
                if len(parts) >= 3:
                    result['avg_identity_ref'] = float(parts[1])
                    result['avg_identity_qry'] = float(parts[2])
            
            # Aligned bases
            if line.startswith('AlignedBases'):
                parts = line.split()
                if len(parts) >= 3:
                    # Extract percentage from format like "128563(100.00%)"
                    match1 = re.search(r'(\d+)\(([0-9.]+)%\)', parts[1])
                    match2 = re.search(r'(\d+)\(([0-9.]+)%\)', parts[2])
                    if match1:
                        result['aligned_bases_ref'] = int(match1.group(1))
                        result['aligned_pct_ref'] = float(match1.group(2))
                    if match2:
                        result['aligned_bases_qry'] = int(match2.group(1))
                        result['aligned_pct_qry'] = float(match2.group(2))
            
            # Total SNPs
            if line.startswith('TotalSNPs'):
                parts = line.split()
                if len(parts) >= 3:
                    result['total_snps'] = int(parts[1])
            
            # Total Indels
            if line.startswith('TotalIndels'):
                parts = line.split()
                if len(parts) >= 3:
                    result['total_indels'] = int(parts[1])
            
            # Breakpoints
            if line.startswith('Breakpoints'):
                parts = line.split()
                if len(parts) >= 3:
                    result['breakpoints'] = int(parts[1])
            
            # Relocations
            if line.startswith('Relocations'):
                parts = line.split()
                if len(parts) >= 3:
                    result['relocations'] = int(parts[1])
            
            # Translocations
            if line.startswith('Translocations'):
                parts = line.split()
                if len(parts) >= 3:
                    result['translocations'] = int(parts[1])
            
            # Inversions
            if line.startswith('Inversions'):
                parts = line.split()
                if len(parts) >= 3:
                    result['inversions'] = int(parts[1])
    
    return result

def main():
    print("=" * 70)
    print("Whole-Plasmid Comparison: Our KPC Plasmids vs Global References")
    print("=" * 70)
    
    plasmids = get_plasmid_files()
    print(f"\nTotal plasmids: {len(plasmids)}")
    
    our_plasmids = {k: v for k, v in plasmids.items() if v['group'] == 'This study'}
    ref_plasmids = {k: v for k, v in plasmids.items() if v['group'] == 'Reference'}
    
    print(f"  Our plasmids: {len(our_plasmids)}")
    print(f"  Reference plasmids: {len(ref_plasmids)}")
    
    for label, info in sorted(plasmids.items()):
        # Get file size
        seq_len = 0
        with open(info['file']) as f:
            for line in f:
                if not line.startswith('>'):
                    seq_len += len(line.strip())
        print(f"  {label:25s}  {seq_len:>10,} bp  ({info['origin']})")
    
    # Run all pairwise comparisons
    all_labels = sorted(plasmids.keys())
    n = len(all_labels)
    
    # Initialize matrices
    identity_matrix = np.full((n, n), 100.0)
    snp_matrix = np.zeros((n, n), dtype=int)
    aligned_matrix = np.full((n, n), 100.0)
    
    results = []
    total_pairs = n * (n - 1) // 2
    pair_count = 0
    
    print(f"\nRunning {total_pairs} pairwise dnadiff comparisons...")
    
    for i in range(n):
        for j in range(i + 1, n):
            pair_count += 1
            label_i = all_labels[i]
            label_j = all_labels[j]
            
            if pair_count % 20 == 0:
                print(f"  Progress: {pair_count}/{total_pairs}...")
            
            prefix = os.path.join(DNADIFF_DIR, f"{label_i}_vs_{label_j}")
            
            result = run_dnadiff(
                plasmids[label_i]['file'],
                plasmids[label_j]['file'],
                prefix
            )
            
            if result:
                avg_id = result.get('avg_identity_ref', 0)
                snps = result.get('total_snps', 0)
                aligned_pct = result.get('aligned_pct_ref', 0)
                
                identity_matrix[i, j] = avg_id
                identity_matrix[j, i] = avg_id
                snp_matrix[i, j] = snps
                snp_matrix[j, i] = snps
                aligned_matrix[i, j] = aligned_pct
                aligned_matrix[j, i] = aligned_pct
                
                results.append({
                    'plasmid_1': label_i,
                    'plasmid_2': label_j,
                    'group_1': plasmids[label_i]['group'],
                    'group_2': plasmids[label_j]['group'],
                    'origin_1': plasmids[label_i]['origin'],
                    'origin_2': plasmids[label_j]['origin'],
                    'avg_identity': avg_id,
                    'total_snps': snps,
                    'total_indels': result.get('total_indels', 0),
                    'aligned_pct_ref': result.get('aligned_pct_ref', 0),
                    'aligned_pct_qry': result.get('aligned_pct_qry', 0),
                    'breakpoints': result.get('breakpoints', 0),
                    'relocations': result.get('relocations', 0),
                    'translocations': result.get('translocations', 0),
                    'inversions': result.get('inversions', 0),
                })
    
    print(f"\nCompleted {len(results)} pairwise comparisons")
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(WORKDIR, 'pairwise_comparison_results.csv'), index=False)
    
    # Save matrices
    id_df = pd.DataFrame(identity_matrix, index=all_labels, columns=all_labels)
    id_df.to_csv(os.path.join(WORKDIR, 'identity_matrix.csv'))
    
    snp_df = pd.DataFrame(snp_matrix, index=all_labels, columns=all_labels)
    snp_df.to_csv(os.path.join(WORKDIR, 'snp_matrix.csv'))
    
    aligned_df = pd.DataFrame(aligned_matrix, index=all_labels, columns=all_labels)
    aligned_df.to_csv(os.path.join(WORKDIR, 'aligned_coverage_matrix.csv'))
    
    # Print summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    # Our plasmids vs each other
    our_pairs = results_df[(results_df['group_1'] == 'This study') & (results_df['group_2'] == 'This study')]
    if len(our_pairs) > 0:
        print(f"\nOUR PLASMIDS (intra-study):")
        print(f"  Pairwise identity: {our_pairs['avg_identity'].min():.4f}% - {our_pairs['avg_identity'].max():.4f}%")
        print(f"  Mean identity: {our_pairs['avg_identity'].mean():.4f}%")
        print(f"  SNP range: {our_pairs['total_snps'].min()} - {our_pairs['total_snps'].max()}")
        print(f"  Mean SNPs: {our_pairs['total_snps'].mean():.1f}")
        print(f"  Aligned coverage: {our_pairs['aligned_pct_ref'].min():.2f}% - {our_pairs['aligned_pct_ref'].max():.2f}%")
    
    # Our plasmids vs references
    cross_pairs = results_df[
        ((results_df['group_1'] == 'This study') & (results_df['group_2'] == 'Reference')) |
        ((results_df['group_1'] == 'Reference') & (results_df['group_2'] == 'This study'))
    ]
    if len(cross_pairs) > 0:
        print(f"\nOUR vs GLOBAL REFERENCES:")
        print(f"  Pairwise identity: {cross_pairs['avg_identity'].min():.4f}% - {cross_pairs['avg_identity'].max():.4f}%")
        print(f"  Mean identity: {cross_pairs['avg_identity'].mean():.4f}%")
        print(f"  SNP range: {cross_pairs['total_snps'].min()} - {cross_pairs['total_snps'].max()}")
        print(f"  Mean SNPs: {cross_pairs['total_snps'].mean():.1f}")
        print(f"  Aligned coverage: {cross_pairs['aligned_pct_ref'].min():.2f}% - {cross_pairs['aligned_pct_ref'].max():.2f}%")
    
    # References vs each other
    ref_pairs = results_df[(results_df['group_1'] == 'Reference') & (results_df['group_2'] == 'Reference')]
    if len(ref_pairs) > 0:
        print(f"\nGLOBAL REFERENCES (inter-reference):")
        print(f"  Pairwise identity: {ref_pairs['avg_identity'].min():.4f}% - {ref_pairs['avg_identity'].max():.4f}%")
        print(f"  Mean identity: {ref_pairs['avg_identity'].mean():.4f}%")
        print(f"  SNP range: {ref_pairs['total_snps'].min()} - {ref_pairs['total_snps'].max()}")
        print(f"  Mean SNPs: {ref_pairs['total_snps'].mean():.1f}")
    
    # Save summary
    summary = {
        'total_plasmids': len(plasmids),
        'our_plasmids': len(our_plasmids),
        'reference_plasmids': len(ref_plasmids),
        'total_pairwise_comparisons': len(results),
        'intra_study': {
            'n_pairs': len(our_pairs),
            'identity_range': [float(our_pairs['avg_identity'].min()), float(our_pairs['avg_identity'].max())] if len(our_pairs) > 0 else [],
            'mean_identity': float(our_pairs['avg_identity'].mean()) if len(our_pairs) > 0 else 0,
            'snp_range': [int(our_pairs['total_snps'].min()), int(our_pairs['total_snps'].max())] if len(our_pairs) > 0 else [],
            'mean_snps': float(our_pairs['total_snps'].mean()) if len(our_pairs) > 0 else 0,
        },
        'study_vs_reference': {
            'n_pairs': len(cross_pairs),
            'identity_range': [float(cross_pairs['avg_identity'].min()), float(cross_pairs['avg_identity'].max())] if len(cross_pairs) > 0 else [],
            'mean_identity': float(cross_pairs['avg_identity'].mean()) if len(cross_pairs) > 0 else 0,
            'snp_range': [int(cross_pairs['total_snps'].min()), int(cross_pairs['total_snps'].max())] if len(cross_pairs) > 0 else [],
            'mean_snps': float(cross_pairs['total_snps'].mean()) if len(cross_pairs) > 0 else 0,
        },
        'inter_reference': {
            'n_pairs': len(ref_pairs),
            'identity_range': [float(ref_pairs['avg_identity'].min()), float(ref_pairs['avg_identity'].max())] if len(ref_pairs) > 0 else [],
            'mean_identity': float(ref_pairs['avg_identity'].mean()) if len(ref_pairs) > 0 else 0,
            'snp_range': [int(ref_pairs['total_snps'].min()), int(ref_pairs['total_snps'].max())] if len(ref_pairs) > 0 else [],
            'mean_snps': float(ref_pairs['total_snps'].mean()) if len(ref_pairs) > 0 else 0,
        }
    }
    
    with open(os.path.join(WORKDIR, 'comparison_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nAll results saved to {WORKDIR}/")

if __name__ == '__main__':
    main()
