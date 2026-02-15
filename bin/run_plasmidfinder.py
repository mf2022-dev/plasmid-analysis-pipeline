#!/usr/bin/env python3
"""
PlasmidFinder Replicon Typing - BLAST-based approach
Runs PlasmidFinder database against all complete LR assemblies
to identify plasmid replicon types (Inc groups)
"""

import subprocess
import os
import csv
import json
from collections import defaultdict

# Configuration
FASTA_DIR = "/home/ubuntu/plasmid_analysis/LR_assemblies"
DB_DIR = "/home/ubuntu/plasmid_analysis/plasmidfinder_db"
OUTPUT_DIR = "/home/ubuntu/plasmid_analysis/02_replicon_typing"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# PlasmidFinder databases - focus on Enterobacteriales for K. pneumoniae
# but also check all databases for comprehensive coverage
DB_FILES = {
    "enterobacteriales": os.path.join(DB_DIR, "enterobacteriales.fsa"),
    "Inc18": os.path.join(DB_DIR, "Inc18.fsa"),
    "NT_Rep": os.path.join(DB_DIR, "NT_Rep.fsa"),
    "Rep1": os.path.join(DB_DIR, "Rep1.fsa"),
    "Rep2": os.path.join(DB_DIR, "Rep2.fsa"),
    "Rep3": os.path.join(DB_DIR, "Rep3.fsa"),
    "RepA_N": os.path.join(DB_DIR, "RepA_N.fsa"),
    "RepL": os.path.join(DB_DIR, "RepL.fsa"),
    "Rep_trans": os.path.join(DB_DIR, "Rep_trans.fsa"),
}

# Thresholds (PlasmidFinder defaults)
MIN_IDENTITY = 80.0  # minimum % identity
MIN_COVERAGE = 60.0  # minimum % coverage

def get_fasta_files():
    """Get all LR FASTA files"""
    files = []
    for f in sorted(os.listdir(FASTA_DIR)):
        if f.endswith('.fasta') and '_LR' in f:
            files.append(os.path.join(FASTA_DIR, f))
    return files

def parse_fasta_contigs(fasta_file):
    """Parse FASTA file to get contig names and lengths"""
    contigs = {}
    current_name = None
    current_seq = []
    with open(fasta_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    contigs[current_name] = len(''.join(current_seq))
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_name:
        contigs[current_name] = len(''.join(current_seq))
    return contigs

def get_db_seq_lengths(db_fasta):
    """Get lengths of sequences in PlasmidFinder database"""
    lengths = {}
    current_name = None
    current_seq = []
    with open(db_fasta) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    lengths[current_name] = len(''.join(current_seq))
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_name:
        lengths[current_name] = len(''.join(current_seq))
    return lengths

def run_blast(query_fasta, db_fasta, db_name):
    """Run BLAST search against PlasmidFinder database"""
    output_file = os.path.join(OUTPUT_DIR, f"blast_{db_name}_temp.txt")
    
    cmd = [
        "blastn",
        "-query", query_fasta,
        "-subject", db_fasta,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-evalue", "1e-5",
        "-perc_identity", str(MIN_IDENTITY),
        "-out", output_file
    ]
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        print(f"  WARNING: BLAST timed out for {db_name}")
        return []
    
    results = []
    if os.path.exists(output_file):
        with open(output_file) as fh:
            for line in fh:
                fields = line.strip().split('\t')
                if len(fields) >= 14:
                    qseqid = fields[0]
                    sseqid = fields[1]
                    pident = float(fields[2])
                    length = int(fields[3])
                    qlen = int(fields[12])
                    slen = int(fields[13])
                    
                    # Calculate coverage relative to subject (database sequence)
                    coverage = (length / slen) * 100
                    
                    if pident >= MIN_IDENTITY and coverage >= MIN_COVERAGE:
                        results.append({
                            'contig': qseqid,
                            'replicon': sseqid,
                            'identity': pident,
                            'coverage': coverage,
                            'alignment_length': length,
                            'contig_length': qlen,
                            'subject_length': slen,
                            'database': db_name
                        })
        os.remove(output_file)
    
    return results

def classify_contig(contig_length):
    """Classify contig as chromosome or plasmid based on size"""
    if contig_length > 1_000_000:
        return "chromosome"
    else:
        return "plasmid"

def main():
    fasta_files = get_fasta_files()
    print(f"Found {len(fasta_files)} LR assemblies for replicon typing")
    print(f"Using PlasmidFinder databases: {len(DB_FILES)} databases")
    print(f"Thresholds: Identity >= {MIN_IDENTITY}%, Coverage >= {MIN_COVERAGE}%")
    print("=" * 80)
    
    all_results = []
    isolate_summaries = []
    
    for fasta_file in fasta_files:
        isolate_name = os.path.basename(fasta_file).replace('.fasta', '')
        print(f"\n{'='*60}")
        print(f"Analyzing: {isolate_name}")
        print(f"{'='*60}")
        
        # Parse contigs
        contigs = parse_fasta_contigs(fasta_file)
        print(f"  Total contigs: {len(contigs)}")
        
        # Classify contigs
        chromosomes = {k: v for k, v in contigs.items() if v > 1_000_000}
        plasmids = {k: v for k, v in contigs.items() if v <= 1_000_000}
        print(f"  Chromosomes: {len(chromosomes)}")
        print(f"  Plasmid candidates: {len(plasmids)}")
        
        # Run BLAST against all databases
        isolate_hits = []
        for db_name, db_fasta in DB_FILES.items():
            hits = run_blast(fasta_file, db_fasta, db_name)
            for hit in hits:
                hit['isolate'] = isolate_name
                hit['contig_type'] = classify_contig(hit['contig_length'])
                isolate_hits.append(hit)
        
        # Deduplicate - keep best hit per replicon per contig
        best_hits = {}
        for hit in isolate_hits:
            key = (hit['contig'], hit['replicon'])
            if key not in best_hits or hit['identity'] > best_hits[key]['identity']:
                best_hits[key] = hit
        
        # Further deduplicate - keep best hit per replicon type (across contigs)
        replicon_best = {}
        for hit in best_hits.values():
            # Extract replicon type name
            rep_name = hit['replicon'].split('_')[0] if '_' in hit['replicon'] else hit['replicon']
            key = (hit['isolate'], rep_name)
            if key not in replicon_best or hit['identity'] > replicon_best[key]['identity']:
                replicon_best[key] = hit
        
        unique_hits = list(best_hits.values())
        all_results.extend(unique_hits)
        
        # Print results
        if unique_hits:
            print(f"\n  REPLICON TYPES FOUND:")
            # Group by contig
            by_contig = defaultdict(list)
            for hit in unique_hits:
                by_contig[hit['contig']].append(hit)
            
            replicon_types = set()
            for contig, hits in sorted(by_contig.items(), key=lambda x: -contigs.get(x[0], 0)):
                contig_size = contigs.get(contig, 0)
                contig_type = "CHROMOSOME" if contig_size > 1_000_000 else f"PLASMID ({contig_size/1000:.1f} kb)"
                print(f"\n    {contig} [{contig_type}]:")
                for hit in sorted(hits, key=lambda x: -x['identity']):
                    print(f"      - {hit['replicon']} ({hit['database']}): {hit['identity']:.1f}% identity, {hit['coverage']:.1f}% coverage")
                    replicon_types.add(hit['replicon'])
            
            # Summary for this isolate
            isolate_summaries.append({
                'isolate': isolate_name,
                'total_contigs': len(contigs),
                'chromosomes': len(chromosomes),
                'plasmids': len(plasmids),
                'replicon_hits': len(unique_hits),
                'replicon_types': '; '.join(sorted(replicon_types)),
                'plasmid_sizes_kb': '; '.join([f"{v/1000:.1f}" for k, v in sorted(plasmids.items(), key=lambda x: -x[1])])
            })
        else:
            print(f"\n  NO REPLICON TYPES FOUND")
            isolate_summaries.append({
                'isolate': isolate_name,
                'total_contigs': len(contigs),
                'chromosomes': len(chromosomes),
                'plasmids': len(plasmids),
                'replicon_hits': 0,
                'replicon_types': 'None',
                'plasmid_sizes_kb': '; '.join([f"{v/1000:.1f}" for k, v in sorted(plasmids.items(), key=lambda x: -x[1])])
            })
    
    # Save detailed results
    detail_file = os.path.join(OUTPUT_DIR, "replicon_typing_detailed.csv")
    with open(detail_file, 'w', newline='') as fh:
        if all_results:
            writer = csv.DictWriter(fh, fieldnames=all_results[0].keys())
            writer.writeheader()
            writer.writerows(all_results)
    print(f"\n\nDetailed results saved to: {detail_file}")
    
    # Save summary
    summary_file = os.path.join(OUTPUT_DIR, "replicon_typing_summary.csv")
    with open(summary_file, 'w', newline='') as fh:
        if isolate_summaries:
            writer = csv.DictWriter(fh, fieldnames=isolate_summaries[0].keys())
            writer.writeheader()
            writer.writerows(isolate_summaries)
    print(f"Summary saved to: {summary_file}")
    
    # Print overall summary
    print(f"\n{'='*80}")
    print(f"OVERALL SUMMARY")
    print(f"{'='*80}")
    print(f"Total isolates analyzed: {len(fasta_files)}")
    print(f"Total replicon hits: {len(all_results)}")
    
    # Count replicon types across all isolates
    all_replicons = defaultdict(int)
    for hit in all_results:
        all_replicons[hit['replicon']] += 1
    
    print(f"\nReplicon type distribution:")
    for rep, count in sorted(all_replicons.items(), key=lambda x: -x[1]):
        print(f"  {rep}: {count}/{len(fasta_files)} isolates ({count/len(fasta_files)*100:.0f}%)")
    
    # Contig-level analysis: which plasmids carry which replicons
    print(f"\nPlasmid-Replicon associations:")
    plasmid_replicons = defaultdict(lambda: defaultdict(list))
    for hit in all_results:
        if hit['contig_type'] == 'plasmid':
            size_kb = hit['contig_length'] / 1000
            size_cat = f"{int(round(size_kb, -1))}kb"
            plasmid_replicons[size_cat][hit['replicon']].append(hit['isolate'])
    
    for size_cat in sorted(plasmid_replicons.keys(), key=lambda x: -int(x.replace('kb',''))):
        print(f"\n  ~{size_cat} plasmid:")
        for rep, isolates in sorted(plasmid_replicons[size_cat].items()):
            print(f"    {rep}: found in {len(isolates)} isolates")

if __name__ == '__main__':
    main()
