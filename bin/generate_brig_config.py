#!/usr/bin/env python3
"""
Generate BRIG configuration XML and run BRIG for circular plasmid comparison.
Alternative: Generate pyCircos/pyGenomeViz circular plots if BRIG is not available.
"""

import argparse
import glob
import os
import sys
import subprocess


def generate_brig_xml(reference, query_dir, output_dir):
    """Generate BRIG configuration XML file."""
    queries = sorted(glob.glob(os.path.join(query_dir, '*.fasta')) +
                     glob.glob(os.path.join(query_dir, '*.fa')))

    colors = [
        '0,0,255', '255,0,0', '0,128,0', '255,165,0', '128,0,128',
        '0,128,128', '255,0,255', '128,128,0', '0,0,128', '128,0,0',
        '0,255,255', '255,255,0', '64,64,64', '192,0,0', '0,192,0'
    ]

    rings = []
    for i, query in enumerate(queries):
        color = colors[i % len(colors)]
        name = os.path.basename(query).replace('.fasta', '').replace('.fa', '')
        rings.append(f'''    <ring>
        <name>{name}</name>
        <colour>{color}</colour>
        <sequence>{os.path.abspath(query)}</sequence>
        <upperIdentity>100</upperIdentity>
        <lowerIdentity>70</lowerIdentity>
    </ring>''')

    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<brig>
    <reference>{os.path.abspath(reference)}</reference>
    <output>{os.path.join(os.path.abspath(output_dir), 'brig_comparison.svg')}</output>
    <blastOptions>-evalue 1e-5 -num_threads 4</blastOptions>
    <imageWidth>3000</imageWidth>
    <imageHeight>3000</imageHeight>
    <title>PlasmidScope - Plasmid Comparison</title>
{''.join(rings)}
</brig>"""

    config_path = os.path.join(output_dir, 'brig_config.xml')
    with open(config_path, 'w') as f:
        f.write(xml)

    return config_path


def generate_pygenomeviz_circular(reference, query_dir, output_dir):
    """Generate circular comparison using pyGenomeViz as BRIG alternative."""
    try:
        from pygenomeviz import GenomeViz
        from Bio import SeqIO
    except ImportError:
        print("pyGenomeViz or BioPython not available. Generating BLAST-based comparison instead.",
              file=sys.stderr)
        return generate_blast_comparison(reference, query_dir, output_dir)

    queries = sorted(glob.glob(os.path.join(query_dir, '*.fasta')) +
                     glob.glob(os.path.join(query_dir, '*.fa')))

    # Create BLAST databases and run comparisons
    for query in queries:
        name = os.path.basename(query).replace('.fasta', '').replace('.fa', '')
        blast_out = os.path.join(output_dir, f'{name}_blast.tsv')
        subprocess.run([
            'blastn', '-query', reference, '-subject', query,
            '-outfmt', '6', '-evalue', '1e-5', '-out', blast_out
        ], check=False)

    print(f"BLAST comparisons generated in {output_dir}")
    return output_dir


def generate_blast_comparison(reference, query_dir, output_dir):
    """Generate BLAST-based comparison TSV files for manual BRIG import."""
    queries = sorted(glob.glob(os.path.join(query_dir, '*.fasta')) +
                     glob.glob(os.path.join(query_dir, '*.fa')))

    os.makedirs(output_dir, exist_ok=True)

    for query in queries:
        name = os.path.basename(query).replace('.fasta', '').replace('.fa', '')
        blast_out = os.path.join(output_dir, f'{name}_vs_reference.tsv')

        cmd = [
            'blastn', '-query', reference, '-subject', query,
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            '-evalue', '1e-5',
            '-out', blast_out
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"  Generated: {blast_out}")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"  Warning: BLAST failed for {name}: {e}", file=sys.stderr)

    # Generate summary
    summary_path = os.path.join(output_dir, 'comparison_summary.txt')
    with open(summary_path, 'w') as f:
        f.write(f"Reference: {reference}\n")
        f.write(f"Queries: {len(queries)}\n\n")
        for query in queries:
            name = os.path.basename(query).replace('.fasta', '').replace('.fa', '')
            f.write(f"  - {name}\n")

    return output_dir


def main():
    parser = argparse.ArgumentParser(description='Generate BRIG configuration for plasmid comparison')
    parser.add_argument('--reference', required=True, help='Reference plasmid FASTA')
    parser.add_argument('--query_dir', required=True, help='Directory with query plasmid FASTAs')
    parser.add_argument('--output', required=True, help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Generate BRIG config
    config = generate_brig_xml(args.reference, args.query_dir, args.output)
    print(f"BRIG config generated: {config}")

    # Also generate BLAST comparisons
    generate_blast_comparison(args.reference, args.query_dir, args.output)
    print(f"BLAST comparisons generated in {args.output}")


if __name__ == '__main__':
    main()
