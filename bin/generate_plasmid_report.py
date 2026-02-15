#!/usr/bin/env python3
"""
Generate a comprehensive HTML report from PlasmidScope pipeline results.
Aggregates all analysis outputs into a single, publication-quality report.
"""

import argparse
import glob
import json
import os
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def parse_mob_suite_results(mob_dir):
    """Parse MOB-suite results from multiple samples."""
    results = []
    for report_file in glob.glob(os.path.join(mob_dir, '**/contig_report.txt'), recursive=True):
        try:
            with open(report_file) as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= len(header):
                        row = dict(zip(header, fields))
                        results.append(row)
        except Exception as e:
            print(f"Warning: Could not parse {report_file}: {e}", file=sys.stderr)
    return results


def parse_plasmidfinder_results(pf_dir):
    """Parse PlasmidFinder results."""
    results = []
    for result_file in glob.glob(os.path.join(pf_dir, '**/results_tab.tsv'), recursive=True):
        try:
            with open(result_file) as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= len(header):
                        row = dict(zip(header, fields))
                        results.append(row)
        except Exception as e:
            print(f"Warning: Could not parse {result_file}: {e}", file=sys.stderr)
    return results


def parse_abricate_results(abricate_dir):
    """Parse ABRicate results from multiple databases."""
    results = defaultdict(list)
    for tsv_file in glob.glob(os.path.join(abricate_dir, '**/*.tsv'), recursive=True):
        db_name = 'unknown'
        if 'card' in tsv_file.lower():
            db_name = 'card'
        elif 'vfdb' in tsv_file.lower():
            db_name = 'vfdb'
        elif 'plasmidfinder' in tsv_file.lower():
            db_name = 'plasmidfinder'

        try:
            with open(tsv_file) as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= len(header):
                        row = dict(zip(header, fields))
                        row['database'] = db_name
                        results[db_name].append(row)
        except Exception as e:
            print(f"Warning: Could not parse {tsv_file}: {e}", file=sys.stderr)
    return results


def parse_isfinder_results(is_dir):
    """Parse ISfinder BLAST results."""
    results = []
    for tsv_file in glob.glob(os.path.join(is_dir, '**/*_is_elements.tsv'), recursive=True):
        try:
            with open(tsv_file) as f:
                for line in f:
                    if line.startswith('No IS') or line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 13:
                        results.append({
                            'query': fields[0],
                            'subject': fields[1],
                            'identity': float(fields[2]),
                            'length': int(fields[3]),
                            'qstart': int(fields[6]),
                            'qend': int(fields[7]),
                            'evalue': float(fields[10]),
                            'bitscore': float(fields[11]),
                            'description': fields[12] if len(fields) > 12 else ''
                        })
        except Exception as e:
            print(f"Warning: Could not parse {tsv_file}: {e}", file=sys.stderr)
    return results


def parse_integron_results(integron_dir):
    """Parse IntegronFinder results."""
    results = []
    for tsv_file in glob.glob(os.path.join(integron_dir, '**/*.integrons'), recursive=True):
        try:
            with open(tsv_file) as f:
                for line in f:
                    if line.startswith('#') or line.startswith('ID'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 8:
                        results.append({
                            'id': fields[0],
                            'type': fields[3] if len(fields) > 3 else 'unknown',
                            'start': fields[4] if len(fields) > 4 else '',
                            'end': fields[5] if len(fields) > 5 else '',
                            'annotation': fields[7] if len(fields) > 7 else ''
                        })
        except Exception as e:
            print(f"Warning: Could not parse {tsv_file}: {e}", file=sys.stderr)
    return results


def generate_summary_tables(output_dir, mob_results, pf_results, abricate_results,
                            is_results, integron_results):
    """Generate summary TSV tables."""
    tables_dir = os.path.join(output_dir, 'summary_tables')
    os.makedirs(tables_dir, exist_ok=True)

    # Plasmid inventory table
    with open(os.path.join(tables_dir, 'plasmid_inventory.tsv'), 'w') as f:
        f.write("sample\tcontig\ttype\tsize_bp\treplicon\tmob_type\tmpf_type\tmobility\n")
        for r in mob_results:
            f.write(f"{r.get('sample_id', 'NA')}\t{r.get('contig_id', 'NA')}\t"
                    f"{r.get('molecule_type', 'NA')}\t{r.get('size', 'NA')}\t"
                    f"{r.get('rep_type(s)', 'NA')}\t{r.get('relaxase_type(s)', 'NA')}\t"
                    f"{r.get('mpf_type', 'NA')}\t{r.get('predicted_mobility', 'NA')}\n")

    # Replicon type summary
    with open(os.path.join(tables_dir, 'replicon_types.tsv'), 'w') as f:
        f.write("replicon_type\tcount\tsamples\n")
        rep_counts = defaultdict(list)
        for r in pf_results:
            rep = r.get('Plasmid', r.get('plasmid', 'unknown'))
            sample = r.get('sample', 'unknown')
            rep_counts[rep].append(sample)
        for rep, samples in sorted(rep_counts.items(), key=lambda x: -len(x[1])):
            f.write(f"{rep}\t{len(samples)}\t{';'.join(set(samples))}\n")

    # AMR gene summary
    with open(os.path.join(tables_dir, 'amr_genes.tsv'), 'w') as f:
        f.write("gene\tproduct\tcontig\tidentity\tcoverage\tdatabase\n")
        for db, hits in abricate_results.items():
            for h in hits:
                f.write(f"{h.get('GENE', 'NA')}\t{h.get('PRODUCT', 'NA')}\t"
                        f"{h.get('SEQUENCE', 'NA')}\t{h.get('%IDENTITY', 'NA')}\t"
                        f"{h.get('%COVERAGE', 'NA')}\t{db}\n")

    # IS element summary
    with open(os.path.join(tables_dir, 'is_elements.tsv'), 'w') as f:
        f.write("is_element\tcontig\tidentity\tlength\tqstart\tqend\tevalue\n")
        for r in is_results:
            f.write(f"{r['subject']}\t{r['query']}\t{r['identity']:.1f}\t"
                    f"{r['length']}\t{r['qstart']}\t{r['qend']}\t{r['evalue']}\n")

    # Integron summary
    with open(os.path.join(tables_dir, 'integrons.tsv'), 'w') as f:
        f.write("id\ttype\tstart\tend\tannotation\n")
        for r in integron_results:
            f.write(f"{r['id']}\t{r['type']}\t{r['start']}\t{r['end']}\t{r['annotation']}\n")

    return tables_dir


def generate_figures(output_dir, mob_results, pf_results, abricate_results, is_results):
    """Generate publication-quality figures."""
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    if not HAS_MATPLOTLIB:
        print("Warning: matplotlib not available, skipping figure generation", file=sys.stderr)
        # Create placeholder
        with open(os.path.join(figures_dir, 'README.txt'), 'w') as f:
            f.write("Install matplotlib to generate figures: pip install matplotlib\n")
        return figures_dir

    # Figure 1: Replicon type distribution
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('PlasmidScope - Plasmid Analysis Summary', fontsize=16, fontweight='bold')

    # Panel A: Replicon types
    rep_counts = defaultdict(int)
    for r in pf_results:
        rep = r.get('Plasmid', r.get('plasmid', 'unknown'))
        rep_counts[rep] += 1
    if rep_counts:
        sorted_reps = sorted(rep_counts.items(), key=lambda x: -x[1])[:15]
        labels, values = zip(*sorted_reps)
        axes[0, 0].barh(range(len(labels)), values, color='#2196F3')
        axes[0, 0].set_yticks(range(len(labels)))
        axes[0, 0].set_yticklabels(labels, fontsize=8)
        axes[0, 0].set_xlabel('Count')
        axes[0, 0].set_title('A. Plasmid Replicon Types', fontweight='bold')
        axes[0, 0].invert_yaxis()

    # Panel B: Mobility prediction
    mob_counts = defaultdict(int)
    for r in mob_results:
        mob = r.get('predicted_mobility', 'unknown')
        mob_counts[mob] += 1
    if mob_counts:
        labels = list(mob_counts.keys())
        values = list(mob_counts.values())
        colors = ['#4CAF50', '#FF9800', '#F44336', '#9E9E9E']
        axes[0, 1].pie(values, labels=labels, autopct='%1.1f%%',
                       colors=colors[:len(labels)], startangle=90)
        axes[0, 1].set_title('B. Plasmid Mobility', fontweight='bold')

    # Panel C: Top AMR genes
    gene_counts = defaultdict(int)
    for db, hits in abricate_results.items():
        if db == 'card':
            for h in hits:
                gene = h.get('GENE', 'unknown')
                gene_counts[gene] += 1
    if gene_counts:
        sorted_genes = sorted(gene_counts.items(), key=lambda x: -x[1])[:15]
        labels, values = zip(*sorted_genes)
        axes[1, 0].barh(range(len(labels)), values, color='#E91E63')
        axes[1, 0].set_yticks(range(len(labels)))
        axes[1, 0].set_yticklabels(labels, fontsize=8)
        axes[1, 0].set_xlabel('Count')
        axes[1, 0].set_title('C. AMR Genes (CARD)', fontweight='bold')
        axes[1, 0].invert_yaxis()

    # Panel D: IS element families
    is_counts = defaultdict(int)
    for r in is_results:
        is_name = r['subject'].split('_')[0] if '_' in r['subject'] else r['subject']
        is_counts[is_name] += 1
    if is_counts:
        sorted_is = sorted(is_counts.items(), key=lambda x: -x[1])[:15]
        labels, values = zip(*sorted_is)
        axes[1, 1].barh(range(len(labels)), values, color='#FF5722')
        axes[1, 1].set_yticks(range(len(labels)))
        axes[1, 1].set_yticklabels(labels, fontsize=8)
        axes[1, 1].set_xlabel('Count')
        axes[1, 1].set_title('D. Insertion Sequences', fontweight='bold')
        axes[1, 1].invert_yaxis()

    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, 'plasmid_summary.png'), dpi=300, bbox_inches='tight')
    plt.close()

    return figures_dir


def generate_html_report(output_dir, mob_results, pf_results, abricate_results,
                         is_results, integron_results, tables_dir, figures_dir):
    """Generate the main HTML report."""

    n_samples = len(set(r.get('sample_id', '') for r in mob_results))
    n_plasmids = sum(1 for r in mob_results if r.get('molecule_type') == 'plasmid')
    n_amr = sum(len(hits) for db, hits in abricate_results.items() if db == 'card')
    n_is = len(is_results)
    n_integrons = len(integron_results)

    # Replicon summary
    rep_counts = defaultdict(int)
    for r in pf_results:
        rep = r.get('Plasmid', r.get('plasmid', 'unknown'))
        rep_counts[rep] += 1

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PlasmidScope Analysis Report</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background: #f8f9fa; color: #333; line-height: 1.6; }}
        .header {{ background: linear-gradient(135deg, #1a237e, #0d47a1); color: white; padding: 40px; text-align: center; }}
        .header h1 {{ font-size: 2.5em; margin-bottom: 10px; }}
        .header p {{ font-size: 1.1em; opacity: 0.9; }}
        .container {{ max-width: 1200px; margin: 0 auto; padding: 20px; }}
        .card {{ background: white; border-radius: 12px; box-shadow: 0 2px 10px rgba(0,0,0,0.08); padding: 30px; margin: 20px 0; }}
        .card h2 {{ color: #1a237e; border-bottom: 2px solid #e3f2fd; padding-bottom: 10px; margin-bottom: 20px; }}
        .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 20px 0; }}
        .stat-box {{ background: linear-gradient(135deg, #e3f2fd, #bbdefb); border-radius: 10px; padding: 20px; text-align: center; }}
        .stat-box .number {{ font-size: 2.5em; font-weight: bold; color: #1565c0; }}
        .stat-box .label {{ font-size: 0.9em; color: #555; margin-top: 5px; }}
        table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
        th {{ background: #1a237e; color: white; padding: 12px; text-align: left; font-size: 0.9em; }}
        td {{ padding: 10px 12px; border-bottom: 1px solid #e0e0e0; font-size: 0.9em; }}
        tr:hover {{ background: #f5f5f5; }}
        .badge {{ display: inline-block; padding: 3px 10px; border-radius: 12px; font-size: 0.8em; font-weight: bold; }}
        .badge-conj {{ background: #c8e6c9; color: #2e7d32; }}
        .badge-mob {{ background: #fff9c4; color: #f57f17; }}
        .badge-nonmob {{ background: #ffcdd2; color: #c62828; }}
        .footer {{ text-align: center; padding: 30px; color: #888; font-size: 0.9em; }}
        img {{ max-width: 100%; height: auto; border-radius: 8px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>PlasmidScope</h1>
        <p>Comprehensive Plasmid Genomic Analysis Report</p>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>

    <div class="container">
        <div class="stats">
            <div class="stat-box">
                <div class="number">{n_samples}</div>
                <div class="label">Samples Analyzed</div>
            </div>
            <div class="stat-box">
                <div class="number">{n_plasmids}</div>
                <div class="label">Plasmids Identified</div>
            </div>
            <div class="stat-box">
                <div class="number">{len(rep_counts)}</div>
                <div class="label">Unique Replicon Types</div>
            </div>
            <div class="stat-box">
                <div class="number">{n_amr}</div>
                <div class="label">AMR Gene Hits</div>
            </div>
            <div class="stat-box">
                <div class="number">{n_is}</div>
                <div class="label">IS Elements</div>
            </div>
            <div class="stat-box">
                <div class="number">{n_integrons}</div>
                <div class="label">Integrons</div>
            </div>
        </div>

        <div class="card">
            <h2>1. Plasmid Inventory</h2>
            <table>
                <tr><th>Sample</th><th>Contig</th><th>Size (bp)</th><th>Replicon</th><th>MOB Type</th><th>Mobility</th></tr>"""

    for r in mob_results[:50]:
        mobility = r.get('predicted_mobility', 'unknown')
        badge_class = 'badge-conj' if 'conj' in mobility.lower() else \
                      'badge-mob' if 'mob' in mobility.lower() else 'badge-nonmob'
        html += f"""
                <tr>
                    <td>{r.get('sample_id', 'NA')}</td>
                    <td>{r.get('contig_id', 'NA')}</td>
                    <td>{r.get('size', 'NA')}</td>
                    <td>{r.get('rep_type(s)', '-')}</td>
                    <td>{r.get('relaxase_type(s)', '-')}</td>
                    <td><span class="badge {badge_class}">{mobility}</span></td>
                </tr>"""

    html += f"""
            </table>
        </div>

        <div class="card">
            <h2>2. Replicon Type Distribution</h2>
            <table>
                <tr><th>Replicon Type</th><th>Count</th><th>Prevalence (%)</th></tr>"""

    total_reps = sum(rep_counts.values())
    for rep, count in sorted(rep_counts.items(), key=lambda x: -x[1]):
        pct = (count / total_reps * 100) if total_reps > 0 else 0
        html += f"""
                <tr><td>{rep}</td><td>{count}</td><td>{pct:.1f}%</td></tr>"""

    html += f"""
            </table>
        </div>

        <div class="card">
            <h2>3. Summary Figures</h2>
            <img src="figures/plasmid_summary.png" alt="Plasmid Analysis Summary">
        </div>

        <div class="card">
            <h2>4. Pipeline Information</h2>
            <p><strong>Pipeline:</strong> PlasmidScope v1.0.0</p>
            <p><strong>Tools used:</strong> Bakta, ABRicate, AMRFinderPlus, PlasmidFinder, MOB-suite, ISfinder, IntegronFinder</p>
            <p><strong>Databases:</strong> CARD, VFDB, PlasmidFinder, ISfinder, mobileOG-db</p>
            <p><strong>Report generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>

    <div class="footer">
        <p>PlasmidScope - Comprehensive Plasmid Genomic Analysis Pipeline</p>
        <p>Report generated automatically. For questions, see documentation.</p>
    </div>
</body>
</html>"""

    report_path = os.path.join(output_dir, 'plasmid_analysis_report.html')
    with open(report_path, 'w') as f:
        f.write(html)

    return report_path


def main():
    parser = argparse.ArgumentParser(description='Generate PlasmidScope analysis report')
    parser.add_argument('--bakta', nargs='*', default=[], help='Bakta summary files')
    parser.add_argument('--abricate', nargs='*', default=[], help='ABRicate result directories')
    parser.add_argument('--amrfinder', nargs='*', default=[], help='AMRFinderPlus result files')
    parser.add_argument('--plasmidfinder', nargs='*', default=[], help='PlasmidFinder result dirs')
    parser.add_argument('--mob_suite', nargs='*', default=[], help='MOB-suite result directories')
    parser.add_argument('--isfinder', nargs='*', default=[], help='ISfinder result files')
    parser.add_argument('--integrons', nargs='*', default=[], help='IntegronFinder result dirs')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Parse all results
    mob_results = []
    for d in args.mob_suite:
        mob_results.extend(parse_mob_suite_results(d))

    pf_results = []
    for d in args.plasmidfinder:
        pf_results.extend(parse_plasmidfinder_results(d))

    abricate_results = defaultdict(list)
    for d in args.abricate:
        for db, hits in parse_abricate_results(d).items():
            abricate_results[db].extend(hits)

    is_results = []
    for d in args.isfinder:
        is_results.extend(parse_isfinder_results(d))

    integron_results = []
    for d in args.integrons:
        integron_results.extend(parse_integron_results(d))

    # Generate outputs
    tables_dir = generate_summary_tables(args.output_dir, mob_results, pf_results,
                                          abricate_results, is_results, integron_results)
    figures_dir = generate_figures(args.output_dir, mob_results, pf_results,
                                   abricate_results, is_results)
    report_path = generate_html_report(args.output_dir, mob_results, pf_results,
                                        abricate_results, is_results, integron_results,
                                        tables_dir, figures_dir)

    print(f"Report generated: {report_path}")
    print(f"Summary tables: {tables_dir}")
    print(f"Figures: {figures_dir}")


if __name__ == '__main__':
    main()
