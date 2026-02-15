#!/usr/bin/env python3
"""
Create BRIG-style circular maps and EasyFig-style linear comparison figures
for the genomic environment of key resistance genes (blaKPC-2, blaNDM-1, blaCTX-M-65)
across 8 ST11 K. pneumoniae isolates.

Uses pyGenomeViz for EasyFig-style linear comparisons
Uses pyCirclize for BRIG-style circular maps
"""

import os
import sys
import json
import csv
import io
import subprocess
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, SimpleLocation

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

BASE = Path("/home/ubuntu/plasmid_analysis")
LR_DIR = BASE / "LR_assemblies"
BAKTA_DIR = BASE / "04_AMR_Plasmid_Association" / "Bakta_Annotations"
OUT_DIR = BASE / "08_BRIG_EasyFig"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Contig mapping: Bakta contig names -> original FASTA contig numbers
# Bakta renames contigs as contig_1, contig_2, etc. matching order in FASTA
# So contig_1 = first sequence, contig_2 = second, etc.

SAMPLES = ['JKPB244', 'JKPB284', 'JKPB381', 'JKPR282', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']

FASTA_MAP = {
    'JKPB244': '1_JKPB244_LR.fasta',
    'JKPB284': '2_JKPB284_LR.fasta',
    'JKPB381': '4_JKPB381_LR.fasta',
    'JKPR282': '5_JKPR282_LR.fasta',
    'KP21802': '6_KP21802_LR.fasta',
    'KP21915': '7_KP21915_LR.fasta',
    'MDRKP121': '8_MDRKP121_LR.fasta',
    'MDRKP122': '9_MDRKP122_LR.fasta',
}

# KPC plasmid: contig_3 (~130kb) for most, contig_2 for JKPR282
KPC_CONTIG = {
    'JKPB244': 3, 'JKPB284': 3, 'JKPB381': 3, 'JKPR282': 2,
    'KP21802': 3, 'KP21915': 3, 'MDRKP121': 3, 'MDRKP122': 3
}

# NDM plasmid: contig_2 (~350kb) for those that have it
NDM_CONTIG = {
    'JKPB244': 2, 'JKPB284': 2, 'KP21802': 2,
    'KP21915': 2, 'MDRKP121': 2, 'MDRKP122': 2
}


def extract_contig_sequence(sample, contig_num):
    """Extract a specific contig sequence from the LR assembly FASTA."""
    fasta_file = LR_DIR / FASTA_MAP[sample]
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if contig_num <= len(records):
        return records[contig_num - 1]
    return None


def get_genes_on_contig(sample, contig_name):
    """Get all annotated genes on a specific contig from TSV."""
    tsv_file = BAKTA_DIR / sample / f"{sample}.tsv"
    genes = []
    
    if tsv_file.exists():
        with open(tsv_file) as f:
            lines = [l for l in f if not l.startswith('# ')]
        reader = csv.DictReader(io.StringIO(''.join(lines)), delimiter='\t')
        for row in reader:
            contig = row.get('#Sequence Id', row.get('Sequence Id', ''))
            if contig == contig_name:
                genes.append({
                    'gene': row.get('Gene', '') or '',
                    'product': row.get('Product', '') or '',
                    'start': int(row.get('Start', 0)),
                    'stop': int(row.get('Stop', 0)),
                    'strand': 1 if row.get('Strand', '+') == '+' else -1,
                    'type': row.get('Type', '')
                })
    
    return genes


def classify_gene_color(gene, product):
    """Classify gene for coloring in the figure."""
    gene_lower = (gene or '').lower()
    product_lower = (product or '').lower()
    
    # Carbapenemases - RED
    if any(x in gene_lower for x in ['blakpc', 'blandm', 'blavim', 'blaimp', 'blaoxa-48']):
        return '#E63946', 'Carbapenemase'
    # Beta-lactamases - ORANGE
    elif gene_lower.startswith('bla'):
        return '#F4A261', 'Beta-lactamase'
    # Aminoglycoside resistance - BLUE
    elif any(x in gene_lower for x in ['aac', 'aph', 'ant', 'rmt', 'arm']):
        return '#457B9D', 'Aminoglycoside resistance'
    # Quinolone resistance - TEAL
    elif any(x in gene_lower for x in ['qnr', 'oqx']):
        return '#2A9D8F', 'Quinolone resistance'
    # Other AMR
    elif any(x in gene_lower for x in ['sul', 'dfr', 'tet', 'mph', 'msr', 'ere', 'cat', 'flo', 'fos', 'ble', 'arr']):
        return '#E9C46A', 'Other AMR'
    # Transposases / IS elements - PURPLE
    elif 'transpos' in product_lower or gene_lower.startswith('is') or 'insertion' in product_lower:
        return '#9B59B6', 'Transposase/IS'
    # Integrases - DARK RED
    elif 'integrase' in product_lower or 'recombinase' in product_lower:
        return '#C0392B', 'Integrase/Recombinase'
    # Conjugation/Transfer - DARK GREEN
    elif any(x in product_lower for x in ['conjugat', 'transfer', 'relaxase', 'type iv']):
        return '#27AE60', 'Conjugation/Transfer'
    # Replication
    elif any(x in product_lower for x in ['replica', 'partitio']):
        return '#3498DB', 'Replication/Partition'
    # Hypothetical
    elif 'hypothetical' in product_lower:
        return '#D5D8DC', 'Hypothetical'
    else:
        return '#AEB6BF', 'Other'


def extract_region(sample, contig_num, center_pos, flank=20000):
    """Extract a genomic region around a gene of interest."""
    record = extract_contig_sequence(sample, contig_num)
    if record is None:
        return None, []
    
    contig_name = f"contig_{contig_num}"
    genes = get_genes_on_contig(sample, contig_name)
    
    start = max(0, center_pos - flank)
    end = min(len(record.seq), center_pos + flank)
    
    # Filter genes in region
    region_genes = []
    for g in genes:
        if g['start'] >= start and g['stop'] <= end:
            region_genes.append({
                **g,
                'start': g['start'] - start,
                'stop': g['stop'] - start,
            })
    
    # Extract subsequence
    sub_seq = record.seq[start:end]
    sub_record = SeqRecord(sub_seq, id=f"{sample}_{contig_name}_{start}-{end}", 
                           description=f"{sample} {contig_name} region {start}-{end}")
    
    return sub_record, region_genes


# ============================================================
# EasyFig-style Linear Comparison (pyGenomeViz)
# ============================================================
def create_easyfig_kpc_environment():
    """Create EasyFig-style linear comparison of blaKPC-2 genomic environment."""
    from pygenomeviz import GenomeViz
    
    print("Creating EasyFig-style KPC-2 environment comparison...")
    
    # KPC-2 positions (center of gene)
    kpc_centers = {
        'JKPB244': (3, 42786), 'JKPB284': (3, 42786), 'JKPB381': (3, 40691),
        'JKPR282': (2, 42786), 'KP21802': (3, 42786), 'KP21915': (3, 42786),
        'MDRKP121': (3, 42785), 'MDRKP122': (3, 42785)
    }
    
    flank = 15000  # 15kb flanking region
    
    gv = GenomeViz(
        fig_width=16,
        fig_track_height=0.7,
        feature_track_ratio=0.3,
        track_align_type="center",
    )
    
    sample_data = {}
    for sample in SAMPLES:
        contig_num, center = kpc_centers[sample]
        record, genes = extract_region(sample, contig_num, center, flank)
        if record:
            sample_data[sample] = (record, genes)
    
    # Add tracks for each sample
    for sample in SAMPLES:
        if sample not in sample_data:
            continue
        record, genes = sample_data[sample]
        
        track = gv.add_feature_track(sample, len(record.seq))
        
        for g in genes:
            color, category = classify_gene_color(g['gene'], g['product'])
            label = g['gene'] if g['gene'] else ''
            
            # Only label significant genes
            show_label = False
            if any(x in (g['gene'] or '').lower() for x in ['bla', 'aac', 'aph', 'rmt', 'arm', 'qnr', 'sul', 'dfr', 'tet', 'mph', 'fos', 'cat', 'ble', 'arr']):
                show_label = True
            if 'transpos' in (g['product'] or '').lower() and g['gene']:
                show_label = True
            
            track.add_feature(
                g['start'], g['stop'], g['strand'],
                label=label if show_label else "",
                text_kws=dict(size=6, rotation=45) if show_label else None,
                facecolor=color,
                linewidth=0.5,
                arrow_shaft_ratio=0.5,
            )
    
    # Add BLAST-based similarity links between adjacent tracks
    track_names = [s for s in SAMPLES if s in sample_data]
    for i in range(len(track_names) - 1):
        s1, s2 = track_names[i], track_names[i+1]
        rec1, _ = sample_data[s1]
        rec2, _ = sample_data[s2]
        
        # Save temp FASTA files for BLAST
        f1 = OUT_DIR / f"temp_{s1}_kpc.fasta"
        f2 = OUT_DIR / f"temp_{s2}_kpc.fasta"
        SeqIO.write(rec1, f1, "fasta")
        SeqIO.write(rec2, f2, "fasta")
        
        # Run BLAST
        blast_out = OUT_DIR / f"blast_{s1}_{s2}_kpc.txt"
        cmd = f"blastn -query {f1} -subject {f2} -outfmt 6 -evalue 1e-5 -max_target_seqs 100 -out {blast_out}"
        subprocess.run(cmd, shell=True, capture_output=True)
        
        # Parse BLAST results and add links
        if blast_out.exists():
            with open(blast_out) as bf:
                for line in bf:
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:
                        identity = float(parts[2])
                        q_start, q_end = int(parts[6]), int(parts[7])
                        s_start, s_end = int(parts[8]), int(parts[9])
                        aln_len = int(parts[3])
                        
                        if aln_len >= 500 and identity >= 80:
                            # Color by identity
                            if identity >= 99:
                                color = "#3498DB"
                                alpha = 0.4
                            elif identity >= 95:
                                color = "#85C1E9"
                                alpha = 0.3
                            elif identity >= 90:
                                color = "#AED6F1"
                                alpha = 0.2
                            else:
                                color = "#D6EAF8"
                                alpha = 0.15
                            
                            gv.add_link(
                                (s1, q_start, q_end),
                                (s2, s_start, s_end),
                                color=color,
                                v=alpha,
                            )
        
        # Clean up temp files
        f1.unlink(missing_ok=True)
        f2.unlink(missing_ok=True)
        blast_out.unlink(missing_ok=True)
    
    # Add legend
    fig = gv.plotfig()
    
    legend_items = [
        ('#E63946', 'Carbapenemase'),
        ('#F4A261', 'Beta-lactamase'),
        ('#457B9D', 'Aminoglycoside resistance'),
        ('#2A9D8F', 'Quinolone resistance'),
        ('#E9C46A', 'Other AMR'),
        ('#9B59B6', 'Transposase/IS'),
        ('#C0392B', 'Integrase/Recombinase'),
        ('#27AE60', 'Conjugation/Transfer'),
        ('#3498DB', 'Replication/Partition'),
        ('#AEB6BF', 'Other'),
    ]
    
    handles = [mpatches.Patch(color=c, label=l) for c, l in legend_items]
    fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=7,
              framealpha=0.9, bbox_to_anchor=(0.5, -0.02))
    
    fig.suptitle('Genomic Environment of blaKPC-2 (~130kb IncFII+IncR Plasmid)\n'
                 'Linear Comparison Across 8 ST11 K. pneumoniae Isolates',
                 fontsize=12, fontweight='bold', y=1.02)
    
    fig.savefig(OUT_DIR / 'EasyFig_KPC2_Environment.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / 'EasyFig_KPC2_Environment.pdf', bbox_inches='tight')
    plt.close()
    print("  EasyFig KPC-2 environment - DONE")


def create_easyfig_ndm_environment():
    """Create EasyFig-style linear comparison of blaNDM-1 genomic environment."""
    from pygenomeviz import GenomeViz
    
    print("Creating EasyFig-style NDM-1 environment comparison...")
    
    # NDM-1 positions (center of gene)
    ndm_centers = {
        'JKPB244': (2, 219318), 'JKPB284': (2, 218422),
        'KP21802': (2, 211939), 'KP21915': (2, 211191),
        'MDRKP121': (2, 217231), 'MDRKP122': (2, 218383)
    }
    
    ndm_samples = [s for s in SAMPLES if s in ndm_centers]
    flank = 15000
    
    gv = GenomeViz(
        fig_width=16,
        fig_track_height=0.7,
        feature_track_ratio=0.3,
        track_align_type="center",
    )
    
    sample_data = {}
    for sample in ndm_samples:
        contig_num, center = ndm_centers[sample]
        record, genes = extract_region(sample, contig_num, center, flank)
        if record:
            sample_data[sample] = (record, genes)
    
    for sample in ndm_samples:
        if sample not in sample_data:
            continue
        record, genes = sample_data[sample]
        
        track = gv.add_feature_track(sample, len(record.seq))
        
        for g in genes:
            color, category = classify_gene_color(g['gene'], g['product'])
            label = g['gene'] if g['gene'] else ''
            
            show_label = False
            if any(x in (g['gene'] or '').lower() for x in ['bla', 'aac', 'aph', 'rmt', 'arm', 'qnr', 'sul', 'dfr', 'tet', 'mph', 'msr', 'fos', 'cat', 'ble', 'arr']):
                show_label = True
            if 'transpos' in (g['product'] or '').lower() and g['gene']:
                show_label = True
            
            track.add_feature(
                g['start'], g['stop'], g['strand'],
                label=label if show_label else "",
                text_kws=dict(size=6, rotation=45) if show_label else None,
                facecolor=color,
                linewidth=0.5,
                arrow_shaft_ratio=0.5,
            )
    
    # Add BLAST links
    track_names = [s for s in ndm_samples if s in sample_data]
    for i in range(len(track_names) - 1):
        s1, s2 = track_names[i], track_names[i+1]
        rec1, _ = sample_data[s1]
        rec2, _ = sample_data[s2]
        
        f1 = OUT_DIR / f"temp_{s1}_ndm.fasta"
        f2 = OUT_DIR / f"temp_{s2}_ndm.fasta"
        SeqIO.write(rec1, f1, "fasta")
        SeqIO.write(rec2, f2, "fasta")
        
        blast_out = OUT_DIR / f"blast_{s1}_{s2}_ndm.txt"
        cmd = f"blastn -query {f1} -subject {f2} -outfmt 6 -evalue 1e-5 -max_target_seqs 100 -out {blast_out}"
        subprocess.run(cmd, shell=True, capture_output=True)
        
        if blast_out.exists():
            with open(blast_out) as bf:
                for line in bf:
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:
                        identity = float(parts[2])
                        q_start, q_end = int(parts[6]), int(parts[7])
                        s_start, s_end = int(parts[8]), int(parts[9])
                        aln_len = int(parts[3])
                        
                        if aln_len >= 500 and identity >= 80:
                            if identity >= 99:
                                color = "#3498DB"
                                alpha = 0.4
                            elif identity >= 95:
                                color = "#85C1E9"
                                alpha = 0.3
                            elif identity >= 90:
                                color = "#AED6F1"
                                alpha = 0.2
                            else:
                                color = "#D6EAF8"
                                alpha = 0.15
                            
                            gv.add_link(
                                (s1, q_start, q_end),
                                (s2, s_start, s_end),
                                color=color,
                                v=alpha,
                            )
        
        f1.unlink(missing_ok=True)
        f2.unlink(missing_ok=True)
        blast_out.unlink(missing_ok=True)
    
    fig = gv.plotfig()
    
    legend_items = [
        ('#E63946', 'Carbapenemase'),
        ('#F4A261', 'Beta-lactamase'),
        ('#457B9D', 'Aminoglycoside resistance'),
        ('#2A9D8F', 'Quinolone resistance'),
        ('#E9C46A', 'Other AMR'),
        ('#9B59B6', 'Transposase/IS'),
        ('#C0392B', 'Integrase/Recombinase'),
        ('#27AE60', 'Conjugation/Transfer'),
        ('#3498DB', 'Replication/Partition'),
        ('#AEB6BF', 'Other'),
    ]
    
    handles = [mpatches.Patch(color=c, label=l) for c, l in legend_items]
    fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=7,
              framealpha=0.9, bbox_to_anchor=(0.5, -0.02))
    
    fig.suptitle('Genomic Environment of blaNDM-1 (~350kb IncFIB(pNDM-Mar)+IncHI1B Plasmid)\n'
                 'Linear Comparison Across 6 ST11 K. pneumoniae Isolates',
                 fontsize=12, fontweight='bold', y=1.02)
    
    fig.savefig(OUT_DIR / 'EasyFig_NDM1_Environment.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / 'EasyFig_NDM1_Environment.pdf', bbox_inches='tight')
    plt.close()
    print("  EasyFig NDM-1 environment - DONE")


# ============================================================
# BRIG-style Circular Map (pyCirclize)
# ============================================================
def create_brig_kpc_plasmid():
    """Create BRIG-style circular map of the KPC plasmid (~130kb)."""
    from pycirclize import Circos
    from pycirclize.utils import ColorCycler
    
    print("Creating BRIG-style circular map for KPC plasmid...")
    
    # Use KP21915 as reference (has full GBFF annotation)
    ref_sample = 'KP21915'
    ref_contig_num = KPC_CONTIG[ref_sample]
    ref_record = extract_contig_sequence(ref_sample, ref_contig_num)
    ref_genes = get_genes_on_contig(ref_sample, f"contig_{ref_contig_num}")
    ref_len = len(ref_record.seq)
    
    # Save reference FASTA
    ref_fasta = OUT_DIR / "ref_kpc_plasmid.fasta"
    SeqIO.write(ref_record, ref_fasta, "fasta")
    
    # Create BLAST databases and run comparisons
    subprocess.run(f"makeblastdb -in {ref_fasta} -dbtype nucl -out {OUT_DIR}/ref_kpc_db", 
                   shell=True, capture_output=True)
    
    blast_results = {}
    for sample in SAMPLES:
        if sample == ref_sample:
            continue
        contig_num = KPC_CONTIG[sample]
        query_record = extract_contig_sequence(sample, contig_num)
        if query_record:
            query_fasta = OUT_DIR / f"query_{sample}_kpc.fasta"
            SeqIO.write(query_record, query_fasta, "fasta")
            
            blast_out = OUT_DIR / f"blast_brig_{sample}_kpc.txt"
            cmd = f"blastn -query {query_fasta} -db {OUT_DIR}/ref_kpc_db -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 1e-5 -out {blast_out}"
            subprocess.run(cmd, shell=True, capture_output=True)
            
            hits = []
            if blast_out.exists():
                with open(blast_out) as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 12:
                            hits.append({
                                'identity': float(parts[2]),
                                's_start': int(parts[8]),
                                's_end': int(parts[9]),
                                'length': int(parts[3])
                            })
            blast_results[sample] = hits
            query_fasta.unlink(missing_ok=True)
            blast_out.unlink(missing_ok=True)
    
    # Create circular figure
    circos = Circos(sectors={f"KPC Plasmid ({ref_sample} ref, {ref_len/1000:.0f}kb)": ref_len})
    sector = list(circos.sectors)[0]
    
    # Track 1: Gene annotations (outermost)
    gene_track = sector.add_track((95, 100))
    gene_track.axis(fc="#F0F0F0", ec="black", lw=0.5)
    
    for g in ref_genes:
        color, _ = classify_gene_color(g['gene'], g['product'])
        gene_track.rect(g['start'], g['stop'], fc=color, ec="none", lw=0)
    
    # Track 2: Gene labels for AMR genes
    label_track = sector.add_track((100, 105))
    for g in ref_genes:
        gene_name = g['gene'] or ''
        if any(x in gene_name.lower() for x in ['bla', 'aac', 'aph', 'rmt', 'arm', 'qnr', 'sul', 'dfr', 'mph', 'fos', 'cat', 'ble']):
            mid = (g['start'] + g['stop']) / 2
            label_track.text(gene_name, mid, size=5, orientation="vertical")
    
    # Tracks 3-9: BLAST identity rings for each sample
    ring_colors = ['#E63946', '#457B9D', '#2A9D8F', '#F4A261', '#264653', '#E9C46A', '#9B59B6']
    sample_list = [s for s in SAMPLES if s != ref_sample]
    
    for idx, sample in enumerate(sample_list):
        inner = 88 - idx * 8
        outer = inner + 6
        if inner < 30:
            break
        
        blast_track = sector.add_track((inner, outer))
        blast_track.axis(fc="#F8F8F8", ec="gray", lw=0.3)
        
        color = ring_colors[idx % len(ring_colors)]
        
        for hit in blast_results.get(sample, []):
            s_start = min(hit['s_start'], hit['s_end'])
            s_end = max(hit['s_start'], hit['s_end'])
            identity = hit['identity']
            
            # Color intensity based on identity
            alpha = max(0.3, identity / 100)
            blast_track.rect(s_start, s_end, fc=color, ec="none", lw=0, alpha=alpha)
        
        # Add sample label
        blast_track.text(sample, ref_len * 0.02, size=5, color=color, fontweight='bold')
    
    # Add tick marks
    tick_track = sector.add_track((25, 28))
    tick_positions = list(range(0, ref_len, 10000))
    for pos in tick_positions:
        tick_track.line([pos, pos], [25, 28], lw=0.3, color="gray")
    
    # Scale bar labels
    scale_track = sector.add_track((20, 24))
    for pos in range(0, ref_len, 20000):
        scale_track.text(f"{pos//1000}kb", pos, size=4, color="gray")
    
    fig = circos.plotfig()
    
    # Add title
    fig.suptitle(f'BRIG-style Circular Map: KPC-2 Plasmid (~{ref_len/1000:.0f}kb)\n'
                 f'Reference: {ref_sample} | Rings: BLAST identity of 7 other isolates',
                 fontsize=11, fontweight='bold', y=0.98)
    
    # Add legend
    legend_items = [(ring_colors[i], sample_list[i]) for i in range(min(len(sample_list), len(ring_colors)))]
    legend_items.insert(0, ('#E63946', 'Carbapenemase'))
    legend_items.insert(1, ('#F4A261', 'Beta-lactamase'))
    legend_items.insert(2, ('#9B59B6', 'Transposase/IS'))
    
    handles = [mpatches.Patch(color=c, label=l) for c, l in legend_items]
    fig.legend(handles=handles, loc='lower right', fontsize=6, framealpha=0.9, ncol=2)
    
    fig.savefig(OUT_DIR / 'BRIG_KPC_Plasmid.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / 'BRIG_KPC_Plasmid.pdf', bbox_inches='tight')
    plt.close()
    print("  BRIG KPC plasmid - DONE")
    
    # Clean up
    for f in OUT_DIR.glob("ref_kpc_db*"):
        f.unlink(missing_ok=True)
    ref_fasta.unlink(missing_ok=True)


def create_brig_ndm_plasmid():
    """Create BRIG-style circular map of the NDM plasmid (~350kb)."""
    from pycirclize import Circos
    
    print("Creating BRIG-style circular map for NDM plasmid...")
    
    ref_sample = 'KP21915'
    ref_contig_num = NDM_CONTIG[ref_sample]
    ref_record = extract_contig_sequence(ref_sample, ref_contig_num)
    ref_genes = get_genes_on_contig(ref_sample, f"contig_{ref_contig_num}")
    ref_len = len(ref_record.seq)
    
    ref_fasta = OUT_DIR / "ref_ndm_plasmid.fasta"
    SeqIO.write(ref_record, ref_fasta, "fasta")
    
    subprocess.run(f"makeblastdb -in {ref_fasta} -dbtype nucl -out {OUT_DIR}/ref_ndm_db",
                   shell=True, capture_output=True)
    
    ndm_samples = [s for s in SAMPLES if s in NDM_CONTIG and s != ref_sample]
    
    blast_results = {}
    for sample in ndm_samples:
        contig_num = NDM_CONTIG[sample]
        query_record = extract_contig_sequence(sample, contig_num)
        if query_record:
            query_fasta = OUT_DIR / f"query_{sample}_ndm.fasta"
            SeqIO.write(query_record, query_fasta, "fasta")
            
            blast_out = OUT_DIR / f"blast_brig_{sample}_ndm.txt"
            cmd = f"blastn -query {query_fasta} -db {OUT_DIR}/ref_ndm_db -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 1e-5 -out {blast_out}"
            subprocess.run(cmd, shell=True, capture_output=True)
            
            hits = []
            if blast_out.exists():
                with open(blast_out) as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 12:
                            hits.append({
                                'identity': float(parts[2]),
                                's_start': int(parts[8]),
                                's_end': int(parts[9]),
                                'length': int(parts[3])
                            })
            blast_results[sample] = hits
            query_fasta.unlink(missing_ok=True)
            blast_out.unlink(missing_ok=True)
    
    circos = Circos(sectors={f"NDM Plasmid ({ref_sample} ref, {ref_len/1000:.0f}kb)": ref_len})
    sector = list(circos.sectors)[0]
    
    gene_track = sector.add_track((95, 100))
    gene_track.axis(fc="#F0F0F0", ec="black", lw=0.5)
    
    for g in ref_genes:
        color, _ = classify_gene_color(g['gene'], g['product'])
        gene_track.rect(g['start'], g['stop'], fc=color, ec="none", lw=0)
    
    label_track = sector.add_track((100, 105))
    for g in ref_genes:
        gene_name = g['gene'] or ''
        if any(x in gene_name.lower() for x in ['bla', 'aac', 'aph', 'rmt', 'arm', 'qnr', 'sul', 'dfr', 'mph', 'msr', 'fos', 'cat', 'ble']):
            mid = (g['start'] + g['stop']) / 2
            label_track.text(gene_name, mid, size=4, orientation="vertical")
    
    ring_colors = ['#E63946', '#457B9D', '#2A9D8F', '#F4A261', '#264653']
    
    for idx, sample in enumerate(ndm_samples):
        inner = 88 - idx * 10
        outer = inner + 8
        if inner < 35:
            break
        
        blast_track = sector.add_track((inner, outer))
        blast_track.axis(fc="#F8F8F8", ec="gray", lw=0.3)
        
        color = ring_colors[idx % len(ring_colors)]
        
        for hit in blast_results.get(sample, []):
            s_start = min(hit['s_start'], hit['s_end'])
            s_end = max(hit['s_start'], hit['s_end'])
            identity = hit['identity']
            alpha = max(0.3, identity / 100)
            blast_track.rect(s_start, s_end, fc=color, ec="none", lw=0, alpha=alpha)
        
        blast_track.text(sample, ref_len * 0.02, size=5, color=color, fontweight='bold')
    
    tick_track = sector.add_track((30, 33))
    for pos in range(0, ref_len, 20000):
        tick_track.line([pos, pos], [30, 33], lw=0.3, color="gray")
    
    scale_track = sector.add_track((25, 29))
    for pos in range(0, ref_len, 50000):
        scale_track.text(f"{pos//1000}kb", pos, size=4, color="gray")
    
    fig = circos.plotfig()
    
    fig.suptitle(f'BRIG-style Circular Map: NDM-1 Plasmid (~{ref_len/1000:.0f}kb)\n'
                 f'Reference: {ref_sample} | Rings: BLAST identity of 5 other NDM+ isolates',
                 fontsize=11, fontweight='bold', y=0.98)
    
    legend_items = [(ring_colors[i], ndm_samples[i]) for i in range(min(len(ndm_samples), len(ring_colors)))]
    legend_items.insert(0, ('#E63946', 'Carbapenemase'))
    legend_items.insert(1, ('#F4A261', 'Beta-lactamase'))
    legend_items.insert(2, ('#9B59B6', 'Transposase/IS'))
    
    handles = [mpatches.Patch(color=c, label=l) for c, l in legend_items]
    fig.legend(handles=handles, loc='lower right', fontsize=6, framealpha=0.9, ncol=2)
    
    fig.savefig(OUT_DIR / 'BRIG_NDM_Plasmid.png', dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / 'BRIG_NDM_Plasmid.pdf', bbox_inches='tight')
    plt.close()
    print("  BRIG NDM plasmid - DONE")
    
    for f in OUT_DIR.glob("ref_ndm_db*"):
        f.unlink(missing_ok=True)
    ref_fasta.unlink(missing_ok=True)


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("BRIG & EasyFig-style Visualization Pipeline")
    print("=" * 60)
    
    # EasyFig-style linear comparisons
    create_easyfig_kpc_environment()
    create_easyfig_ndm_environment()
    
    # BRIG-style circular maps
    create_brig_kpc_plasmid()
    create_brig_ndm_plasmid()
    
    print("\n" + "=" * 60)
    print(f"All figures saved to: {OUT_DIR}")
    print("Files generated:")
    for f in sorted(OUT_DIR.glob('*')):
        if f.is_file() and not f.name.startswith('temp'):
            print(f"  {f.name}")
