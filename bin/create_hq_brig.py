#!/usr/bin/env python3
"""
High-Quality BRIG-style Circular Plasmid Maps
with ALL genes annotated around the outer ring.
Publication-ready for Lancet Microbe.
Uses pyCirclize for circular genome visualization.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
from pycirclize import Circos
from pycirclize.parser import Genbank
from pathlib import Path
from Bio import SeqIO
import subprocess, tempfile

BASE = Path("/home/ubuntu/plasmid_analysis")
ANNOT_DIR = BASE / "04_AMR_Plasmid_Association" / "Bakta_Annotations"
LR_DIR = BASE / "LR_assemblies"
OUT_DIR = BASE / "08_BRIG_EasyFig"
OUT_DIR.mkdir(exist_ok=True)

SAMPLES = ['JKPB244', 'JKPB284', 'JKPB381', 'JKPR282', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']
NDM_SAMPLES = ['JKPB244', 'JKPB284', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']

# Color scheme
COLORS = {
    'carbapenemase': '#D32F2F',
    'beta_lactamase': '#FF6F00',
    'aminoglycoside': '#1565C0',
    'quinolone': '#2E7D32',
    'macrolide': '#6A1B9A',
    'sulfonamide': '#00838F',
    'trimethoprim': '#4E342E',
    'phenicol': '#F57F17',
    'fosfomycin': '#AD1457',
    'tetracycline': '#E65100',
    'bleomycin': '#795548',
    'other_amr': '#795548',
    'transposase': '#90A4AE',
    'integrase': '#EF5350',
    'replication': '#42A5F5',
    'conjugation': '#26A69A',
    'mercury': '#78909C',
    'hypothetical': '#E0E0E0',
    'other': '#CFD8DC',
}

# BLAST identity color gradient
IDENTITY_COLORS = [
    (100, '#2E7D32'),  # Dark green
    (95, '#4CAF50'),   # Green
    (90, '#8BC34A'),   # Light green
    (85, '#CDDC39'),   # Lime
    (80, '#FFEB3B'),   # Yellow
    (70, '#FF9800'),   # Orange
    (0, '#F44336'),    # Red
]

def get_identity_color(identity):
    for threshold, color in IDENTITY_COLORS:
        if identity >= threshold:
            return color
    return '#F44336'

def classify_gene(gene, product):
    g = (gene or '').lower()
    p = (product or '').lower()
    
    if any(x in g for x in ['blakpc', 'blandm', 'blaoxa', 'blavim', 'blaimp']):
        return 'carbapenemase', True
    if any(x in g for x in ['blactx', 'blashv', 'blatem', 'blaper']):
        return 'beta_lactamase', True
    if any(x in g for x in ['rmtb', 'arma', 'aph(', 'aac(', 'ant(', 'aada']):
        return 'aminoglycoside', True
    if any(x in g for x in ['qnr', 'oqx']):
        return 'quinolone', True
    if any(x in g for x in ['mph(', 'msr(', 'ere(']):
        return 'macrolide', True
    if any(x in g for x in ['sul1', 'sul2', 'sul3']):
        return 'sulfonamide', True
    if any(x in g for x in ['dfra', 'dfrb']):
        return 'trimethoprim', True
    if any(x in g for x in ['cata', 'cmla', 'flor']):
        return 'phenicol', True
    if any(x in g for x in ['fosa', 'fosb']):
        return 'fosfomycin', True
    if any(x in g for x in ['tet(']):
        return 'tetracycline', True
    if g == 'ble':
        return 'bleomycin', True
    if 'tnp' in g or 'transpos' in p:
        return 'transposase', False
    if 'integr' in p or 'recombinas' in p or 'resolvase' in p.lower():
        return 'integrase', False
    if any(x in p for x in ['replicat', 'partition', 'plasmid stable', 'plasmid segreg']):
        return 'replication', False
    if any(x in p for x in ['conjug', 'transfer', 'fertility', 'pilus']):
        return 'conjugation', False
    if 'mercury' in p or (g and g.startswith('mer') and len(g) <= 4):
        return 'mercury', False
    if 'hypothetical' in p:
        return 'hypothetical', False
    return 'other', False

def parse_tsv(tsv_path, target_contig):
    genes = []
    with open(tsv_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            seq_id, feat_type, start, stop, strand, locus, gene, product = parts[:8]
            if seq_id != target_contig:
                continue
            genes.append({
                'start': int(start),
                'stop': int(stop),
                'strand': 1 if strand == '+' else -1,
                'gene': gene if gene else None,
                'product': product,
                'type': feat_type,
            })
    return sorted(genes, key=lambda x: x['start'])

# Map sample names to file prefixes
SAMPLE_FILE_MAP = {
    'JKPB244': '1_JKPB244_LR.fasta',
    'JKPB284': '2_JKPB284_LR.fasta',
    'JKPB381': '4_JKPB381_LR.fasta',
    'JKPR282': '5_JKPR282_LR.fasta',
    'KP21802': '6_KP21802_LR.fasta',
    'KP21915': '7_KP21915_LR.fasta',
    'MDRKP121': '8_MDRKP121_LR.fasta',
    'MDRKP122': '9_MDRKP122_LR.fasta',
}

def get_contig_seq(sample, contig_num):
    """Get contig sequence. Contig IDs in FASTA are just numbers (1, 2, 3...)."""
    fname = SAMPLE_FILE_MAP.get(sample)
    if not fname:
        return None
    fasta = LR_DIR / fname
    if fasta.exists():
        for rec in SeqIO.parse(str(fasta), 'fasta'):
            # Contig IDs are just numbers: '1', '2', '3', etc.
            if rec.id == str(contig_num):
                return rec
    return None

def run_blast_comparison(ref_seq, query_seq):
    """Run BLAST between reference and query sequences, return identity segments."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as ref_f:
        ref_f.write(f">ref\n{str(ref_seq.seq)}\n")
        ref_path = ref_f.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as qry_f:
        qry_f.write(f">query\n{str(query_seq.seq)}\n")
        qry_path = qry_f.name
    
    # Make BLAST database
    subprocess.run(['makeblastdb', '-in', ref_path, '-dbtype', 'nucl', '-out', ref_path + '.db'],
                  capture_output=True)
    
    # Run BLAST
    result = subprocess.run(
        ['blastn', '-query', qry_path, '-db', ref_path + '.db',
         '-outfmt', '6 qstart qend sstart send pident length',
         '-evalue', '1e-10', '-max_target_seqs', '1'],
        capture_output=True, text=True
    )
    
    segments = []
    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 6:
                sstart, send = int(parts[2]), int(parts[3])
                pident = float(parts[4])
                segments.append((min(sstart, send), max(sstart, send), pident))
    
    # Clean up
    for f in [ref_path, qry_path, ref_path + '.db.ndb', ref_path + '.db.nhr',
              ref_path + '.db.nin', ref_path + '.db.not', ref_path + '.db.nsq',
              ref_path + '.db.ntf', ref_path + '.db.nto']:
        try:
            os.unlink(f)
        except:
            pass
    
    return segments


def create_brig_map(contig_num, ref_sample, comparison_samples, title, output_name, plasmid_size_label):
    """Create a BRIG-style circular map with gene annotations and BLAST rings."""
    
    print(f"  Getting reference sequence from {ref_sample}...")
    ref_seq = get_contig_seq(ref_sample, contig_num)
    if ref_seq is None:
        print(f"  ERROR: No sequence found for {ref_sample} contig_{contig_num}")
        return
    
    ref_len = len(ref_seq)
    ref_genes = parse_tsv(str(ANNOT_DIR / ref_sample / f"{ref_sample}.tsv"), f"contig_{contig_num}")
    
    print(f"  Reference: {ref_sample}, {ref_len:,} bp, {len(ref_genes)} genes")
    
    # Create Circos plot
    sectors = {f"pKPC ({plasmid_size_label})": ref_len}
    circos = Circos(sectors, space=5)
    
    sector = circos.sectors[0]
    
    # Track 1: Outer gene annotation track (forward strand)
    fwd_track = sector.add_track((95, 100))
    fwd_track.axis(fc='#FAFAFA', ec='#BDBDBD', lw=0.3)
    
    # Track 2: Outer gene annotation track (reverse strand)
    rev_track = sector.add_track((89, 94))
    rev_track.axis(fc='#FAFAFA', ec='#BDBDBD', lw=0.3)
    
    # Draw genes on forward and reverse tracks
    amr_label_positions = []
    
    for g in ref_genes:
        cat, is_amr = classify_gene(g['gene'], g['product'])
        color = COLORS.get(cat, COLORS['other'])
        
        if cat == 'hypothetical':
            alpha = 0.3
        elif is_amr:
            alpha = 1.0
        else:
            alpha = 0.7
        
        if g['strand'] == 1:
            fwd_track.rect(g['start'], g['stop'], fc=color, ec='#424242', lw=0.2, alpha=alpha)
        else:
            rev_track.rect(g['start'], g['stop'], fc=color, ec='#424242', lw=0.2, alpha=alpha)
        
        # Collect AMR gene labels
        label = g['gene'] if g['gene'] else None
        if label and cat != 'hypothetical':
            gene_mid = (g['start'] + g['stop']) / 2
            amr_label_positions.append({
                'pos': gene_mid,
                'label': label,
                'is_amr': is_amr,
                'cat': cat,
                'strand': g['strand'],
            })
    
    # Track 3-N: BLAST comparison rings
    ring_colors = ['#1976D2', '#388E3C', '#F57C00', '#7B1FA2', '#C62828', '#00838F', '#5D4037', '#455A64']
    
    comp_samples = [s for s in comparison_samples if s != ref_sample]
    
    for i, sample in enumerate(comp_samples):
        print(f"  Running BLAST: {ref_sample} vs {sample}...")
        query_seq = get_contig_seq(sample, contig_num)
        if query_seq is None:
            print(f"    Skipping {sample} - no sequence")
            continue
        
        segments = run_blast_comparison(ref_seq, query_seq)
        
        r_outer = 87 - i * 8
        r_inner = r_outer - 6
        
        blast_track = sector.add_track((r_inner, r_outer))
        blast_track.axis(fc='#F5F5F5', ec='#E0E0E0', lw=0.2)
        
        # Draw BLAST hits colored by identity
        for sstart, send, pident in segments:
            color = get_identity_color(pident)
            blast_track.rect(sstart, send, fc=color, ec='none', alpha=0.85)
        
        # Sample label on the ring
        blast_track.text(sample, x=ref_len * 0.98, fontsize=5, color='#212121',
                        adjust_rotation=True)
    
    # GC content track (innermost)
    gc_r_outer = 87 - len(comp_samples) * 8 - 2
    gc_r_inner = gc_r_outer - 8
    
    if gc_r_inner > 15:
        gc_track = sector.add_track((gc_r_inner, gc_r_outer))
        gc_track.axis(fc='#FAFAFA', ec='#E0E0E0', lw=0.2)
        
        # Calculate GC content in windows
        window_size = max(ref_len // 500, 100)
        seq_str = str(ref_seq.seq).upper()
        gc_values = []
        positions = []
        for pos in range(0, ref_len - window_size, window_size):
            window = seq_str[pos:pos + window_size]
            gc = (window.count('G') + window.count('C')) / len(window) * 100
            gc_values.append(gc)
            positions.append(pos + window_size // 2)
        
        if gc_values:
            gc_mean = np.mean(gc_values)
            # Plot GC content as line
            vmin, vmax = min(gc_values), max(gc_values)
            for j in range(len(positions) - 1):
                gc_val = gc_values[j]
                if gc_val >= gc_mean:
                    color = '#2E7D32'
                else:
                    color = '#C62828'
                gc_track.rect(positions[j], positions[j+1], fc=color, ec='none', alpha=0.5)
    
    # Add gene labels around the outside
    label_track = sector.add_track((101, 115))
    
    # Filter and space labels to avoid overlap
    min_spacing = ref_len * 0.015  # Minimum spacing between labels
    
    # Sort: AMR genes first (priority), then by position
    sorted_labels = sorted(amr_label_positions, key=lambda x: (-x['is_amr'], x['pos']))
    placed_positions = []
    
    for lb in sorted_labels:
        pos = lb['pos']
        label = lb['label']
        is_amr = lb['is_amr']
        cat = lb['cat']
        
        # Check overlap
        too_close = False
        for pp in placed_positions:
            if abs(pos - pp) < min_spacing:
                too_close = True
                break
        
        if too_close and not is_amr:
            continue  # Skip non-AMR if overlapping
        
        if is_amr:
            fontsize = 7
            fontweight = 'bold'
            color = COLORS.get(cat, '#D32F2F')
        else:
            fontsize = 4.5
            fontweight = 'normal'
            color = '#616161'
        
        label_track.text(label, x=pos, fontsize=fontsize, color=color,
                        r=108, adjust_rotation=True, fontweight=fontweight,
                        fontstyle='italic')
        placed_positions.append(pos)
    
    # Add tick marks
    sector.add_track((100, 101)).axis(fc='none', ec='none')
    
    # Scale ticks every 10kb
    for pos in range(0, ref_len, 10000):
        sector.add_track((100, 101))
    
    # Create figure
    fig = circos.plotfig(figsize=(16, 16), dpi=300)
    
    # Add title
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98, fontfamily='serif')
    
    # Add legend
    legend_items = [
        ('Carbapenemase', COLORS['carbapenemase']),
        ('Beta-lactamase', COLORS['beta_lactamase']),
        ('Aminoglycoside', COLORS['aminoglycoside']),
        ('Quinolone', COLORS['quinolone']),
        ('Macrolide', COLORS['macrolide']),
        ('Sulfonamide', COLORS['sulfonamide']),
        ('Trimethoprim', COLORS['trimethoprim']),
        ('Phenicol', COLORS['phenicol']),
        ('Fosfomycin', COLORS['fosfomycin']),
        ('Transposase/IS', COLORS['transposase']),
        ('Integrase', COLORS['integrase']),
        ('Replication', COLORS['replication']),
        ('Conjugation', COLORS['conjugation']),
        ('Mercury res.', COLORS['mercury']),
        ('Hypothetical', COLORS['hypothetical']),
    ]
    
    # Identity legend
    identity_items = [
        ('>=100%', '#2E7D32'),
        ('>=95%', '#4CAF50'),
        ('>=90%', '#8BC34A'),
        ('>=85%', '#CDDC39'),
        ('>=80%', '#FFEB3B'),
        ('>=70%', '#FF9800'),
        ('<70%', '#F44336'),
    ]
    
    handles1 = [mpatches.Patch(facecolor=c, edgecolor='#424242', lw=0.5, label=n) for n, c in legend_items]
    handles2 = [mpatches.Patch(facecolor=c, edgecolor='#424242', lw=0.5, label=n) for n, c in identity_items]
    
    leg1 = fig.legend(handles=handles1, loc='lower left', bbox_to_anchor=(0.02, 0.02),
                     ncol=3, fontsize=6.5, title='Gene Categories', title_fontsize=8,
                     frameon=True, fancybox=True, edgecolor='#BDBDBD')
    
    leg2 = fig.legend(handles=handles2, loc='lower right', bbox_to_anchor=(0.98, 0.02),
                     ncol=2, fontsize=6.5, title='BLAST Identity', title_fontsize=8,
                     frameon=True, fancybox=True, edgecolor='#BDBDBD')
    
    # Ring labels
    fig.text(0.02, 0.5, 'Outer rings: Forward/Reverse CDS\nInner rings: BLAST comparison\nInnermost: GC content',
            fontsize=7, va='center', color='#757575', fontstyle='italic')
    
    for ext in ['png', 'pdf']:
        outpath = OUT_DIR / f"{output_name}.{ext}"
        fig.savefig(str(outpath), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_name}")


# ============================================================
# MAIN
# ============================================================
print("=" * 60)
print("High-Quality BRIG-style Circular Maps")
print("=" * 60)

# 1. KPC plasmid BRIG map
print("\n1. Creating KPC plasmid BRIG map...")
create_brig_map(
    contig_num=3,
    ref_sample='KP21915',
    comparison_samples=SAMPLES,
    title='BRIG Map: blaKPC-2-carrying IncFII+IncR Plasmid (~131 kb)\n'
          'Reference: KP21915 | Comparison: 7 ST11 K. pneumoniae Isolates',
    output_name='HQ_BRIG_KPC_Plasmid',
    plasmid_size_label='~131 kb'
)

# 2. NDM plasmid BRIG map
print("\n2. Creating NDM plasmid BRIG map...")
create_brig_map(
    contig_num=2,
    ref_sample='KP21915',
    comparison_samples=NDM_SAMPLES,
    title='BRIG Map: blaNDM-1-carrying IncFIB(pNDM-Mar)+IncHI1B Plasmid (~354 kb)\n'
          'Reference: KP21915 | Comparison: 5 ST11 K. pneumoniae Isolates',
    output_name='HQ_BRIG_NDM_Plasmid',
    plasmid_size_label='~354 kb'
)

print("\n" + "=" * 60)
print("All BRIG figures complete!")
print("=" * 60)
