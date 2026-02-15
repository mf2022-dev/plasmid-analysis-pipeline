#!/usr/bin/env python3
"""
High-Quality EasyFig-style Linear Comparison Maps
with ALL genes annotated for KPC-2 and NDM-1 plasmids.
Publication-ready for Lancet Microbe.
"""
import os, csv, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow, FancyBboxPatch
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as pe
import numpy as np
from pathlib import Path
from Bio import SeqIO

BASE = Path("/home/ubuntu/plasmid_analysis")
ANNOT_DIR = BASE / "04_AMR_Plasmid_Association" / "Bakta_Annotations"
LR_DIR = BASE / "LR_assemblies"
OUT_DIR = BASE / "08_BRIG_EasyFig"
OUT_DIR.mkdir(exist_ok=True)

SAMPLES = ['JKPB244', 'JKPB284', 'JKPB381', 'JKPR282', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']

# Color scheme - Lancet-quality
COLORS = {
    'carbapenemase': '#D32F2F',      # Deep red
    'beta_lactamase': '#FF6F00',     # Amber
    'aminoglycoside': '#1565C0',     # Blue
    'quinolone': '#2E7D32',          # Green
    'macrolide': '#6A1B9A',          # Purple
    'sulfonamide': '#00838F',        # Teal
    'trimethoprim': '#4E342E',       # Brown
    'phenicol': '#F57F17',           # Yellow-dark
    'fosfomycin': '#AD1457',         # Pink
    'tetracycline': '#E65100',       # Deep orange
    'other_amr': '#795548',          # Brown
    'transposase': '#78909C',        # Blue-grey
    'integrase': '#C62828',          # Dark red
    'replication': '#0277BD',        # Light blue
    'conjugation': '#00695C',        # Dark teal
    'mercury': '#546E7A',            # Grey-blue
    'hypothetical': '#E0E0E0',       # Light grey
    'other': '#BDBDBD',             # Grey
}

# AMR gene classification
def classify_gene(gene, product):
    g = (gene or '').lower()
    p = (product or '').lower()
    
    if any(x in g for x in ['blakpc', 'blandm', 'blaoxa', 'blavim', 'blaimp']):
        return 'carbapenemase', True
    if any(x in g for x in ['blactx', 'blashv', 'blatem', 'blaper']):
        return 'beta_lactamase', True
    if any(x in g for x in ['rmtb', 'arma', 'aph(', 'aac(', 'ant(', 'aada']):
        return 'aminoglycoside', True
    if any(x in g for x in ['qnr', 'oqx', 'aac(6']):
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
    if 'ble' == g:
        return 'other_amr', True
    if 'tnp' in g or 'transpos' in p:
        return 'transposase', False
    if 'integr' in p or 'recombinas' in p:
        return 'integrase', False
    if 'replicat' in p or 'partition' in p or 'parM' in g or 'parA' in g or 'parB' in g:
        return 'replication', False
    if 'conjug' in p or 'transfer' in p or 'tra' == g[:3] if g else False:
        return 'conjugation', False
    if 'mercury' in p or 'mer' == g[:3] if g else False:
        return 'mercury', False
    if 'hypothetical' in p:
        return 'hypothetical', False
    return 'other', False

def parse_tsv(tsv_path, target_contig):
    """Parse Bakta TSV and return all genes on the target contig."""
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


def get_contig_size(sample, contig_num):
    """Get contig size from LR assembly FASTA."""
    fasta = LR_DIR / f"{sample}_LR.fasta"
    if fasta.exists():
        for i, rec in enumerate(SeqIO.parse(str(fasta), 'fasta'), 1):
            if i == contig_num:
                return len(rec)
    return None


def draw_gene_arrow(ax, x_start, x_end, y, height, strand, color, alpha=1.0):
    """Draw a gene as a directional arrow."""
    gene_len = abs(x_end - x_start)
    arrow_head = min(gene_len * 0.3, 800)  # Arrow head size
    
    if strand == 1:  # Forward
        body_end = x_end - arrow_head
        if body_end < x_start:
            body_end = x_start
        # Body
        rect = plt.Rectangle((x_start, y - height/2), body_end - x_start, height,
                             facecolor=color, edgecolor='black', linewidth=0.3, alpha=alpha)
        ax.add_patch(rect)
        # Arrow head
        triangle = plt.Polygon([
            [body_end, y - height/2 - height*0.15],
            [x_end, y],
            [body_end, y + height/2 + height*0.15]
        ], facecolor=color, edgecolor='black', linewidth=0.3, alpha=alpha)
        ax.add_patch(triangle)
    else:  # Reverse
        body_start = x_start + arrow_head
        if body_start > x_end:
            body_start = x_end
        # Body
        rect = plt.Rectangle((body_start, y - height/2), x_end - body_start, height,
                             facecolor=color, edgecolor='black', linewidth=0.3, alpha=alpha)
        ax.add_patch(rect)
        # Arrow head
        triangle = plt.Polygon([
            [body_start, y - height/2 - height*0.15],
            [x_start, y],
            [body_start, y + height/2 + height*0.15]
        ], facecolor=color, edgecolor='black', linewidth=0.3, alpha=alpha)
        ax.add_patch(triangle)


def create_easyfig_plasmid(plasmid_type, contig_num, samples_list, title, output_name):
    """Create a high-quality EasyFig-style linear map with ALL genes annotated."""
    
    # Collect data for all samples
    all_data = {}
    max_len = 0
    for sample in samples_list:
        tsv = ANNOT_DIR / sample / f"{sample}.tsv"
        if not tsv.exists():
            continue
        genes = parse_tsv(str(tsv), f"contig_{contig_num}")
        if not genes:
            continue
        contig_size = get_contig_size(sample, contig_num)
        if contig_size is None:
            contig_size = max(g['stop'] for g in genes) + 1000
        all_data[sample] = {'genes': genes, 'size': contig_size}
        max_len = max(max_len, contig_size)
    
    if not all_data:
        print(f"  No data found for {plasmid_type}")
        return
    
    n_samples = len(all_data)
    track_height = 1.8  # Height per sample track
    gap = 0.6  # Gap between tracks
    label_space = 1.2  # Space for gene labels above/below
    
    fig_height = n_samples * (track_height + gap + label_space) + 3
    fig_width = 28  # Wide figure for full plasmid
    
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height), dpi=300)
    
    # Title
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20, fontfamily='serif')
    
    gene_height = 0.35
    
    sample_list = [s for s in samples_list if s in all_data]
    
    for idx, sample in enumerate(sample_list):
        data = all_data[sample]
        genes = data['genes']
        size = data['size']
        
        y_center = (n_samples - 1 - idx) * (track_height + gap)
        
        # Draw backbone line
        ax.plot([0, size], [y_center, y_center], color='#424242', linewidth=1.0, zorder=1)
        
        # Sample label
        ax.text(-max_len * 0.02, y_center, sample, fontsize=10, fontweight='bold',
                ha='right', va='center', fontfamily='sans-serif')
        
        # Draw ALL genes
        label_positions_above = []
        label_positions_below = []
        
        for g in genes:
            cat, is_amr = classify_gene(g['gene'], g['product'])
            color = COLORS.get(cat, COLORS['other'])
            alpha = 1.0 if cat != 'hypothetical' else 0.5
            
            draw_gene_arrow(ax, g['start'], g['stop'], y_center, gene_height,
                          g['strand'], color, alpha)
            
            # Determine label
            label = g['gene'] if g['gene'] else None
            if not label and cat == 'hypothetical':
                label = None  # Skip hypothetical labels
            elif not label:
                # Use abbreviated product name
                prod = g['product']
                if len(prod) > 20:
                    # Try to get a short name
                    words = prod.split()
                    if len(words) > 2:
                        label = None  # Too long, skip
                    else:
                        label = prod[:18]
                else:
                    label = prod
            
            if label:
                gene_mid = (g['start'] + g['stop']) / 2
                fontsize = 5.5 if is_amr else 4.5
                fontweight = 'bold' if is_amr else 'normal'
                fontstyle = 'italic' if is_amr else 'normal'
                
                # Alternate labels above and below to avoid overlap
                if g['strand'] == 1:
                    y_label = y_center + gene_height/2 + 0.15
                    va = 'bottom'
                    label_positions_above.append((gene_mid, label))
                else:
                    y_label = y_center - gene_height/2 - 0.15
                    va = 'top'
                    label_positions_below.append((gene_mid, label))
                
                text_color = '#B71C1C' if is_amr else '#212121'
                
                ax.text(gene_mid, y_label, label, fontsize=fontsize,
                       ha='center', va=va, rotation=45,
                       fontweight=fontweight, fontstyle=fontstyle,
                       color=text_color, zorder=10,
                       path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
        
        # Scale markers every 10kb
        for pos in range(0, int(size) + 1, 10000):
            ax.plot([pos, pos], [y_center - 0.08, y_center + 0.08], color='#616161', linewidth=0.5)
            if idx == n_samples - 1:  # Only on bottom track
                ax.text(pos, y_center - gene_height - 0.8, f'{pos//1000}kb',
                       fontsize=5, ha='center', va='top', color='#616161')
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLORS['carbapenemase'], edgecolor='black', linewidth=0.5, label='Carbapenemase'),
        mpatches.Patch(facecolor=COLORS['beta_lactamase'], edgecolor='black', linewidth=0.5, label='Beta-lactamase'),
        mpatches.Patch(facecolor=COLORS['aminoglycoside'], edgecolor='black', linewidth=0.5, label='Aminoglycoside resistance'),
        mpatches.Patch(facecolor=COLORS['quinolone'], edgecolor='black', linewidth=0.5, label='Quinolone resistance'),
        mpatches.Patch(facecolor=COLORS['macrolide'], edgecolor='black', linewidth=0.5, label='Macrolide resistance'),
        mpatches.Patch(facecolor=COLORS['sulfonamide'], edgecolor='black', linewidth=0.5, label='Sulfonamide resistance'),
        mpatches.Patch(facecolor=COLORS['trimethoprim'], edgecolor='black', linewidth=0.5, label='Trimethoprim resistance'),
        mpatches.Patch(facecolor=COLORS['phenicol'], edgecolor='black', linewidth=0.5, label='Phenicol resistance'),
        mpatches.Patch(facecolor=COLORS['fosfomycin'], edgecolor='black', linewidth=0.5, label='Fosfomycin resistance'),
        mpatches.Patch(facecolor=COLORS['other_amr'], edgecolor='black', linewidth=0.5, label='Other AMR'),
        mpatches.Patch(facecolor=COLORS['transposase'], edgecolor='black', linewidth=0.5, label='Transposase/IS element'),
        mpatches.Patch(facecolor=COLORS['integrase'], edgecolor='black', linewidth=0.5, label='Integrase/Recombinase'),
        mpatches.Patch(facecolor=COLORS['replication'], edgecolor='black', linewidth=0.5, label='Replication/Partition'),
        mpatches.Patch(facecolor=COLORS['conjugation'], edgecolor='black', linewidth=0.5, label='Conjugation/Transfer'),
        mpatches.Patch(facecolor=COLORS['mercury'], edgecolor='black', linewidth=0.5, label='Mercury resistance'),
        mpatches.Patch(facecolor=COLORS['hypothetical'], edgecolor='black', linewidth=0.5, label='Hypothetical protein'),
        mpatches.Patch(facecolor=COLORS['other'], edgecolor='black', linewidth=0.5, label='Other function'),
    ]
    
    legend = ax.legend(handles=legend_elements, loc='lower center',
                      bbox_to_anchor=(0.5, -0.08), ncol=6, fontsize=6,
                      frameon=True, fancybox=True, shadow=False,
                      edgecolor='#BDBDBD', facecolor='white')
    
    ax.set_xlim(-max_len * 0.08, max_len * 1.02)
    y_min = -gap - label_space
    y_max = (n_samples - 1) * (track_height + gap) + track_height + label_space
    ax.set_ylim(y_min, y_max)
    ax.set_aspect('auto')
    ax.axis('off')
    
    plt.tight_layout()
    
    # Save
    for ext in ['png', 'pdf']:
        outpath = OUT_DIR / f"{output_name}.{ext}"
        fig.savefig(str(outpath), dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved: {output_name}.png/.pdf")


# ============================================================
# MAIN
# ============================================================
print("=" * 60)
print("High-Quality EasyFig-style Visualization Pipeline")
print("=" * 60)

# KPC plasmid (contig_3, ~131kb)
print("\n1. Creating KPC-2 plasmid linear map (all genes)...")
create_easyfig_plasmid(
    'KPC', 3, SAMPLES,
    'Genomic Organization of blaKPC-2-carrying IncFII+IncR Plasmid (~131 kb)\nLinear Comparison Across 8 ST11 K. pneumoniae Isolates',
    'HQ_EasyFig_KPC_Plasmid'
)

# NDM plasmid (contig_2, ~354kb) - only NDM+ isolates
NDM_SAMPLES = ['JKPB244', 'JKPB284', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']
print("\n2. Creating NDM-1 plasmid linear map (all genes)...")
create_easyfig_plasmid(
    'NDM', 2, NDM_SAMPLES,
    'Genomic Organization of blaNDM-1-carrying IncFIB(pNDM-Mar)+IncHI1B Plasmid (~354 kb)\nLinear Comparison Across 6 ST11 K. pneumoniae Isolates',
    'HQ_EasyFig_NDM_Plasmid'
)

print("\nDone!")
