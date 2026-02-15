#!/usr/bin/env python3
"""
High-Quality EasyFig-style Linear Comparison Maps v2
Publication-ready for Lancet Microbe.
- ALL genes shown as arrows
- AMR genes prominently labeled with colored boxes
- Zoomed-in panels for key resistance regions
- Clean, professional design
"""
import os, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
NDM_SAMPLES = ['JKPB244', 'JKPB284', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']

# Professional color palette
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
    'hypothetical': '#E8E8E8',
    'other': '#CFD8DC',
}

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

def get_contig_size(sample, contig_num):
    fasta = LR_DIR / f"{sample}_LR.fasta"
    if fasta.exists():
        for i, rec in enumerate(SeqIO.parse(str(fasta), 'fasta'), 1):
            if i == contig_num:
                return len(rec)
    return None

def draw_gene_arrow(ax, x_start, x_end, y, height, strand, color, alpha=1.0, lw=0.3):
    """Draw a gene as a directional arrow with proper proportions."""
    gene_len = abs(x_end - x_start)
    arrow_head = min(gene_len * 0.25, 1500)
    
    if strand == 1:
        body_end = max(x_end - arrow_head, x_start)
        rect = plt.Rectangle((x_start, y - height/2), body_end - x_start, height,
                             facecolor=color, edgecolor='#424242', linewidth=lw, alpha=alpha, zorder=3)
        ax.add_patch(rect)
        triangle = plt.Polygon([
            [body_end, y - height*0.65],
            [x_end, y],
            [body_end, y + height*0.65]
        ], facecolor=color, edgecolor='#424242', linewidth=lw, alpha=alpha, zorder=3)
        ax.add_patch(triangle)
    else:
        body_start = min(x_start + arrow_head, x_end)
        rect = plt.Rectangle((body_start, y - height/2), x_end - body_start, height,
                             facecolor=color, edgecolor='#424242', linewidth=lw, alpha=alpha, zorder=3)
        ax.add_patch(rect)
        triangle = plt.Polygon([
            [body_start, y - height*0.65],
            [x_start, y],
            [body_start, y + height*0.65]
        ], facecolor=color, edgecolor='#424242', linewidth=lw, alpha=alpha, zorder=3)
        ax.add_patch(triangle)


def create_full_plasmid_map(contig_num, samples_list, title, output_name):
    """Create full plasmid map with all genes and prominent AMR labels."""
    
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
        print(f"  No data found")
        return
    
    n_samples = len(all_data)
    track_spacing = 3.0
    gene_height = 0.45
    
    fig_height = n_samples * track_spacing + 4
    fig_width = 32
    
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height), dpi=300)
    fig.patch.set_facecolor('white')
    
    ax.set_title(title, fontsize=16, fontweight='bold', pad=25,
                fontfamily='serif', color='#212121')
    
    sample_list = [s for s in samples_list if s in all_data]
    
    for idx, sample in enumerate(sample_list):
        data = all_data[sample]
        genes = data['genes']
        size = data['size']
        
        y_center = (n_samples - 1 - idx) * track_spacing
        
        # Backbone
        ax.plot([0, size], [y_center, y_center], color='#9E9E9E', linewidth=0.8, zorder=1)
        
        # Sample label with box
        ax.text(-max_len * 0.015, y_center, sample, fontsize=11, fontweight='bold',
                ha='right', va='center', fontfamily='sans-serif',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5', edgecolor='#BDBDBD', linewidth=0.5))
        
        # Collect AMR gene positions for prominent labeling
        amr_labels = []
        
        for g in genes:
            cat, is_amr = classify_gene(g['gene'], g['product'])
            color = COLORS.get(cat, COLORS['other'])
            
            if cat == 'hypothetical':
                alpha = 0.35
                lw = 0.15
                h = gene_height * 0.6
            elif is_amr:
                alpha = 1.0
                lw = 0.5
                h = gene_height * 1.2
            elif cat in ['transposase', 'integrase']:
                alpha = 0.7
                lw = 0.3
                h = gene_height * 0.85
            else:
                alpha = 0.65
                lw = 0.25
                h = gene_height * 0.8
            
            draw_gene_arrow(ax, g['start'], g['stop'], y_center, h,
                          g['strand'], color, alpha, lw)
            
            # Collect labels
            label = g['gene'] if g['gene'] else None
            if label and cat != 'hypothetical':
                gene_mid = (g['start'] + g['stop']) / 2
                amr_labels.append({
                    'x': gene_mid,
                    'label': label,
                    'is_amr': is_amr,
                    'cat': cat,
                    'strand': g['strand'],
                })
        
        # Draw labels - AMR genes get prominent labels, others get small text
        used_positions = []
        
        # First pass: AMR genes (prominent)
        for lb in sorted(amr_labels, key=lambda x: -x['is_amr']):
            x = lb['x']
            label = lb['label']
            is_amr = lb['is_amr']
            cat = lb['cat']
            strand = lb['strand']
            
            # Check for overlap with existing labels
            too_close = False
            for ux in used_positions:
                if abs(x - ux) < max_len * 0.015:
                    too_close = True
                    break
            
            if is_amr:
                # AMR genes: large bold labels with colored background
                y_offset = gene_height * 1.8 if strand == 1 else -gene_height * 1.8
                y_label = y_center + y_offset
                va = 'bottom' if strand == 1 else 'top'
                
                color = COLORS.get(cat, '#795548')
                
                ax.annotate(label, xy=(x, y_center), xytext=(x, y_label),
                           fontsize=7.5, fontweight='bold', fontstyle='italic',
                           ha='center', va=va, color='white', zorder=15,
                           bbox=dict(boxstyle='round,pad=0.2', facecolor=color,
                                    edgecolor='#424242', linewidth=0.5, alpha=0.95),
                           arrowprops=dict(arrowstyle='-', color=color, lw=0.8))
                
                used_positions.append(x)
            elif not too_close and cat not in ['hypothetical']:
                # Non-AMR named genes: small italic labels
                y_offset = gene_height * 1.3 if strand == 1 else -gene_height * 1.3
                y_label = y_center + y_offset
                va = 'bottom' if strand == 1 else 'top'
                
                ax.text(x, y_label, label, fontsize=4.0, ha='center', va=va,
                       rotation=50, color='#616161', fontstyle='italic',
                       path_effects=[pe.withStroke(linewidth=1.0, foreground='white')],
                       zorder=8)
                used_positions.append(x)
        
        # Scale bar ticks
        for pos in range(0, int(size) + 1, 10000):
            ax.plot([pos, pos], [y_center - 0.1, y_center + 0.1],
                   color='#9E9E9E', linewidth=0.4, zorder=2)
    
    # Bottom scale bar
    bottom_y = -track_spacing * 0.5
    for pos in range(0, int(max_len) + 1, 10000):
        ax.text(pos, bottom_y, f'{pos//1000}', fontsize=5, ha='center', va='top',
               color='#757575')
    ax.text(max_len / 2, bottom_y - 0.6, 'Position (kb)', fontsize=8, ha='center',
           va='top', color='#424242')
    
    # Legend
    legend_items = [
        ('Carbapenemase', 'carbapenemase'),
        ('Beta-lactamase', 'beta_lactamase'),
        ('Aminoglycoside', 'aminoglycoside'),
        ('Quinolone', 'quinolone'),
        ('Macrolide', 'macrolide'),
        ('Sulfonamide', 'sulfonamide'),
        ('Trimethoprim', 'trimethoprim'),
        ('Phenicol', 'phenicol'),
        ('Fosfomycin', 'fosfomycin'),
        ('Bleomycin', 'bleomycin'),
        ('Transposase/IS', 'transposase'),
        ('Integrase/Recombinase', 'integrase'),
        ('Replication/Partition', 'replication'),
        ('Conjugation/Transfer', 'conjugation'),
        ('Mercury resistance', 'mercury'),
        ('Hypothetical', 'hypothetical'),
        ('Other', 'other'),
    ]
    
    handles = [mpatches.Patch(facecolor=COLORS[k], edgecolor='#424242', linewidth=0.5, label=n)
               for n, k in legend_items]
    
    legend = ax.legend(handles=handles, loc='lower center',
                      bbox_to_anchor=(0.5, -0.06), ncol=6, fontsize=7,
                      frameon=True, fancybox=True, shadow=False,
                      edgecolor='#BDBDBD', facecolor='white',
                      handlelength=1.2, handleheight=0.8)
    
    ax.set_xlim(-max_len * 0.07, max_len * 1.02)
    y_min = bottom_y - 1.5
    y_max = (n_samples - 1) * track_spacing + track_spacing
    ax.set_ylim(y_min, y_max)
    ax.axis('off')
    
    plt.tight_layout()
    
    for ext in ['png', 'pdf']:
        outpath = OUT_DIR / f"{output_name}.{ext}"
        fig.savefig(str(outpath), dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved: {output_name}")


def create_zoomed_region(contig_num, samples_list, region_start, region_end,
                         title, output_name):
    """Create a zoomed-in view of a specific resistance region with ALL genes labeled."""
    
    all_data = {}
    for sample in samples_list:
        tsv = ANNOT_DIR / sample / f"{sample}.tsv"
        if not tsv.exists():
            continue
        genes = parse_tsv(str(tsv), f"contig_{contig_num}")
        # Filter to region
        region_genes = [g for g in genes 
                       if g['stop'] >= region_start and g['start'] <= region_end]
        if region_genes:
            all_data[sample] = region_genes
    
    if not all_data:
        print(f"  No data for zoomed region")
        return
    
    n_samples = len(all_data)
    region_len = region_end - region_start
    track_spacing = 2.8
    gene_height = 0.5
    
    fig_height = n_samples * track_spacing + 3.5
    fig_width = 24
    
    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height), dpi=300)
    fig.patch.set_facecolor('white')
    
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20,
                fontfamily='serif', color='#212121')
    
    sample_list = [s for s in samples_list if s in all_data]
    
    for idx, sample in enumerate(sample_list):
        genes = all_data[sample]
        y_center = (n_samples - 1 - idx) * track_spacing
        
        # Backbone
        ax.plot([region_start, region_end], [y_center, y_center],
               color='#9E9E9E', linewidth=1.0, zorder=1)
        
        # Sample label
        ax.text(region_start - region_len * 0.02, y_center, sample,
               fontsize=11, fontweight='bold', ha='right', va='center',
               fontfamily='sans-serif',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5',
                        edgecolor='#BDBDBD', linewidth=0.5))
        
        used_above = []
        used_below = []
        
        for g in genes:
            cat, is_amr = classify_gene(g['gene'], g['product'])
            color = COLORS.get(cat, COLORS['other'])
            
            if cat == 'hypothetical':
                alpha = 0.4
                h = gene_height * 0.65
                lw = 0.2
            elif is_amr:
                alpha = 1.0
                h = gene_height * 1.3
                lw = 0.6
            else:
                alpha = 0.75
                h = gene_height * 0.9
                lw = 0.35
            
            draw_gene_arrow(ax, g['start'], g['stop'], y_center, h,
                          g['strand'], color, alpha, lw)
            
            # Label ALL genes in zoomed view
            label = g['gene'] if g['gene'] else None
            if not label and cat != 'hypothetical':
                # Use short product name
                prod = g['product']
                if len(prod) <= 25:
                    label = prod
            
            if label:
                gene_mid = (g['start'] + g['stop']) / 2
                
                if is_amr:
                    # Prominent AMR label
                    y_offset = gene_height * 2.2 if g['strand'] == 1 else -gene_height * 2.2
                    y_label = y_center + y_offset
                    va = 'bottom' if g['strand'] == 1 else 'top'
                    
                    ax.annotate(label, xy=(gene_mid, y_center),
                               xytext=(gene_mid, y_label),
                               fontsize=9, fontweight='bold', fontstyle='italic',
                               ha='center', va=va, color='white', zorder=15,
                               bbox=dict(boxstyle='round,pad=0.25', facecolor=COLORS.get(cat, '#795548'),
                                        edgecolor='#212121', linewidth=0.6, alpha=0.95),
                               arrowprops=dict(arrowstyle='->', color=COLORS.get(cat, '#795548'),
                                             lw=1.0, connectionstyle='arc3,rad=0.1'))
                else:
                    # Regular gene label
                    if g['strand'] == 1:
                        # Check overlap above
                        y_offset = gene_height * 1.5
                        level = 0
                        for ux in used_above:
                            if abs(gene_mid - ux) < region_len * 0.03:
                                level += 1
                        y_label = y_center + y_offset + level * 0.5
                        va = 'bottom'
                        used_above.append(gene_mid)
                    else:
                        y_offset = -gene_height * 1.5
                        level = 0
                        for ux in used_below:
                            if abs(gene_mid - ux) < region_len * 0.03:
                                level += 1
                        y_label = y_center + y_offset - level * 0.5
                        va = 'top'
                        used_below.append(gene_mid)
                    
                    fontsize = 6.5 if cat not in ['hypothetical'] else 5
                    ax.text(gene_mid, y_label, label, fontsize=fontsize,
                           ha='center', va=va, rotation=40,
                           color='#424242', fontstyle='italic',
                           path_effects=[pe.withStroke(linewidth=1.2, foreground='white')],
                           zorder=10)
        
        # Scale ticks every 2kb in zoomed view
        for pos in range(int(region_start // 2000) * 2000, int(region_end) + 1, 2000):
            if region_start <= pos <= region_end:
                ax.plot([pos, pos], [y_center - 0.08, y_center + 0.08],
                       color='#9E9E9E', linewidth=0.4, zorder=2)
    
    # Bottom scale
    bottom_y = -track_spacing * 0.4
    for pos in range(int(region_start // 2000) * 2000, int(region_end) + 1, 2000):
        if region_start <= pos <= region_end:
            ax.text(pos, bottom_y, f'{pos//1000}kb', fontsize=6, ha='center',
                   va='top', color='#757575')
    
    # Legend
    legend_items = [
        ('Carbapenemase', 'carbapenemase'),
        ('Beta-lactamase', 'beta_lactamase'),
        ('Aminoglycoside', 'aminoglycoside'),
        ('Quinolone', 'quinolone'),
        ('Macrolide', 'macrolide'),
        ('Sulfonamide', 'sulfonamide'),
        ('Trimethoprim', 'trimethoprim'),
        ('Phenicol', 'phenicol'),
        ('Fosfomycin', 'fosfomycin'),
        ('Bleomycin', 'bleomycin'),
        ('Transposase/IS', 'transposase'),
        ('Integrase/Recombinase', 'integrase'),
        ('Mercury resistance', 'mercury'),
        ('Other', 'other'),
    ]
    handles = [mpatches.Patch(facecolor=COLORS[k], edgecolor='#424242', linewidth=0.5, label=n)
               for n, k in legend_items]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.06),
             ncol=5, fontsize=7, frameon=True, fancybox=True,
             edgecolor='#BDBDBD', facecolor='white')
    
    ax.set_xlim(region_start - region_len * 0.06, region_end + region_len * 0.02)
    y_min = bottom_y - 1.5
    y_max = (n_samples - 1) * track_spacing + track_spacing
    ax.set_ylim(y_min, y_max)
    ax.axis('off')
    
    plt.tight_layout()
    
    for ext in ['png', 'pdf']:
        outpath = OUT_DIR / f"{output_name}.{ext}"
        fig.savefig(str(outpath), dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
    plt.close()
    print(f"  Saved: {output_name}")


# ============================================================
# MAIN
# ============================================================
print("=" * 60)
print("High-Quality EasyFig v2 - Publication Ready")
print("=" * 60)

# 1. Full KPC plasmid map
print("\n1. Full KPC plasmid map (contig_3, ~131kb)...")
create_full_plasmid_map(
    3, SAMPLES,
    'Genomic Organization of blaKPC-2-carrying IncFII+IncR Plasmid (~131 kb)\n'
    'Linear Comparison Across 8 ST11 K. pneumoniae Isolates',
    'HQ_EasyFig_KPC_Full'
)

# 2. Full NDM plasmid map
print("\n2. Full NDM plasmid map (contig_2, ~354kb)...")
create_full_plasmid_map(
    2, NDM_SAMPLES,
    'Genomic Organization of blaNDM-1-carrying IncFIB(pNDM-Mar)+IncHI1B Plasmid (~354 kb)\n'
    'Linear Comparison Across 6 ST11 K. pneumoniae Isolates',
    'HQ_EasyFig_NDM_Full'
)

# 3. Zoomed KPC region (~25-70kb covering the AMR cluster)
print("\n3. Zoomed KPC resistance region (25-75kb)...")
create_zoomed_region(
    3, SAMPLES, 25000, 75000,
    'Zoomed View: AMR Gene Cluster on blaKPC-2 Plasmid (25-75 kb)\n'
    'blaKPC-2, blaSHV-12, blaTEM-1, rmtB1, catA2, fosA3, Mercury Resistance',
    'HQ_EasyFig_KPC_Zoomed'
)

# 4. Zoomed NDM region (~200-240kb covering the NDM cluster)
print("\n4. Zoomed NDM resistance region (200-250kb)...")
create_zoomed_region(
    2, NDM_SAMPLES, 200000, 250000,
    'Zoomed View: AMR Gene Cluster on blaNDM-1 Plasmid (200-250 kb)\n'
    'blaNDM-1, armA, ble, sul1, dfrA5, aph(3\')-Ia, mph(E), msr(E)',
    'HQ_EasyFig_NDM_Zoomed'
)

print("\n" + "=" * 60)
print("All figures complete!")
print("=" * 60)
