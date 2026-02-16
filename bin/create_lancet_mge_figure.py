#!/usr/bin/env python3
"""
Publication-grade multi-panel figure: Genetic environment of AMR gene cassettes
in ST11 K. pneumoniae (KP21915).

Color scheme:
  - RED: All antimicrobial resistance genes
  - BLUE: All IS elements / transposases
  - ORANGE: Integrases / integrons
  - GREY: Mercury resistance / other accessory genes
  - LIGHT GREY: Hypothetical proteins / other CDS

Transposon cassettes are delineated with bracket lines and shaded backgrounds.

Tools used (for citation):
  - Bakta v1.9.4 (Schwengers et al., 2021, Microb Genom)
  - AMRFinderPlus v3.12.8 (Feldgarden et al., 2021, Sci Rep)
  - Python matplotlib v3.x (Hunter, 2007, Comput Sci Eng)
  - pyGenomeViz v1.6.1 (Shimoyama, 2022) — for reference
  - ISfinder database (Siguier et al., 2006, Nucleic Acids Res)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================
OUT_DIR = Path("/home/ubuntu/plasmid_analysis/08_BRIG_EasyFig")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Apply a clean style first, then set fonts
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 8,
    'axes.linewidth': 0.5,
    'grid.linewidth': 0,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'pdf.fonttype': 42,  # TrueType fonts in PDF (required by journals)
    'ps.fonttype': 42,
})

# Color palette — clean, high-contrast, journal-ready
C_AMR = '#C62828'          # Deep red for resistance genes
C_AMR_EDGE = '#8E0000'     # Darker red edge
C_IS = '#1565C0'           # Deep blue for IS elements
C_IS_EDGE = '#0D47A1'      # Darker blue edge
C_INTEGRASE = '#E65100'    # Deep orange for integrases
C_INT_EDGE = '#BF360C'     # Darker orange edge
C_MERCURY = '#546E7A'      # Blue-grey for mercury resistance
C_MER_EDGE = '#37474F'
C_OTHER = '#BDBDBD'        # Light grey for other genes
C_OTHER_EDGE = '#9E9E9E'
C_HYPO = '#E0E0E0'         # Very light grey for hypothetical
C_HYPO_EDGE = '#BDBDBD'
C_BACKBONE = '#212121'     # Near-black backbone
C_TN_BG = '#FFF3E0'        # Very light orange for transposon cassette background
C_TN_BORDER = '#FF8F00'    # Amber for transposon cassette border

# ============================================================
# PARSE GENES
# ============================================================
def parse_genes(tsv_path, target_contig):
    genes = []
    with open(tsv_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            if parts[0] != target_contig:
                continue
            genes.append({
                'start': int(parts[2]), 'end': int(parts[3]),
                'strand': parts[4], 'gene': parts[6],
                'product': parts[7],
                'dbxref': parts[8] if len(parts) > 8 else ""
            })
    return sorted(genes, key=lambda x: x['start'])

TSV = "/home/ubuntu/plasmid_analysis/04_AMR_Plasmid_Association/Bakta_Annotations/KP21915/KP21915.tsv"
kpc_genes = parse_genes(TSV, "contig_3")
ndm_genes = parse_genes(TSV, "contig_2")

# ============================================================
# GENE CLASSIFICATION
# ============================================================
def classify(g):
    gene = g['gene'].lower()
    product = g['product'].lower()
    
    # AMR genes
    amr_prefixes = ['bla', 'aph(', 'aac(', 'ant(', 'rmtb', 'arma', 'fosa', 'cata', 'catb',
                     'sul1', 'sul2', 'dfra', 'dfrb', 'mph(', 'msr(', 'erm', 'qnr', 'oqx',
                     'tet(', 'mcr', 'flor', 'cmla', 'ble', 'qac']
    for prefix in amr_prefixes:
        if gene.startswith(prefix) or prefix in gene:
            return 'amr'
    if any(k in product for k in ['beta-lactamase', 'aminoglycoside', 'chloramphenicol',
                                   'fosfomycin', 'sulfonamide', 'dihydrofolate', 'bleomycin',
                                   'macrolide', 'ribosomal protection', 'phosphotransferase']):
        if 'transpos' not in product:
            return 'amr'
    
    # IS elements / Transposases
    if gene == 'tnp' or 'transpos' in product or 'tn3' in product.lower():
        return 'is_element'
    if any(k in product.lower() for k in ['is1 ', 'is5 ', 'is6 ', 'is26', 'is110', 'is481',
                                            'is903', 'is4321', 'isncy', 'is91', 'iskpn',
                                            'isec', 'iscfr', 'iskra', 'is5075', 'is15',
                                            'mutator family']):
        return 'is_element'
    
    # Integrases / Integrons
    if 'inti' in gene or 'integrase' in product or 'integron' in product:
        return 'integrase'
    
    # Recombinase / Resolvase
    if 'recombinas' in product or 'resolvas' in product or 'invertase' in product:
        return 'recombinase'
    
    # Mercury resistance
    if gene.startswith('mer') or 'mercury' in product:
        return 'mercury'
    
    # Hypothetical
    if 'hypothetical' in product:
        return 'hypothetical'
    
    return 'other'

def get_color(cat):
    return {
        'amr': (C_AMR, C_AMR_EDGE),
        'is_element': (C_IS, C_IS_EDGE),
        'integrase': (C_INTEGRASE, C_INT_EDGE),
        'recombinase': (C_INTEGRASE, C_INT_EDGE),
        'mercury': (C_MERCURY, C_MER_EDGE),
        'hypothetical': (C_HYPO, C_HYPO_EDGE),
        'other': (C_OTHER, C_OTHER_EDGE),
    }.get(cat, (C_OTHER, C_OTHER_EDGE))

def get_is_name(g):
    """Extract the IS element family name from product description."""
    product = g['product']
    is_names = ['IS1B', 'IS15DIV', 'IS15', 'IS903B', 'IS5075', 'ISKpn27', 'ISKpn26',
                'ISKpn14', 'ISKpn19', 'ISKpn21', 'ISEc29', 'IS4321', 'ISVsa3',
                'ISKra4', 'ISCfr3', 'ISKpn8']
    for name in is_names:
        if name in product:
            return name
    if 'Tn3' in product:
        return 'Tn3'
    if 'Mutator' in product:
        return 'MuA'
    if 'IS91' in product:
        return 'IS91'
    return 'IS'

def get_label(g):
    """Get concise label for a gene."""
    cat = classify(g)
    gene = g['gene']
    
    if cat == 'amr':
        if gene:
            return gene
        # Try to extract from product
        prod = g['product']
        if 'beta-lactamase' in prod.lower():
            return 'bla'
        return prod[:15]
    elif cat == 'is_element':
        return get_is_name(g)
    elif cat == 'integrase':
        if 'intI1' in gene:
            return 'intI1'
        return 'int'
    elif cat == 'recombinase':
        if 'resolv' in g['product'].lower():
            return 'res'
        if 'invertase' in g['product'].lower():
            return 'hin'
        return 'rec'
    elif cat == 'mercury':
        return gene if gene else 'mer'
    elif cat == 'hypothetical':
        return ''
    else:
        if gene:
            return gene
        return ''

# ============================================================
# DRAWING FUNCTIONS
# ============================================================
def draw_gene_arrow(ax, start, end, strand, y, height, facecolor, edgecolor, alpha=1.0, lw=0.8):
    """Draw a gene as a clean directional arrow."""
    gene_len = end - start
    arrow_head = min(gene_len * 0.2, 500)
    
    if strand == '+':
        body_end = end - arrow_head
        xs = [start, body_end, body_end, end, body_end, body_end, start, start]
        ys = [y - height/2, y - height/2, y - height*0.75/2, y,
              y + height*0.75/2, y + height/2, y + height/2, y - height/2]
    else:
        body_start = start + arrow_head
        xs = [end, body_start, body_start, start, body_start, body_start, end, end]
        ys = [y - height/2, y - height/2, y - height*0.75/2, y,
              y + height*0.75/2, y + height/2, y + height/2, y - height/2]
    
    ax.fill(xs, ys, fc=facecolor, ec=edgecolor, linewidth=lw, alpha=alpha, zorder=4)

def draw_transposon_bracket(ax, start, end, y, height, label="Tn", color=C_TN_BORDER):
    """Draw a bracket above/below genes to indicate a transposon cassette."""
    # Shaded background
    ax.fill_between([start, end], y - height * 2.2, y + height * 2.2,
                    alpha=0.08, color=color, zorder=1)
    
    # Top bracket line
    bracket_y = y + height * 2.0
    tick = height * 0.3
    ax.plot([start, start, end, end], 
            [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
            color=color, linewidth=1.2, zorder=5, solid_capstyle='round')
    
    # Label above bracket
    ax.text((start + end) / 2, bracket_y + tick * 1.5, label,
            fontsize=7, fontweight='bold', fontstyle='italic', color=color,
            ha='center', va='bottom', zorder=6)

def draw_cassette(ax, genes, region_start, region_end, title, panel_label,
                  transposon_regions=None, y=0, gene_height=0.5):
    """Draw one complete cassette panel."""
    region_len = region_end - region_start
    region_genes = [g for g in genes if g['end'] >= region_start and g['start'] <= region_end]
    
    # Draw transposon cassette backgrounds first
    if transposon_regions:
        for tn in transposon_regions:
            draw_transposon_bracket(ax, tn['start'], tn['end'], y, gene_height, tn['label'])
    
    # Backbone line
    ax.plot([region_start, region_end], [y, y],
            color=C_BACKBONE, linewidth=1.5, zorder=2, solid_capstyle='round')
    
    # Position tick marks every 5 kb
    tick_interval = 5000
    first_tick = ((region_start // tick_interval) + 1) * tick_interval
    for pos in range(first_tick, region_end, tick_interval):
        ax.plot([pos, pos], [y - gene_height * 0.3, y + gene_height * 0.3],
                color='#757575', linewidth=0.4, zorder=2)
        ax.text(pos, y - gene_height * 1.5, f'{pos//1000}',
                fontsize=4.5, color='#757575', ha='center', va='top')
    
    # Draw genes
    amr_labels = []
    is_labels = []
    other_labels = []
    
    for g in region_genes:
        cat = classify(g)
        fc, ec = get_color(cat)
        s = max(g['start'], region_start)
        e = min(g['end'], region_end)
        
        if cat == 'hypothetical' and (e - s) < region_len * 0.004:
            continue
        
        alpha = 0.35 if cat == 'hypothetical' else 0.92
        lw = 0.3 if cat == 'hypothetical' else 0.8
        
        draw_gene_arrow(ax, s, e, g['strand'], y, gene_height, fc, ec, alpha=alpha, lw=lw)
        
        label = get_label(g)
        if not label:
            continue
        
        mid = (s + e) / 2
        
        if cat == 'amr':
            amr_labels.append((mid, label, s, e))
        elif cat in ['is_element', 'integrase', 'recombinase']:
            is_labels.append((mid, label, s, e))
        elif cat == 'mercury':
            other_labels.append((mid, label, s, e))
        elif cat not in ['hypothetical', 'other']:
            other_labels.append((mid, label, s, e))
    
    # Place AMR labels ABOVE (red, bold, italic)
    _place_labels(ax, amr_labels, y, gene_height, region_len,
                  fontsize=7.5, fontweight='bold', fontstyle='italic',
                  color=C_AMR, above=True, offset_mult=1.2)
    
    # Place IS/integrase labels BELOW (blue, bold)
    _place_labels(ax, is_labels, y, gene_height, region_len,
                  fontsize=6, fontweight='bold', fontstyle='normal',
                  color=C_IS, above=False, offset_mult=1.2)
    
    # Place mercury/other labels above (grey, italic)
    _place_labels(ax, other_labels, y, gene_height, region_len,
                  fontsize=5.5, fontweight='normal', fontstyle='italic',
                  color='#546E7A', above=True, offset_mult=2.5)
    
    # Scale bar
    scale_len = 1000 if region_len < 15000 else 5000
    scale_x = region_start + region_len * 0.02
    scale_y = y - gene_height * 4.5
    ax.plot([scale_x, scale_x + scale_len], [scale_y, scale_y],
            color='black', linewidth=1.5, zorder=5)
    ax.plot([scale_x, scale_x], [scale_y - 0.05, scale_y + 0.05],
            color='black', linewidth=1.0, zorder=5)
    ax.plot([scale_x + scale_len, scale_x + scale_len], [scale_y - 0.05, scale_y + 0.05],
            color='black', linewidth=1.0, zorder=5)
    ax.text(scale_x + scale_len / 2, scale_y - 0.12, f'{scale_len // 1000} kb',
            fontsize=6, ha='center', va='top', fontweight='bold')
    
    # Panel label (A, B, C...)
    ax.text(region_start - region_len * 0.03, y + gene_height * 3.5,
            panel_label, fontsize=16, fontweight='bold', ha='right', va='center',
            fontfamily='serif')
    
    # Title
    ax.set_title(title, fontsize=9, fontweight='bold', fontstyle='italic',
                 pad=8, loc='center')
    
    # Position label (kb)
    ax.text(region_end + region_len * 0.01, y - gene_height * 1.5, 'kb',
            fontsize=5, color='#757575', ha='left', va='top')
    
    ax.set_xlim(region_start - region_len * 0.06, region_end + region_len * 0.06)
    ax.set_ylim(y - gene_height * 6, y + gene_height * 6)
    ax.axis('off')

def _place_labels(ax, labels, y, gene_height, region_len,
                  fontsize, fontweight, fontstyle, color, above, offset_mult):
    """Place labels with anti-overlap logic."""
    if not labels:
        return
    
    min_gap = region_len * 0.035
    placed = []
    
    # Sort by position
    labels_sorted = sorted(labels, key=lambda x: x[0])
    
    for mid, label, s, e in labels_sorted:
        # Check overlap with placed labels
        tier = 0
        for pmid in placed:
            if abs(mid - pmid[0]) < min_gap and pmid[1] == tier:
                tier += 1
        
        if above:
            label_y = y + gene_height * (offset_mult + tier * 1.0)
            va = 'bottom'
        else:
            label_y = y - gene_height * (offset_mult + tier * 1.0)
            va = 'top'
        
        # Draw a thin connecting line from gene to label
        gene_edge_y = y + gene_height * 0.5 if above else y - gene_height * 0.5
        ax.plot([mid, mid], [gene_edge_y, label_y], color=color, linewidth=0.3,
                alpha=0.5, zorder=3, linestyle=':')
        
        ax.text(mid, label_y, label, fontsize=fontsize, fontweight=fontweight,
                fontstyle=fontstyle, color=color, ha='center', va=va,
                rotation=30, rotation_mode='anchor', zorder=6)
        
        placed.append((mid, tier))

# ============================================================
# DEFINE CASSETTES AND TRANSPOSON REGIONS
# ============================================================
cassettes = [
    {
        'genes': kpc_genes,
        'start': 2500, 'end': 13000,
        'title': 'pKPC (~131 kb, IncFII+IncR) — blaCTX-M-65 / fosA3 resistance region',
        'label': 'A',
        'transposons': [
            {'start': 3300, 'end': 7600, 'label': 'Tn (IS1B–IS15DIV–blaCTX-M-65–IS903B)'},
            {'start': 9200, 'end': 12100, 'label': 'Tn (IS15DIV–fosA3–IS15DIV)'},
        ]
    },
    {
        'genes': kpc_genes,
        'start': 31000, 'end': 52000,
        'title': 'pKPC (~131 kb, IncFII+IncR) — blaSHV-12 / blaKPC-2 / Integrase resistance region',
        'label': 'B',
        'transposons': [
            {'start': 32300, 'end': 37600, 'label': 'Tn (IS15DIV–blaSHV-12)'},
            {'start': 39700, 'end': 43300, 'label': 'Tn3-ISKpn27–blaKPC-2'},
            {'start': 49000, 'end': 51100, 'label': 'Tn (IS15DIV–Integrase)'},
        ]
    },
    {
        'genes': kpc_genes,
        'start': 51000, 'end': 67000,
        'title': 'pKPC (~131 kb, IncFII+IncR) — Mercury resistance / blaTEM-1 / rmtB1 region',
        'label': 'C',
        'transposons': [
            {'start': 52000, 'end': 57000, 'label': 'Tn (mer operon–IS5075)'},
            {'start': 58900, 'end': 65200, 'label': 'Tn (IS15DIV–blaTEM-1–rmtB1–IS15)'},
        ]
    },
    {
        'genes': kpc_genes,
        'start': 71000, 'end': 82000,
        'title': 'pKPC (~131 kb, IncFII+IncR) — catA2 resistance region',
        'label': 'D',
        'transposons': [
            {'start': 72000, 'end': 81200, 'label': 'Tn (IS15DIV–ISCfr3–ISKpn26–catA2–IS15DIV)'},
        ]
    },
    {
        'genes': ndm_genes,
        'start': 196000, 'end': 220000,
        'title': 'pNDM (~354 kb, IncFIB+IncHI1B) — Class 1 integron / blaNDM-1 / armA / msr(E) / mph(E) region',
        'label': 'E',
        'transposons': [
            {'start': 198000, 'end': 205600, 'label': 'Class 1 integron (IS4321–sul1–dfrA5–intI1)'},
            {'start': 205700, 'end': 213000, 'label': 'Tn (IS15DIV–aph(3\')-Ia–ble–blaNDM-1–IS15DIV)'},
            {'start': 213200, 'end': 219000, 'label': 'Tn (armA–ISEc29–msr(E)–mph(E)–ISKpn21)'},
        ]
    },
    {
        'genes': ndm_genes,
        'start': 220000, 'end': 243000,
        'title': 'pNDM (~354 kb, IncFIB+IncHI1B) — qnrS1 / aph(3\')-VI / sul2 region',
        'label': 'F',
        'transposons': [
            {'start': 222700, 'end': 230100, 'label': 'Tn (IS15DIV–qnrS1–ISKpn19)'},
            {'start': 231700, 'end': 242100, 'label': 'Tn (aph(3\')-VI–int–MuA–IS15DIV–sul2–ISVsa3–res–IS15DIV)'},
        ]
    },
]

# ============================================================
# CREATE FIGURE
# ============================================================
n_panels = len(cassettes)
fig, axes = plt.subplots(n_panels, 1, figsize=(22, n_panels * 3.0 + 1.5))

# Main title
fig.suptitle(
    'Genetic Environment of Antimicrobial Resistance Gene Cassettes\n'
    'in ST11 Klebsiella pneumoniae KP21915',
    fontsize=13, fontweight='bold', y=0.995, fontfamily='serif'
)

for i, cas in enumerate(cassettes):
    draw_cassette(axes[i], cas['genes'], cas['start'], cas['end'],
                  cas['title'], cas['label'],
                  transposon_regions=cas.get('transposons'))

# Legend
legend_elements = [
    mpatches.Patch(facecolor=C_AMR, edgecolor=C_AMR_EDGE, linewidth=0.8,
                   label='Antimicrobial resistance gene'),
    mpatches.Patch(facecolor=C_IS, edgecolor=C_IS_EDGE, linewidth=0.8,
                   label='IS element / Transposase'),
    mpatches.Patch(facecolor=C_INTEGRASE, edgecolor=C_INT_EDGE, linewidth=0.8,
                   label='Integrase / Integron'),
    mpatches.Patch(facecolor=C_MERCURY, edgecolor=C_MER_EDGE, linewidth=0.8,
                   label='Mercury resistance'),
    mpatches.Patch(facecolor=C_OTHER, edgecolor=C_OTHER_EDGE, linewidth=0.8,
                   label='Other CDS'),
    mpatches.Patch(facecolor=C_HYPO, edgecolor=C_HYPO_EDGE, linewidth=0.8,
                   label='Hypothetical protein'),
    mpatches.Patch(facecolor=C_TN_BG, edgecolor=C_TN_BORDER, linewidth=1.2,
                   label='Transposon cassette', linestyle='-'),
]

fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=8,
           frameon=True, fancybox=False, edgecolor='#9E9E9E',
           bbox_to_anchor=(0.5, 0.002), handlelength=1.5, handleheight=1.0,
           columnspacing=1.5)

plt.tight_layout(rect=[0.03, 0.05, 0.98, 0.97], h_pad=2.5)

# Save in both formats
for fmt in ['png', 'pdf', 'svg']:
    outpath = OUT_DIR / f"Figure_MGE_Cassettes_Lancet.{fmt}"
    fig.savefig(str(outpath), dpi=300, bbox_inches='tight', facecolor='white',
                transparent=False)
    print(f"Saved: {outpath}")

plt.close()
print("\nDone! Lancet-quality MGE cassette figure created.")
