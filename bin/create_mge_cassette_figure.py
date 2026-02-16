#!/usr/bin/env python3
"""
Create a high-quality multi-panel figure showing ALL AMR gene cassettes
with complete MGE context (IS elements, transposons, integrases, integrons)
on both KPC and NDM plasmids from KP21915 reference genome.

Publication-quality for Lancet-level journals.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
from pathlib import Path

OUT_DIR = Path("/home/ubuntu/plasmid_analysis/08_BRIG_EasyFig")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Parse all genes from KP21915 TSV
def parse_genes(tsv_path, target_contig):
    genes = []
    with open(tsv_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            contig = parts[0]
            if contig != target_contig:
                continue
            feat_type = parts[1]
            start = int(parts[2])
            end = int(parts[3])
            strand = parts[4]
            gene = parts[6] if parts[6] else ""
            product = parts[7] if len(parts) > 7 else ""
            dbxref = parts[8] if len(parts) > 8 else ""
            genes.append({
                'start': start, 'end': end, 'strand': strand,
                'gene': gene, 'product': product, 'dbxref': dbxref
            })
    return sorted(genes, key=lambda x: x['start'])

TSV = "/home/ubuntu/plasmid_analysis/04_AMR_Plasmid_Association/Bakta_Annotations/KP21915/KP21915.tsv"
kpc_genes = parse_genes(TSV, "contig_3")
ndm_genes = parse_genes(TSV, "contig_2")

# Classification functions
def classify_gene(g):
    """Classify a gene into categories for coloring."""
    gene = g['gene'].lower()
    product = g['product'].lower()
    
    # AMR - Carbapenemase
    if any(k in gene for k in ['blakpc', 'blandm', 'blaoxa', 'blavim', 'blaimp']):
        return 'carbapenemase'
    # AMR - Beta-lactamase
    if any(k in gene for k in ['blactx', 'blashv', 'blatem']):
        return 'beta_lactamase'
    # AMR - Aminoglycoside
    if any(k in gene for k in ['aph(', 'aac(', 'ant(', 'rmtb', 'arma']) or 'aminoglycoside' in product:
        return 'aminoglycoside'
    # AMR - Quinolone
    if any(k in gene for k in ['qnr', 'oqx']):
        return 'quinolone'
    # AMR - Macrolide
    if any(k in gene for k in ['mph(', 'msr(', 'erm']):
        return 'macrolide'
    # AMR - Sulfonamide
    if 'sul' in gene and any(c.isdigit() for c in gene):
        return 'sulfonamide'
    # AMR - Trimethoprim
    if 'dfra' in gene or 'dfr' in gene:
        return 'trimethoprim'
    # AMR - Phenicol
    if any(k in gene for k in ['cata', 'catb', 'flor', 'cmla']):
        return 'phenicol'
    # AMR - Fosfomycin
    if 'fosa' in gene:
        return 'fosfomycin'
    # AMR - Bleomycin
    if gene == 'ble' or 'bleomycin' in product:
        return 'bleomycin'
    # AMR - QacE
    if 'qac' in gene:
        return 'qac_efflux'
    # IS elements / Transposases
    if 'tnp' in gene or 'transpos' in product.lower():
        return 'is_element'
    # Tn3 transposase
    if 'tn3' in product.lower():
        return 'is_element'
    # Integrase / Integron
    if 'inti' in gene or 'integrase' in product.lower() or 'integron' in product.lower():
        return 'integrase'
    # Recombinase / Resolvase
    if 'recombinas' in product.lower() or 'resolvas' in product.lower():
        return 'recombinase'
    # Mercury resistance
    if gene.startswith('mer') or 'mercury' in product.lower():
        return 'mercury'
    # Replication
    if 'replicat' in product.lower() or gene.startswith('rep'):
        return 'replication'
    # Conjugation / Transfer
    if 'conjug' in product.lower() or 'transfer' in product.lower() or gene.startswith('tra') or gene.startswith('trb'):
        return 'conjugation'
    # Toxin-antitoxin
    if 'toxin' in product.lower() or 'antitoxin' in product.lower():
        return 'ta_system'
    # Hypothetical
    if 'hypothetical' in product.lower():
        return 'hypothetical'
    # Other
    return 'other'

# Color scheme - Lancet-quality
COLORS = {
    'carbapenemase': '#D32F2F',      # Deep red
    'beta_lactamase': '#F44336',     # Red
    'aminoglycoside': '#FF9800',     # Orange
    'quinolone': '#FFC107',          # Amber
    'macrolide': '#9C27B0',          # Purple
    'sulfonamide': '#3F51B5',        # Indigo
    'trimethoprim': '#2196F3',       # Blue
    'phenicol': '#00BCD4',           # Cyan
    'fosfomycin': '#E91E63',         # Pink
    'bleomycin': '#795548',          # Brown
    'qac_efflux': '#607D8B',        # Blue grey
    'is_element': '#FFD54F',         # Light amber/yellow
    'integrase': '#FF6F00',          # Dark amber
    'recombinase': '#F9A825',        # Yellow-orange
    'mercury': '#78909C',            # Grey-blue
    'replication': '#81C784',        # Light green
    'conjugation': '#A5D6A7',        # Lighter green
    'ta_system': '#CE93D8',          # Light purple
    'hypothetical': '#E0E0E0',       # Light grey
    'other': '#BDBDBD',             # Grey
}

def get_label(g):
    """Get a concise label for a gene."""
    gene = g['gene']
    product = g['product']
    cat = classify_gene(g)
    
    if gene and gene != 'tnp':
        return gene
    
    # For IS elements, extract IS family name
    if cat == 'is_element':
        for token in ['IS1B', 'IS15DIV', 'IS15', 'IS903B', 'IS5075', 'ISKpn27', 'ISKpn26', 
                       'ISKpn14', 'ISCfr3', 'ISKpn19', 'ISKpn21', 'ISEc29', 'IS4321',
                       'ISVsa3', 'ISKra4']:
            if token in product:
                return token
        # Try to extract IS name from product
        if 'IS' in product:
            parts = product.split()
            for p in parts:
                if p.startswith('IS') and len(p) > 2:
                    return p
        if 'Tn3' in product:
            return 'Tn3'
        if 'Mutator' in product:
            return 'MuA'
        return 'tnp'
    
    if cat == 'integrase':
        if 'intI1' in gene:
            return 'intI1'
        return 'int'
    
    if cat == 'recombinase':
        return 'res'
    
    if cat == 'mercury':
        return gene if gene else 'mer'
    
    if cat == 'hypothetical':
        return ''
    
    # Shorten long product names
    if len(product) > 20:
        words = product.split()
        return words[0][:12] if words else ''
    return product[:12] if product else ''

def draw_gene_arrow(ax, start, end, strand, y, height, color, alpha=1.0, linewidth=0.5):
    """Draw a gene as a directional arrow."""
    gene_len = end - start
    arrow_head = min(gene_len * 0.25, 400)  # Arrow head size
    
    if strand == '+':
        # Forward arrow
        body_end = end - arrow_head
        xs = [start, body_end, body_end, end, body_end, body_end, start]
        ys = [y - height/2, y - height/2, y - height*0.7/2, y, 
              y + height*0.7/2, y + height/2, y + height/2]
    else:
        # Reverse arrow
        body_start = start + arrow_head
        xs = [end, body_start, body_start, start, body_start, body_start, end]
        ys = [y - height/2, y - height/2, y - height*0.7/2, y, 
              y + height*0.7/2, y + height/2, y + height/2]
    
    ax.fill(xs, ys, fc=color, ec='#333333', linewidth=linewidth, alpha=alpha, zorder=3)

def draw_cassette_panel(ax, genes, region_start, region_end, title, y_center=0, 
                        gene_height=0.6, show_scale=True, panel_label=None):
    """Draw one cassette panel with all genes annotated."""
    region_len = region_end - region_start
    
    # Filter genes in region
    region_genes = [g for g in genes if g['end'] >= region_start and g['start'] <= region_end]
    
    # Draw backbone line
    ax.plot([region_start, region_end], [y_center, y_center], 
            color='#424242', linewidth=2, zorder=1)
    
    # Draw genes
    label_positions = []
    for g in region_genes:
        cat = classify_gene(g)
        color = COLORS.get(cat, '#BDBDBD')
        s = max(g['start'], region_start)
        e = min(g['end'], region_end)
        
        # Skip very tiny hypothetical proteins
        if cat == 'hypothetical' and (e - s) < region_len * 0.005:
            continue
        
        alpha = 0.4 if cat == 'hypothetical' else 0.9
        lw = 0.3 if cat == 'hypothetical' else 0.8
        
        draw_gene_arrow(ax, s, e, g['strand'], y_center, gene_height, color, alpha=alpha, linewidth=lw)
        
        # Add label
        label = get_label(g)
        if not label:
            continue
        
        mid = (s + e) / 2
        
        # Determine label style based on category
        is_amr = cat in ['carbapenemase', 'beta_lactamase', 'aminoglycoside', 'quinolone', 
                         'macrolide', 'sulfonamide', 'trimethoprim', 'phenicol', 'fosfomycin', 'bleomycin']
        is_mge = cat in ['is_element', 'integrase', 'recombinase']
        
        if is_amr:
            fontsize = 7
            fontweight = 'bold'
            label_y = y_center + gene_height * 1.3
            color_text = COLORS[cat]
            fontstyle = 'italic'
            rotation = 35
        elif is_mge:
            fontsize = 5.5
            fontweight = 'bold'
            label_y = y_center - gene_height * 1.3
            color_text = '#BF360C' if cat == 'integrase' else '#F57F17'
            fontstyle = 'normal'
            rotation = 35
        elif cat == 'mercury':
            fontsize = 5.5
            fontweight = 'normal'
            label_y = y_center + gene_height * 1.1
            color_text = '#546E7A'
            fontstyle = 'italic'
            rotation = 35
        elif cat == 'hypothetical':
            continue
        else:
            fontsize = 5
            fontweight = 'normal'
            label_y = y_center + gene_height * 1.1
            color_text = '#616161'
            fontstyle = 'italic'
            rotation = 35
        
        # Check for label overlap
        too_close = False
        for lp in label_positions:
            if abs(mid - lp) < region_len * 0.025:
                too_close = True
                # Offset the label
                label_y += gene_height * 0.6 if label_y > y_center else -gene_height * 0.6
                break
        
        ax.text(mid, label_y, label, fontsize=fontsize, fontweight=fontweight,
                fontstyle=fontstyle, color=color_text, ha='center', va='bottom' if label_y > y_center else 'top',
                rotation=rotation, rotation_mode='anchor', zorder=5)
        label_positions.append(mid)
    
    # Scale bar
    if show_scale:
        scale_len = 1000  # 1 kb
        if region_len > 20000:
            scale_len = 5000  # 5 kb
        scale_x = region_start + region_len * 0.02
        scale_y = y_center - gene_height * 2.5
        ax.plot([scale_x, scale_x + scale_len], [scale_y, scale_y], 
                color='black', linewidth=1.5, zorder=5)
        ax.text(scale_x + scale_len/2, scale_y - 0.15, f'{scale_len//1000} kb',
                fontsize=6, ha='center', va='top', fontweight='bold')
    
    # Panel label (A, B, C, etc.)
    if panel_label:
        ax.text(region_start - region_len * 0.02, y_center + gene_height * 2.5,
                panel_label, fontsize=14, fontweight='bold', ha='right', va='center')
    
    # Title
    ax.text((region_start + region_end) / 2, y_center + gene_height * 3.2,
            title, fontsize=9, fontweight='bold', ha='center', va='bottom',
            style='italic')
    
    ax.set_xlim(region_start - region_len * 0.05, region_end + region_len * 0.05)
    ax.set_ylim(y_center - gene_height * 4, y_center + gene_height * 5)
    ax.axis('off')


# ============================================================
# CREATE THE MULTI-PANEL FIGURE
# ============================================================

# Define cassette regions
cassettes = [
    # KPC plasmid cassettes
    {
        'genes': kpc_genes,
        'start': 2500, 'end': 13000,
        'title': 'KPC Plasmid — Cassette 1: blaCTX-M-65 + fosA3 region',
        'label': 'A',
        'plasmid': 'pKPC (~131 kb, IncFII+IncR)'
    },
    {
        'genes': kpc_genes,
        'start': 31000, 'end': 52000,
        'title': 'KPC Plasmid — Cassette 2: blaSHV-12 + blaKPC-2 + Integrase region',
        'label': 'B',
        'plasmid': 'pKPC (~131 kb, IncFII+IncR)'
    },
    {
        'genes': kpc_genes,
        'start': 51000, 'end': 67000,
        'title': 'KPC Plasmid — Cassette 3: Mercury resistance + blaTEM-1 + rmtB1 region',
        'label': 'C',
        'plasmid': 'pKPC (~131 kb, IncFII+IncR)'
    },
    {
        'genes': kpc_genes,
        'start': 71000, 'end': 82000,
        'title': 'KPC Plasmid — Cassette 4: catA2 region',
        'label': 'D',
        'plasmid': 'pKPC (~131 kb, IncFII+IncR)'
    },
    # NDM plasmid cassettes
    {
        'genes': ndm_genes,
        'start': 196000, 'end': 220000,
        'title': 'NDM Plasmid — Cassette 5: Class 1 integron + blaNDM-1 + armA + msr(E)/mph(E) region',
        'label': 'E',
        'plasmid': 'pNDM (~354 kb, IncFIB+IncHI1B)'
    },
    {
        'genes': ndm_genes,
        'start': 220000, 'end': 243000,
        'title': 'NDM Plasmid — Cassette 6: qnrS1 + aph(3\')-VI + sul2 region',
        'label': 'F',
        'plasmid': 'pNDM (~354 kb, IncFIB+IncHI1B)'
    },
]

# Create figure
fig, axes = plt.subplots(len(cassettes), 1, figsize=(20, len(cassettes) * 3.2))
fig.suptitle('Genetic Environment of AMR Gene Cassettes in ST11 K. pneumoniae (KP21915)\n'
             'Mobile Genetic Elements Surrounding Key Resistance Determinants',
             fontsize=14, fontweight='bold', y=0.98)

for i, cas in enumerate(cassettes):
    ax = axes[i]
    draw_cassette_panel(ax, cas['genes'], cas['start'], cas['end'],
                       cas['title'], panel_label=cas['label'])

# Add legend at bottom
legend_elements = [
    mpatches.Patch(facecolor=COLORS['carbapenemase'], edgecolor='#333', label='Carbapenemase', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['beta_lactamase'], edgecolor='#333', label='Beta-lactamase (ESBL)', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['aminoglycoside'], edgecolor='#333', label='Aminoglycoside resistance', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['quinolone'], edgecolor='#333', label='Quinolone resistance', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['macrolide'], edgecolor='#333', label='Macrolide resistance', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['sulfonamide'], edgecolor='#333', label='Sulfonamide resistance', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['trimethoprim'], edgecolor='#333', label='Trimethoprim resistance', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['phenicol'], edgecolor='#333', label='Phenicol resistance', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['fosfomycin'], edgecolor='#333', label='Fosfomycin resistance', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['bleomycin'], edgecolor='#333', label='Bleomycin resistance', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['is_element'], edgecolor='#333', label='IS element / Transposase', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['integrase'], edgecolor='#333', label='Integrase / Integron', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['recombinase'], edgecolor='#333', label='Recombinase / Resolvase', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['mercury'], edgecolor='#333', label='Mercury resistance', linewidth=0.5),
    mpatches.Patch(facecolor=COLORS['other'], edgecolor='#333', label='Other / Hypothetical', linewidth=0.5),
]

fig.legend(handles=legend_elements, loc='lower center', ncol=5, fontsize=7.5,
           frameon=True, fancybox=True, shadow=False, edgecolor='#BDBDBD',
           bbox_to_anchor=(0.5, 0.005))

plt.tight_layout(rect=[0.02, 0.06, 0.98, 0.96])

# Save
for fmt in ['png', 'pdf']:
    outpath = OUT_DIR / f"HQ_MGE_Cassettes_AllRegions.{fmt}"
    fig.savefig(str(outpath), dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {outpath}")

plt.close()
print("\nDone! Multi-panel MGE cassette figure created.")
