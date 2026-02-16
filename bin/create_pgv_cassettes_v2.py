#!/usr/bin/env python3
"""
A+++++ Publication-grade multi-panel figure: Genetic environment of AMR gene
cassettes in ST11 Klebsiella pneumoniae KP21915.

Tool: pyGenomeViz v1.6.1 (Shimoyama, 2022)
Citation: Shimoyama Y. pyGenomeViz: A genome visualization python package
          for comparative genomics. https://github.com/moshi4/pyGenomeViz
"""
import matplotlib
matplotlib.use('Agg')

from pygenomeviz import GenomeViz
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

OUT_DIR = Path("/home/ubuntu/plasmid_analysis/08_BRIG_EasyFig")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# COLORS — High contrast, journal-ready
# ============================================================
C_AMR       = '#C62828'   # Deep red
C_IS        = '#1565C0'   # Deep blue
C_INTEGRASE = '#E65100'   # Deep orange
C_MERCURY   = '#37474F'   # Dark teal
C_OTHER     = '#B0BEC5'   # Blue-grey light
C_HYPO      = '#CFD8DC'   # Very light blue-grey

# ============================================================
# PARSE GENES FROM BAKTA TSV
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
                'strand': 1 if parts[4] == '+' else -1,
                'gene': parts[6],
                'product': parts[7],
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
    
    amr_kw = ['bla', 'aph(', 'aac(', 'ant(', 'rmtb', 'arma', 'fosa', 'cata', 'catb',
              'sul1', 'sul2', 'dfra', 'dfrb', 'mph(', 'msr(', 'erm', 'qnr', 'oqx',
              'tet(', 'mcr', 'flor', 'cmla', 'ble', 'qac']
    for kw in amr_kw:
        if kw in gene:
            return 'amr'
    if any(k in product for k in ['beta-lactamase', 'aminoglycoside', 'chloramphenicol',
                                   'fosfomycin', 'sulfonamide', 'dihydrofolate', 'bleomycin',
                                   'macrolide', 'phosphotransferase']):
        if 'transpos' not in product:
            return 'amr'
    
    if gene == 'tnp' or 'transpos' in product or 'tn3' in product.lower():
        return 'is'
    if any(k in product.lower() for k in ['is1 ', 'is5 ', 'is6 ', 'is26', 'is110', 'is481',
                                            'is903', 'is4321', 'isncy', 'is91', 'iskpn',
                                            'isec', 'iscfr', 'iskra', 'is5075', 'is15',
                                            'mutator family']):
        return 'is'
    
    if 'inti' in gene or 'integrase' in product or 'integron' in product:
        return 'integrase'
    if 'recombinas' in product or 'resolvas' in product or 'invertase' in product:
        return 'integrase'
    
    if gene.startswith('mer') or 'mercury' in product or 'tellurium' in product:
        return 'mercury'
    
    if 'hypothetical' in product:
        return 'hypo'
    
    return 'other'

def get_color(cat):
    return {'amr': C_AMR, 'is': C_IS, 'integrase': C_INTEGRASE,
            'mercury': C_MERCURY, 'other': C_OTHER, 'hypo': C_HYPO}.get(cat, C_OTHER)

def get_label(g):
    cat = classify(g)
    gene = g['gene']
    product = g['product']
    
    if cat == 'amr':
        return gene if gene else product[:15]
    elif cat == 'is':
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
        if gene == 'tnp':
            return 'tnp'
        return 'IS'
    elif cat == 'integrase':
        if 'intI1' in gene or 'intI1' in product:
            return 'intI1'
        if 'resolv' in product.lower():
            return 'res'
        return 'int'
    elif cat == 'mercury':
        return gene if gene else 'mer'
    return ''

# ============================================================
# DEFINE CASSETTE REGIONS
# ============================================================
cassettes = [
    {
        'name': 'A', 'genes': kpc_genes, 'start': 2500, 'end': 13000,
        'title': 'pKPC (IncFII+IncR, ~131 kb) — blaCTX-M-65 / fosA3 region',
        'tn_regions': [(3300, 7600, 'IS1B–blaCTX-M-65–IS903B'),
                       (9200, 12100, 'IS15DIV–fosA3–IS15DIV')],
    },
    {
        'name': 'B', 'genes': kpc_genes, 'start': 31000, 'end': 52000,
        'title': 'pKPC (IncFII+IncR, ~131 kb) — blaSHV-12 / blaKPC-2 region',
        'tn_regions': [(32300, 37600, 'IS15DIV–blaSHV-12'),
                       (39700, 43300, 'Tn3–ISKpn27–blaKPC-2'),
                       (49000, 51100, 'IS15DIV–int')],
    },
    {
        'name': 'C', 'genes': kpc_genes, 'start': 51000, 'end': 67000,
        'title': 'pKPC (IncFII+IncR, ~131 kb) — mer operon / blaTEM-1 / rmtB1 region',
        'tn_regions': [(52000, 57000, 'mer operon–IS5075'),
                       (58900, 65200, 'IS15DIV–blaTEM-1–rmtB1–IS15')],
    },
    {
        'name': 'D', 'genes': kpc_genes, 'start': 71000, 'end': 82000,
        'title': 'pKPC (IncFII+IncR, ~131 kb) — catA2 region',
        'tn_regions': [(72000, 81200, 'IS15DIV–ISCfr3–ISKpn26–catA2–IS15DIV')],
    },
    {
        'name': 'E', 'genes': ndm_genes, 'start': 196000, 'end': 220000,
        'title': 'pNDM (IncFIB+IncHI1B, ~354 kb) — Class 1 integron / blaNDM-1 / armA region',
        'tn_regions': [(198000, 205600, 'Class 1 integron (sul1–dfrA5–intI1)'),
                       (205700, 213000, 'IS15DIV–aph(3\')-Ia–ble–blaNDM-1'),
                       (213200, 219000, 'armA–ISEc29–msr(E)–mph(E)')],
    },
    {
        'name': 'F', 'genes': ndm_genes, 'start': 220000, 'end': 243000,
        'title': 'pNDM (IncFIB+IncHI1B, ~354 kb) — qnrS1 / aph(3\')-VI / sul2 region',
        'tn_regions': [(222700, 230100, 'IS15DIV–qnrS1–ISKpn19'),
                       (231700, 242100, 'aph(3\')-VI–MuA–sul2–ISVsa3')],
    },
]

# ============================================================
# CREATE EACH PANEL
# ============================================================
panel_figs = []

for cas in cassettes:
    region_start = cas['start']
    region_end = cas['end']
    region_len = region_end - region_start
    
    region_genes = [g for g in cas['genes'] 
                    if g['end'] >= region_start and g['start'] <= region_end]
    
    gv = GenomeViz(
        fig_width=24,
        fig_track_height=1.2,
        feature_track_ratio=0.3,
        theme='light',
        show_axis=False,
    )
    
    track = gv.add_feature_track(
        cas['title'],
        region_end - region_start,
        labelsize=13,
    )
    
    for g in region_genes:
        cat = classify(g)
        color = get_color(cat)
        label = get_label(g)
        s = max(g['start'], region_start) - region_start
        e = min(g['end'], region_end) - region_start
        
        if cat == 'hypo' and (e - s) < region_len * 0.003:
            continue
        
        # ENHANCED FONT SIZES for A+++++ quality
        if cat == 'amr':
            text_kws = dict(size=14, color=C_AMR, style='italic', weight='bold')
            plotstyle = 'bigarrow'
        elif cat == 'is':
            text_kws = dict(size=10, color=C_IS, style='normal', weight='bold')
            plotstyle = 'bigarrow'
        elif cat == 'integrase':
            text_kws = dict(size=10, color=C_INTEGRASE, style='italic', weight='bold')
            plotstyle = 'bigarrow'
        elif cat == 'mercury':
            text_kws = dict(size=9, color=C_MERCURY, style='italic', weight='normal')
            plotstyle = 'bigarrow'
        elif cat == 'hypo':
            label = ''
            text_kws = dict(size=6, color='#BDBDBD')
            plotstyle = 'bigarrow'
        else:
            if label:
                text_kws = dict(size=8, color='#546E7A', style='italic')
            else:
                text_kws = dict(size=6, color='#BDBDBD')
            plotstyle = 'bigarrow'
        
        alpha = 0.25 if cat == 'hypo' else 0.92
        
        track.add_feature(
            s, e,
            strand=g['strand'],
            plotstyle=plotstyle,
            arrow_shaft_ratio=0.5,
            label=label,
            text_kws=text_kws,
            fc=color,
            ec='#212121' if cat in ['amr', 'is', 'integrase'] else '#BDBDBD',
            lw=1.0 if cat in ['amr', 'is', 'integrase'] else 0.3,
            alpha=alpha,
        )
    
    if region_len <= 12000:
        gv.set_scale_bar(scale_size_label=(1000, "1 kb"), labelsize=11)
    else:
        gv.set_scale_bar(scale_size_label=(5000, "5 kb"), labelsize=11)
    
    fig = gv.plotfig()
    
    # Add transposon cassette brackets
    ax_list = fig.axes
    if ax_list:
        ax = ax_list[0]
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        y_range = ylim[1] - ylim[0]
        
        for tn_start, tn_end, tn_label in cas['tn_regions']:
            ts = tn_start - region_start
            te = tn_end - region_start
            
            # Shaded background for transposon cassette
            ax.axvspan(ts, te, alpha=0.07, color='#FF8F00', zorder=0)
            
            # Bracket at top
            bracket_y = ylim[1] - y_range * 0.03
            tick_h = y_range * 0.04
            ax.plot([ts, ts, te, te],
                    [bracket_y - tick_h, bracket_y, bracket_y, bracket_y - tick_h],
                    color='#FF8F00', linewidth=2.0, zorder=10, clip_on=False,
                    solid_capstyle='round')
            
            # Bracket label
            ax.text((ts + te) / 2, bracket_y + tick_h * 0.3, tn_label,
                    fontsize=8, fontweight='bold', fontstyle='italic',
                    color='#E65100', ha='center', va='bottom', zorder=11,
                    clip_on=False)
    
    # Panel label
    if ax_list:
        ax = ax_list[0]
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.text(xlim[0] - (xlim[1] - xlim[0]) * 0.04, ylim[1],
                cas['name'], fontsize=24, fontweight='bold',
                fontfamily='serif', ha='right', va='top', zorder=12,
                clip_on=False)
    
    # Save individual panel PNG at 300 DPI
    panel_png = OUT_DIR / f"Panel_{cas['name']}_cassette_v2.png"
    fig.savefig(str(panel_png), dpi=300, bbox_inches='tight', facecolor='white')
    
    # Save individual panel PDF (vector)
    panel_pdf = OUT_DIR / f"Panel_{cas['name']}_cassette_v2.pdf"
    fig.savefig(str(panel_pdf), bbox_inches='tight', facecolor='white')
    
    plt.close(fig)
    print(f"Panel {cas['name']} saved: {panel_png}")

# ============================================================
# COMBINE ALL PANELS INTO ONE FIGURE using PIL (memory efficient)
# ============================================================
from PIL import Image, ImageDraw, ImageFont
import numpy as np

panels = []
for cas in cassettes:
    p = Image.open(str(OUT_DIR / f"Panel_{cas['name']}_cassette_v2.png"))
    panels.append(p)

max_width = max(p.width for p in panels)
total_height = sum(p.height for p in panels)

# Title and legend space
title_h = 120
legend_h = 100
gap = 10
combined_h = title_h + total_height + len(panels) * gap + legend_h

combined = Image.new('RGB', (max_width, combined_h), 'white')
draw = ImageDraw.Draw(combined)

# Title text
try:
    font_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf", 42)
    font_sub = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf", 30)
    font_leg = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 24)
except:
    font_title = ImageFont.load_default()
    font_sub = ImageFont.load_default()
    font_leg = ImageFont.load_default()

title_text = "Genetic Environment of Antimicrobial Resistance Gene Cassettes"
sub_text = "ST11 Klebsiella pneumoniae KP21915"

# Center title
bbox_t = draw.textbbox((0, 0), title_text, font=font_title)
tw = bbox_t[2] - bbox_t[0]
draw.text(((max_width - tw) // 2, 15), title_text, fill='#212121', font=font_title)

bbox_s = draw.textbbox((0, 0), sub_text, font=font_sub)
sw = bbox_s[2] - bbox_s[0]
draw.text(((max_width - sw) // 2, 70), sub_text, fill='#424242', font=font_sub)

# Paste panels
y_pos = title_h
for panel in panels:
    x_offset = (max_width - panel.width) // 2
    combined.paste(panel, (x_offset, y_pos))
    y_pos += panel.height + gap

# Legend
legend_y = y_pos + 5
legend_items = [
    (C_AMR, "Antimicrobial resistance gene"),
    (C_IS, "IS element / Transposase"),
    (C_INTEGRASE, "Integrase / Recombinase"),
    (C_MERCURY, "Heavy metal resistance"),
    ('#B0BEC5', "Other CDS"),
]

x_pos = 100
for color, label in legend_items:
    # Draw color box
    draw.rectangle([x_pos, legend_y + 5, x_pos + 30, legend_y + 30], fill=color, outline='#424242')
    draw.text((x_pos + 38, legend_y + 5), label, fill='#212121', font=font_leg)
    x_pos += 38 + len(label) * 14 + 40

# Transposon bracket legend
draw.rectangle([x_pos, legend_y + 5, x_pos + 30, legend_y + 30], fill='#FFF3E0', outline='#FF8F00', width=2)
draw.text((x_pos + 38, legend_y + 5), "Transposon cassette", fill='#212121', font=font_leg)

# Save combined
combined_png = OUT_DIR / "Figure_MGE_Cassettes_pyGenomeViz_v2.png"
combined.save(str(combined_png), dpi=(300, 300))
print(f"\nCombined figure saved: {combined_png}")
print(f"Size: {combined.width}x{combined.height} px")

# Close panels to free memory
for p in panels:
    p.close()
combined.close()

print("\n=== ALL DONE ===")
print("Tool: pyGenomeViz v1.6.1")
print("Citation: Shimoyama Y. pyGenomeViz: A genome visualization python package")
print("          for comparative genomics. https://github.com/moshi4/pyGenomeViz")
