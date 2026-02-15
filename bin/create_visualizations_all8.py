#!/usr/bin/env python3
"""
Create publication-ready visualizations for ST11 K. pneumoniae plasmid analysis.
All 8 genomes included. Designed for Lancet-level quality.
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
import numpy as np
from collections import defaultdict, Counter
from pathlib import Path
import json
import csv
import io

# Apply style first, then set fonts
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.2
})

BASE = Path("/home/ubuntu/plasmid_analysis")
OUT_DIR = BASE / "07_Visualizations"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = ['JKPB244', 'JKPB284', 'JKPB381', 'JKPR282', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']

# Color palette - professional, journal-quality
COLORS = {
    'KPC': '#E63946',      # Red
    'NDM': '#457B9D',      # Blue
    'OXA': '#F4A261',      # Orange
    'CTX-M': '#2A9D8F',    # Teal
    'SHV': '#264653',      # Dark blue
    'TEM': '#E9C46A',      # Yellow
    'chromosome': '#6C757D',  # Gray
    'plasmid': '#E63946',     # Red
    'present': '#2A9D8F',     # Teal
    'absent': '#F8F9FA',      # Light gray
}

def load_amr_data():
    """Load AMR data from CSV."""
    amr_csv = BASE / "04_AMR_Plasmid_Association" / "amr_genes_all8_samples.csv"
    df = pd.read_csv(amr_csv)
    df['Gene'] = df['Gene'].fillna('')
    df['Product'] = df['Product'].fillna('')
    return df

def load_mge_data():
    """Load MGE data from CSV."""
    mge_csv = BASE / "04_AMR_Plasmid_Association" / "mge_elements_all8_samples.csv"
    df = pd.read_csv(mge_csv)
    df['Gene'] = df['Gene'].fillna('')
    df['Product'] = df['Product'].fillna('')
    return df

def load_integron_data():
    """Load IntegronFinder results."""
    integron_dir = BASE / "05_MGE_Analysis" / "IntegronFinder_Results"
    results = {}
    if integron_dir.exists():
        for sample_dir in integron_dir.iterdir():
            if sample_dir.is_dir():
                sample = sample_dir.name
                summary_file = sample_dir / f"{sample}.summary"
                if summary_file.exists():
                    with open(summary_file) as f:
                        content = f.read()
                    results[sample] = content
    return results

# ============================================================
# FIGURE 1: Carbapenemase Distribution (Stacked Bar)
# ============================================================
def figure1_carbapenemase_distribution(amr_df):
    """Create a carbapenemase distribution chart across all 8 isolates."""
    fig, ax = plt.subplots(figsize=(10, 5))
    
    # Key beta-lactamase genes
    target_genes = {
        'blaKPC': 'KPC',
        'blaNDM': 'NDM',
        'blaOXA': 'OXA',
        'blaCTX-M': 'CTX-M',
        'blaSHV': 'SHV',
        'blaTEM': 'TEM'
    }
    
    # Build presence/absence matrix
    matrix = {}
    for prefix, label in target_genes.items():
        matrix[label] = []
        for sample in SAMPLES:
            sample_df = amr_df[amr_df['Sample'] == sample]
            present = any(sample_df['Gene'].str.startswith(prefix))
            matrix[label].append(1 if present else 0)
    
    # Plot as grouped bar chart
    x = np.arange(len(SAMPLES))
    width = 0.12
    offsets = np.arange(len(target_genes)) - (len(target_genes) - 1) / 2
    
    gene_colors = ['#E63946', '#457B9D', '#F4A261', '#2A9D8F', '#264653', '#E9C46A']
    
    for i, (label, values) in enumerate(matrix.items()):
        bars = ax.bar(x + offsets[i] * width, values, width, 
                      label=label, color=gene_colors[i], edgecolor='white', linewidth=0.5)
    
    ax.set_xlabel('Isolate', fontweight='bold')
    ax.set_ylabel('Presence (1) / Absence (0)', fontweight='bold')
    ax.set_title('Distribution of Key Beta-Lactamase Genes Across ST11 K. pneumoniae Isolates',
                 fontweight='bold', fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(SAMPLES, rotation=45, ha='right')
    ax.set_ylim(0, 1.3)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['Absent', 'Present'])
    ax.legend(loc='upper right', ncol=3, framealpha=0.9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    fig.savefig(OUT_DIR / 'Figure1_Carbapenemase_Distribution.png')
    fig.savefig(OUT_DIR / 'Figure1_Carbapenemase_Distribution.pdf')
    plt.close()
    print("Figure 1: Carbapenemase Distribution - DONE")

# ============================================================
# FIGURE 2: AMR Gene Heatmap
# ============================================================
def figure2_amr_heatmap(amr_df):
    """Create a comprehensive AMR gene presence/absence heatmap."""
    # Focus on clinically relevant AMR genes
    key_genes = [
        'blaKPC-2', 'blaNDM-1', 'blaOXA-1', 'blaCTX-M-65', 'blaSHV-12', 'blaTEM-1',
        'rmtB', 'armA', 'aac(6\')-Ib-cr', 'aph(3\')-Ia', 'aadA2',
        'qnrS1', 'oqxA', 'oqxB',
        'fosA', 'fosA3',
        'sul1', 'sul2', 'dfrA5', 'dfrA12',
        'tet(A)', 'tet(D)',
        'mph(A)', 'mph(E)', 'msr(E)', 'ere(A)',
        'catA2', 'floR',
        'ble',
        'arr-3'
    ]
    
    # Build matrix
    matrix_data = {}
    for gene in key_genes:
        row = []
        for sample in SAMPLES:
            sample_df = amr_df[amr_df['Sample'] == sample]
            # Check exact match or starts-with for variants
            present = any(sample_df['Gene'].str.contains(gene.replace('(', '\\(').replace(')', '\\)').replace('\'', '\\\''), regex=True, na=False))
            if not present:
                # Try simpler match
                present = any(sample_df['Gene'] == gene)
            row.append(1 if present else 0)
        matrix_data[gene] = row
    
    # Create DataFrame
    heatmap_df = pd.DataFrame(matrix_data, index=SAMPLES).T
    
    # Remove genes not found in any sample
    heatmap_df = heatmap_df[heatmap_df.sum(axis=1) > 0]
    
    # Categorize genes by drug class
    drug_classes = {
        'Carbapenems': ['blaKPC-2', 'blaNDM-1'],
        'ESBL': ['blaCTX-M-65', 'blaSHV-12', 'blaTEM-1'],
        'Other β-lactams': ['blaOXA-1'],
        'Aminoglycosides': ['rmtB', 'armA', 'aac(6\')-Ib-cr', 'aph(3\')-Ia', 'aadA2'],
        'Quinolones': ['qnrS1', 'oqxA', 'oqxB'],
        'Fosfomycin': ['fosA', 'fosA3'],
        'Sulfonamides/Trimethoprim': ['sul1', 'sul2', 'dfrA5', 'dfrA12'],
        'Tetracyclines': ['tet(A)', 'tet(D)'],
        'Macrolides': ['mph(A)', 'mph(E)', 'msr(E)', 'ere(A)'],
        'Phenicols': ['catA2', 'floR'],
        'Bleomycin': ['ble'],
        'Rifampicin': ['arr-3']
    }
    
    # Reorder by drug class
    ordered_genes = []
    class_labels = []
    for cls, genes in drug_classes.items():
        for g in genes:
            if g in heatmap_df.index:
                ordered_genes.append(g)
                class_labels.append(cls)
    
    heatmap_df = heatmap_df.loc[[g for g in ordered_genes if g in heatmap_df.index]]
    
    fig, ax = plt.subplots(figsize=(10, max(8, len(heatmap_df) * 0.35)))
    
    # Custom colormap
    cmap = matplotlib.colors.ListedColormap(['#F0F0F0', '#2A9D8F'])
    
    sns.heatmap(heatmap_df, cmap=cmap, linewidths=0.5, linecolor='white',
                cbar=False, ax=ax, square=True,
                xticklabels=True, yticklabels=True)
    
    ax.set_xlabel('Isolate', fontweight='bold')
    ax.set_ylabel('AMR Gene', fontweight='bold')
    ax.set_title('AMR Gene Presence/Absence Across ST11 K. pneumoniae Isolates',
                 fontweight='bold', fontsize=11, pad=15)
    
    # Add legend
    present_patch = mpatches.Patch(color='#2A9D8F', label='Present')
    absent_patch = mpatches.Patch(color='#F0F0F0', label='Absent')
    ax.legend(handles=[present_patch, absent_patch], loc='upper right',
              bbox_to_anchor=(1.15, 1.0))
    
    plt.tight_layout()
    fig.savefig(OUT_DIR / 'Figure2_AMR_Heatmap.png')
    fig.savefig(OUT_DIR / 'Figure2_AMR_Heatmap.pdf')
    plt.close()
    print("Figure 2: AMR Heatmap - DONE")

# ============================================================
# FIGURE 3: Plasmid Architecture Comparison
# ============================================================
def figure3_plasmid_architecture():
    """Create a plasmid architecture comparison across all 8 isolates."""
    # Load enumeration data
    enum_file = BASE / "01_plasmid_enumeration" / "plasmid_enumeration_summary.csv"
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plasmid data from our analysis (contig sizes from FNA files)
    # We'll extract from the Bakta annotations
    plasmid_data = {}
    for sample in SAMPLES:
        sample_dir = BASE / "04_AMR_Plasmid_Association" / "Bakta_Annotations" / sample
        fna_file = sample_dir / f"{sample}.fna"
        json_file = sample_dir / f"{sample}.json"
        
        contigs = []
        if json_file.exists():
            with open(json_file) as f:
                data = json.load(f)
            for seq in data.get('sequences', []):
                contigs.append({'name': seq['id'], 'length': seq.get('length', 0)})
        elif fna_file.exists():
            current_contig = None
            current_size = 0
            with open(fna_file) as f:
                for line in f:
                    if line.startswith('>'):
                        if current_contig:
                            contigs.append({'name': current_contig, 'length': current_size})
                        current_contig = line.strip().split()[0][1:]
                        current_size = 0
                    else:
                        current_size += len(line.strip())
                if current_contig:
                    contigs.append({'name': current_contig, 'length': current_size})
        
        # Sort by size descending
        contigs.sort(key=lambda x: x['length'], reverse=True)
        plasmid_data[sample] = contigs
    
    # Plot horizontal bars for each sample
    y_positions = np.arange(len(SAMPLES))
    bar_height = 0.6
    
    colors_by_type = {
        'chromosome': '#6C757D',
        'large_plasmid': '#E63946',
        'medium_plasmid': '#457B9D',
        'small_plasmid': '#2A9D8F',
        'small_replicon': '#E9C46A'
    }
    
    for i, sample in enumerate(SAMPLES):
        contigs = plasmid_data.get(sample, [])
        x_offset = 0
        for c in contigs:
            size_kb = c['length'] / 1000
            if c['length'] > 4000000:
                color = colors_by_type['chromosome']
                label_type = 'Chromosome'
            elif c['length'] > 200000:
                color = colors_by_type['large_plasmid']
                label_type = 'Large plasmid'
            elif c['length'] > 80000:
                color = colors_by_type['medium_plasmid']
                label_type = 'Medium plasmid'
            elif c['length'] > 30000:
                color = colors_by_type['small_plasmid']
                label_type = 'Small plasmid'
            else:
                color = colors_by_type['small_replicon']
                label_type = 'Small replicon'
            
            # Plot bar (use log scale for visibility)
            bar_width = np.log10(max(size_kb, 1)) * 50
            ax.barh(i, bar_width, height=bar_height, left=x_offset,
                    color=color, edgecolor='white', linewidth=0.5)
            
            # Add size label for significant contigs
            if size_kb > 50:
                ax.text(x_offset + bar_width/2, i, f'{size_kb:.0f}kb',
                       ha='center', va='center', fontsize=6, fontweight='bold', color='white')
            
            x_offset += bar_width + 2
    
    ax.set_yticks(y_positions)
    ax.set_yticklabels(SAMPLES)
    ax.set_xlabel('Relative Size (log scale)', fontweight='bold')
    ax.set_title('Genome Architecture: Chromosome and Plasmid Content of ST11 K. pneumoniae Isolates',
                 fontweight='bold', fontsize=11)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Legend
    legend_patches = [
        mpatches.Patch(color=colors_by_type['chromosome'], label='Chromosome'),
        mpatches.Patch(color=colors_by_type['large_plasmid'], label='Large plasmid (>200kb)'),
        mpatches.Patch(color=colors_by_type['medium_plasmid'], label='Medium plasmid (80-200kb)'),
        mpatches.Patch(color=colors_by_type['small_plasmid'], label='Small plasmid (30-80kb)'),
        mpatches.Patch(color=colors_by_type['small_replicon'], label='Small replicon (<30kb)')
    ]
    ax.legend(handles=legend_patches, loc='lower right', framealpha=0.9)
    
    plt.tight_layout()
    fig.savefig(OUT_DIR / 'Figure3_Plasmid_Architecture.png')
    fig.savefig(OUT_DIR / 'Figure3_Plasmid_Architecture.pdf')
    plt.close()
    print("Figure 3: Plasmid Architecture - DONE")

# ============================================================
# FIGURE 4: MGE Distribution
# ============================================================
def figure4_mge_distribution(mge_df):
    """Create MGE distribution visualization."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel A: MGE counts by type per sample
    mge_counts = defaultdict(lambda: defaultdict(int))
    for _, row in mge_df.iterrows():
        sample = row['Sample']
        product = str(row.get('Product', '') or '')
        gene = str(row.get('Gene', '') or '')
        
        if 'transposase' in product.lower():
            mge_counts[sample]['Transposase'] += 1
        elif 'integrase' in product.lower():
            mge_counts[sample]['Integrase'] += 1
        elif 'recombinase' in product.lower():
            mge_counts[sample]['Recombinase'] += 1
        elif 'IS' in gene:
            mge_counts[sample]['IS Element'] += 1
        elif 'phage' in product.lower():
            mge_counts[sample]['Phage-related'] += 1
        elif 'conjugat' in product.lower() or 'relaxase' in product.lower():
            mge_counts[sample]['Conjugation'] += 1
        else:
            mge_counts[sample]['Other MGE'] += 1
    
    mge_types = ['Transposase', 'Integrase', 'Recombinase', 'IS Element', 'Phage-related', 'Conjugation', 'Other MGE']
    mge_colors = ['#E63946', '#457B9D', '#2A9D8F', '#264653', '#E9C46A', '#F4A261', '#ADB5BD']
    
    x = np.arange(len(SAMPLES))
    bottom = np.zeros(len(SAMPLES))
    
    for mge_type, color in zip(mge_types, mge_colors):
        values = [mge_counts[s].get(mge_type, 0) for s in SAMPLES]
        ax1.bar(x, values, bottom=bottom, label=mge_type, color=color, edgecolor='white', linewidth=0.3)
        bottom += values
    
    ax1.set_xlabel('Isolate', fontweight='bold')
    ax1.set_ylabel('Number of MGE Elements', fontweight='bold')
    ax1.set_title('A. MGE Element Distribution by Type', fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(SAMPLES, rotation=45, ha='right')
    ax1.legend(loc='upper right', fontsize=7, framealpha=0.9)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Panel B: IS element families
    is_families = defaultdict(lambda: defaultdict(int))
    for _, row in mge_df.iterrows():
        gene = str(row.get('Gene', '') or '')
        sample = row['Sample']
        if gene.startswith('IS'):
            # Extract IS family name
            family = gene.split('_')[0] if '_' in gene else gene
            is_families[sample][family] += 1
    
    # Get top IS families
    all_families = Counter()
    for sample_data in is_families.values():
        all_families.update(sample_data)
    
    top_families = [f for f, _ in all_families.most_common(8)]
    
    if top_families:
        is_colors = plt.cm.Set2(np.linspace(0, 1, len(top_families)))
        bottom2 = np.zeros(len(SAMPLES))
        
        for fam, color in zip(top_families, is_colors):
            values = [is_families[s].get(fam, 0) for s in SAMPLES]
            ax2.bar(x, values, bottom=bottom2, label=fam, color=color, edgecolor='white', linewidth=0.3)
            bottom2 += values
        
        ax2.set_xlabel('Isolate', fontweight='bold')
        ax2.set_ylabel('Count', fontweight='bold')
        ax2.set_title('B. IS Element Families', fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(SAMPLES, rotation=45, ha='right')
        ax2.legend(loc='upper right', fontsize=7, framealpha=0.9)
    else:
        ax2.text(0.5, 0.5, 'No IS families detected', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('B. IS Element Families', fontweight='bold')
    
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    fig.savefig(OUT_DIR / 'Figure4_MGE_Distribution.png')
    fig.savefig(OUT_DIR / 'Figure4_MGE_Distribution.pdf')
    plt.close()
    print("Figure 4: MGE Distribution - DONE")

# ============================================================
# FIGURE 5: AMR Gene Class Distribution
# ============================================================
def figure5_amr_class_distribution(amr_df):
    """Create AMR gene class distribution visualization."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Classify AMR genes by drug class
    def classify_amr(gene, product):
        gene = str(gene or '').lower()
        product = str(product or '').lower()
        
        if 'carbapenem' in product or gene.startswith('blakpc') or gene.startswith('blandm') or gene.startswith('blavim') or gene.startswith('blaimp'):
            return 'Carbapenems'
        elif 'beta-lactam' in product or gene.startswith('bla'):
            return 'Beta-lactams'
        elif 'aminoglycoside' in product or gene.startswith('aac') or gene.startswith('aph') or gene.startswith('ant') or gene.startswith('rmt') or gene.startswith('arm'):
            return 'Aminoglycosides'
        elif 'quinolone' in product or gene.startswith('qnr') or gene.startswith('oqx'):
            return 'Quinolones'
        elif 'macrolide' in product or gene.startswith('mph') or gene.startswith('msr') or gene.startswith('ere') or gene.startswith('erm'):
            return 'Macrolides'
        elif 'tetracycline' in product or gene.startswith('tet'):
            return 'Tetracyclines'
        elif 'sulfonamide' in product or gene.startswith('sul') or gene.startswith('dfr'):
            return 'Sulfonamides/Trimethoprim'
        elif 'fosfomycin' in product or gene.startswith('fos'):
            return 'Fosfomycin'
        elif 'chloramphenicol' in product or 'phenicol' in product or gene.startswith('cat') or gene.startswith('flo') or gene.startswith('cml'):
            return 'Phenicols'
        elif 'efflux' in product or 'pump' in product:
            return 'Efflux pumps'
        elif 'resistance' in product:
            return 'Other resistance'
        else:
            return 'Other'
    
    amr_df_copy = amr_df.copy()
    amr_df_copy['Drug_Class'] = amr_df_copy.apply(lambda r: classify_amr(r['Gene'], r['Product']), axis=1)
    
    # Panel A: Stacked bar by drug class per sample
    class_counts = amr_df_copy.groupby(['Sample', 'Drug_Class']).size().unstack(fill_value=0)
    
    # Order drug classes by total count
    class_order = class_counts.sum().sort_values(ascending=False).index.tolist()
    class_counts = class_counts[class_order]
    
    class_colors = plt.cm.tab20(np.linspace(0, 1, len(class_order)))
    
    class_counts.loc[SAMPLES].plot(kind='bar', stacked=True, ax=ax1, color=class_colors, edgecolor='white', linewidth=0.3)
    ax1.set_xlabel('Isolate', fontweight='bold')
    ax1.set_ylabel('Number of AMR Genes', fontweight='bold')
    ax1.set_title('A. AMR Gene Distribution by Drug Class', fontweight='bold')
    ax1.set_xticklabels(SAMPLES, rotation=45, ha='right')
    ax1.legend(loc='upper right', fontsize=6, framealpha=0.9, ncol=2)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Panel B: Pie chart of overall drug class distribution
    overall = amr_df_copy['Drug_Class'].value_counts()
    # Combine small categories
    threshold = len(amr_df_copy) * 0.03
    major = overall[overall >= threshold]
    minor_sum = overall[overall < threshold].sum()
    if minor_sum > 0:
        major['Other (combined)'] = minor_sum
    
    pie_colors = plt.cm.tab20(np.linspace(0, 1, len(major)))
    wedges, texts, autotexts = ax2.pie(major.values, labels=None, autopct='%1.0f%%',
                                        colors=pie_colors, pctdistance=0.75,
                                        startangle=90)
    ax2.legend(major.index, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7)
    ax2.set_title('B. Overall AMR Gene Class Distribution', fontweight='bold')
    
    for autotext in autotexts:
        autotext.set_fontsize(7)
    
    plt.tight_layout()
    fig.savefig(OUT_DIR / 'Figure5_AMR_Class_Distribution.png')
    fig.savefig(OUT_DIR / 'Figure5_AMR_Class_Distribution.pdf')
    plt.close()
    print("Figure 5: AMR Class Distribution - DONE")

# ============================================================
# FIGURE 6: Dual Carbapenemase Summary
# ============================================================
def figure6_dual_carbapenemase(amr_df):
    """Create a summary figure showing dual carbapenemase producer status."""
    fig, ax = plt.subplots(figsize=(10, 4))
    
    carb_genes = ['blaKPC', 'blaNDM', 'blaOXA', 'blaVIM', 'blaIMP']
    carb_labels = ['KPC', 'NDM', 'OXA', 'VIM', 'IMP']
    
    # Build matrix
    matrix = []
    for sample in SAMPLES:
        row = []
        sample_df = amr_df[amr_df['Sample'] == sample]
        for prefix in carb_genes:
            present = any(sample_df['Gene'].str.startswith(prefix, na=False))
            row.append(1 if present else 0)
        matrix.append(row)
    
    matrix = np.array(matrix)
    
    # Only show genes that appear in at least one sample
    mask = matrix.sum(axis=0) > 0
    matrix = matrix[:, mask]
    labels = [l for l, m in zip(carb_labels, mask) if m]
    
    # Custom colormap
    cmap = matplotlib.colors.ListedColormap(['#F0F0F0', '#E63946'])
    
    im = ax.imshow(matrix.T, cmap=cmap, aspect='auto')
    
    ax.set_xticks(range(len(SAMPLES)))
    ax.set_xticklabels(SAMPLES, rotation=45, ha='right')
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontweight='bold')
    
    # Add text annotations
    for i in range(len(labels)):
        for j in range(len(SAMPLES)):
            text = '+' if matrix[j, i] else '-'
            color = 'white' if matrix[j, i] else '#CCC'
            ax.text(j, i, text, ha='center', va='center', fontweight='bold',
                   fontsize=14, color=color)
    
    # Add dual producer annotation
    dual_count = 0
    for j, sample in enumerate(SAMPLES):
        sample_df = amr_df[amr_df['Sample'] == sample]
        has_kpc = any(sample_df['Gene'].str.startswith('blaKPC', na=False))
        has_ndm = any(sample_df['Gene'].str.startswith('blaNDM', na=False))
        if has_kpc and has_ndm:
            dual_count += 1
            ax.text(j, -0.6, 'DUAL', ha='center', va='center', fontsize=7,
                   fontweight='bold', color='#E63946',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='#FFF3F3', edgecolor='#E63946', linewidth=0.5))
    
    ax.set_title(f'Carbapenemase Gene Distribution: {dual_count}/8 Dual KPC+NDM Producers',
                 fontweight='bold', fontsize=11, pad=20)
    
    ax.set_xlim(-0.5, len(SAMPLES) - 0.5)
    ax.set_ylim(len(labels) - 0.5, -0.8)
    
    # Legend
    present_patch = mpatches.Patch(color='#E63946', label='Present')
    absent_patch = mpatches.Patch(color='#F0F0F0', label='Absent')
    ax.legend(handles=[present_patch, absent_patch], loc='upper right')
    
    plt.tight_layout()
    fig.savefig(OUT_DIR / 'Figure6_Dual_Carbapenemase.png')
    fig.savefig(OUT_DIR / 'Figure6_Dual_Carbapenemase.pdf')
    plt.close()
    print("Figure 6: Dual Carbapenemase Summary - DONE")

# ============================================================
# MAIN
# ============================================================
def main():
    print("Loading data...")
    amr_df = load_amr_data()
    mge_df = load_mge_data()
    
    print(f"AMR records: {len(amr_df)}")
    print(f"MGE records: {len(mge_df)}")
    print()
    
    print("Creating publication-ready figures...")
    figure1_carbapenemase_distribution(amr_df)
    figure2_amr_heatmap(amr_df)
    figure3_plasmid_architecture()
    figure4_mge_distribution(mge_df)
    figure5_amr_class_distribution(amr_df)
    figure6_dual_carbapenemase(amr_df)
    
    print(f"\nAll figures saved to: {OUT_DIR}")
    print("Files generated:")
    for f in sorted(OUT_DIR.glob('Figure*')):
        print(f"  {f.name}")

if __name__ == '__main__':
    main()
