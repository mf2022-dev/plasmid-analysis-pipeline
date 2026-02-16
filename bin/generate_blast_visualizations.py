#!/usr/bin/env python3
"""
Publication-quality visualizations for global KPC-2 cassette dissemination analysis.
Generates Lancet-level figures for manuscript.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns
import json
import os

# ============================================================
# STYLE CONFIGURATION - Lancet-level publication quality
# ============================================================
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'axes.edgecolor': '#333333',
    'grid.linewidth': 0.4,
    'grid.alpha': 0.5,
})

# Color palette - consistent with previous figures
COLORS = {
    'kp': '#C0392B',        # Red for K. pneumoniae
    'ecoli': '#2980B9',     # Blue for E. coli
    'serratia': '#27AE60',  # Green for Serratia
    'proteus': '#F39C12',   # Orange for Proteus
    'pseudo': '#8E44AD',    # Purple for Pseudomonas
    'klebsiella_other': '#E74C3C',  # Light red for other Klebsiella
    'citrobacter': '#16A085',       # Teal for Citrobacter
    'enterobacter': '#D35400',      # Dark orange for Enterobacter
    'other': '#95A5A6',     # Gray for others
    'plasmid': '#C0392B',   # Red for plasmid
    'chromosome': '#2980B9', # Blue for chromosome
    'china': '#E74C3C',
    'taiwan': '#3498DB',
    'sudan': '#2ECC71',
    'japan': '#9B59B6',
    'usa': '#F1C40F',
    'india': '#E67E22',
    'other_country': '#95A5A6',
}

OUTDIR = '/home/ubuntu/plasmid_analysis/10_BLAST_Results/figures'
os.makedirs(OUTDIR, exist_ok=True)

def save_fig(fig, name):
    """Save figure in PNG, PDF, and SVG formats."""
    for fmt in ['png', 'pdf', 'svg']:
        fig.savefig(f'{OUTDIR}/{name}.{fmt}', format=fmt, dpi=300, bbox_inches='tight')
    print(f"  Saved: {name} (PNG/PDF/SVG)")

def load_data():
    """Load enriched BLAST results."""
    df = pd.read_csv('/home/ubuntu/plasmid_analysis/10_BLAST_Results/blast_kpc_enriched_results.csv')
    # Deduplicate by accession (keep best hit per accession)
    df_dedup = df.sort_values('bit_score', ascending=False).drop_duplicates('accession')
    return df, df_dedup

# ============================================================
# FIGURE 1: Species Distribution Bar Chart
# ============================================================
def fig_species_distribution(df_dedup):
    print("\nGenerating Figure 1: Species Distribution...")
    
    species_counts = df_dedup['species'].value_counts()
    
    # Group rare species
    top_species = species_counts[species_counts >= 2]
    other_count = species_counts[species_counts < 2].sum()
    if other_count > 0:
        top_species['Other species'] = other_count
    
    # Color mapping
    color_map = {
        'Klebsiella pneumoniae': COLORS['kp'],
        'Escherichia coli': COLORS['ecoli'],
        'Serratia marcescens': COLORS['serratia'],
        'Proteus mirabilis': COLORS['proteus'],
        'Pseudomonas aeruginosa': COLORS['pseudo'],
        'Klebsiella aerogenes': COLORS['klebsiella_other'],
        'Citrobacter freundii': COLORS['citrobacter'],
        'Enterobacter asburiae': COLORS['enterobacter'],
        'Enterobacter hormaechei': '#D35400',
        'Citrobacter koseri': '#1ABC9C',
        'Citrobacter braakii': '#148F77',
        'Cloning vector': '#BDC3C7',
        'Other species': COLORS['other'],
    }
    
    colors = [color_map.get(s, COLORS['other']) for s in top_species.index]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.barh(range(len(top_species)), top_species.values, color=colors, edgecolor='white', linewidth=0.5)
    
    ax.set_yticks(range(len(top_species)))
    ylabels = []
    for s in top_species.index:
        if s != 'Other species' and s != 'Cloning vector':
            sp = s.replace(' ', '\\ ')
            ylabels.append(f'$\\it{{{sp}}}$')
        else:
            ylabels.append(s)
    ax.set_yticklabels(ylabels, fontsize=9)
    ax.set_xlabel('Number of BLAST Hits (≥90% identity, ≥80% coverage)', fontsize=11)
    ax.set_title('Bacterial Species Carrying Tn4401-like $\\it{bla}$$_{KPC-2}$ Cassette\n(n = 490 unique accessions)', fontsize=12, fontweight='bold')
    
    # Add count labels
    for i, (v, s) in enumerate(zip(top_species.values, top_species.index)):
        pct = v / len(df_dedup) * 100
        ax.text(v + 2, i, f'{v} ({pct:.1f}%)', va='center', fontsize=8, color='#333333')
    
    ax.set_xlim(0, max(top_species.values) * 1.2)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Species_Distribution')
    plt.close()

# ============================================================
# FIGURE 2: Geographic Distribution
# ============================================================
def fig_geographic_distribution(df_dedup):
    print("\nGenerating Figure 2: Geographic Distribution...")
    
    geo = df_dedup['geographic_origin'].value_counts()
    geo_specified = geo[geo.index != 'Not specified']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [2, 1]})
    
    # Left: Bar chart
    color_map = {
        'China': '#E74C3C',
        'Taiwan': '#3498DB',
        'Sudan': '#2ECC71',
        'Japan': '#9B59B6',
        'USA': '#F1C40F',
        'India': '#E67E22',
        'Egypt': '#1ABC9C',
        'United Kingdom': '#34495E',
        'Singapore': '#E91E63',
        'Indonesia': '#00BCD4',
        'South Korea': '#FF5722',
    }
    
    colors = [color_map.get(c, COLORS['other_country']) for c in geo_specified.index]
    
    bars = ax1.barh(range(len(geo_specified)), geo_specified.values, color=colors, edgecolor='white', linewidth=0.5)
    ax1.set_yticks(range(len(geo_specified)))
    ax1.set_yticklabels(geo_specified.index, fontsize=10)
    ax1.set_xlabel('Number of Accessions', fontsize=11)
    ax1.set_title('Geographic Distribution of Tn4401-like\n$\\it{bla}$$_{KPC-2}$ Cassette (n = 467/490)', fontsize=12, fontweight='bold')
    
    for i, v in enumerate(geo_specified.values):
        pct = v / geo_specified.sum() * 100
        ax1.text(v + 2, i, f'{v} ({pct:.1f}%)', va='center', fontsize=9)
    
    ax1.set_xlim(0, max(geo_specified.values) * 1.2)
    ax1.invert_yaxis()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Right: Pie chart for continent-level grouping
    continent_map = {
        'China': 'East Asia',
        'Taiwan': 'East Asia',
        'Japan': 'East Asia',
        'South Korea': 'East Asia',
        'Singapore': 'Southeast Asia',
        'Indonesia': 'Southeast Asia',
        'India': 'South Asia',
        'Sudan': 'Africa',
        'Egypt': 'Africa',
        'USA': 'Americas',
        'United Kingdom': 'Europe',
    }
    
    df_geo = df_dedup[df_dedup['geographic_origin'] != 'Not specified'].copy()
    df_geo['continent'] = df_geo['geographic_origin'].map(continent_map).fillna('Other')
    continent_counts = df_geo['continent'].value_counts()
    
    continent_colors = {
        'East Asia': '#E74C3C',
        'Africa': '#2ECC71',
        'Americas': '#F1C40F',
        'Europe': '#34495E',
        'South Asia': '#E67E22',
        'Southeast Asia': '#9B59B6',
        'Other': '#95A5A6',
    }
    
    colors_pie = [continent_colors.get(c, '#95A5A6') for c in continent_counts.index]
    
    wedges, texts, autotexts = ax2.pie(
        continent_counts.values, 
        labels=continent_counts.index,
        colors=colors_pie,
        autopct='%1.1f%%',
        startangle=90,
        textprops={'fontsize': 9},
        pctdistance=0.75,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1}
    )
    for autotext in autotexts:
        autotext.set_fontsize(8)
    ax2.set_title('Continental Distribution', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Geographic_Distribution')
    plt.close()

# ============================================================
# FIGURE 3: Source Type (Plasmid vs Chromosome)
# ============================================================
def fig_source_type(df_dedup):
    print("\nGenerating Figure 3: Source Type Distribution...")
    
    source_counts = df_dedup['source_type'].value_counts()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Pie chart
    colors_src = ['#C0392B', '#2980B9', '#3498DB']
    wedges, texts, autotexts = ax1.pie(
        source_counts.values,
        labels=source_counts.index,
        colors=colors_src[:len(source_counts)],
        autopct='%1.1f%%',
        startangle=90,
        textprops={'fontsize': 10},
        pctdistance=0.75,
        wedgeprops={'edgecolor': 'white', 'linewidth': 2}
    )
    for autotext in autotexts:
        autotext.set_fontsize(10)
        autotext.set_fontweight('bold')
    ax1.set_title('Genomic Location of Tn4401-like\n$\\it{bla}$$_{KPC-2}$ Cassette', fontsize=12, fontweight='bold')
    
    # Bar chart: species breakdown by source type
    species_source = pd.crosstab(df_dedup['species'], df_dedup['source_type'])
    top_sp = df_dedup['species'].value_counts().head(8).index
    species_source_top = species_source.loc[species_source.index.isin(top_sp)]
    species_source_top = species_source_top.loc[top_sp]
    
    species_source_top.plot(kind='barh', stacked=True, ax=ax2, 
                            color=['#2980B9', '#3498DB', '#C0392B'],
                            edgecolor='white', linewidth=0.5)
    ax2.set_xlabel('Number of Accessions', fontsize=11)
    ax2.set_title('Source Type by Species', fontsize=12, fontweight='bold')
    ylabels2 = []
    for s in species_source_top.index:
        sp = s.replace(' ', '\\ ')
        ylabels2.append(f'$\\it{{{sp}}}$')
    ax2.set_yticklabels(ylabels2, fontsize=9)
    ax2.legend(title='Source', fontsize=8, title_fontsize=9)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.invert_yaxis()
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Source_Type_Distribution')
    plt.close()

# ============================================================
# FIGURE 4: Sequence Identity Distribution
# ============================================================
def fig_identity_distribution(df_dedup):
    print("\nGenerating Figure 4: Sequence Identity Distribution...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    ax1.hist(df_dedup['percent_identity'], bins=50, color='#C0392B', edgecolor='white', 
             linewidth=0.5, alpha=0.85)
    ax1.axvline(x=99, color='#2C3E50', linestyle='--', linewidth=1, alpha=0.7, label='99% identity')
    ax1.set_xlabel('Percent Identity (%)', fontsize=11)
    ax1.set_ylabel('Number of Accessions', fontsize=11)
    ax1.set_title('Sequence Identity Distribution\nof Tn4401-like $\\it{bla}$$_{KPC-2}$ Cassette Hits', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Identity by species
    top_species = df_dedup['species'].value_counts().head(8).index
    data_for_box = [df_dedup[df_dedup['species'] == sp]['percent_identity'].values for sp in top_species]
    
    bp = ax2.boxplot(data_for_box, vert=True, patch_artist=True, widths=0.6)
    
    species_colors = ['#C0392B', '#2980B9', '#27AE60', '#F39C12', '#8E44AD', '#E74C3C', '#16A085', '#D35400']
    for patch, color in zip(bp['boxes'], species_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax2.set_xticklabels([s.split()[-1] for s in top_species], fontsize=8, rotation=45, ha='right', 
                         fontstyle='italic')
    ax2.set_ylabel('Percent Identity (%)', fontsize=11)
    ax2.set_title('Identity by Species', fontsize=12, fontweight='bold')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Identity_Distribution')
    plt.close()

# ============================================================
# FIGURE 5: Comprehensive Summary Panel
# ============================================================
def fig_summary_panel(df_dedup):
    print("\nGenerating Figure 5: Comprehensive Summary Panel...")
    
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)
    
    # Panel A: Species (top-left)
    ax_a = fig.add_subplot(gs[0, 0])
    species_counts = df_dedup['species'].value_counts().head(10)
    colors_sp = ['#C0392B', '#2980B9', '#27AE60', '#F39C12', '#8E44AD', '#E74C3C', '#16A085', '#D35400', '#BDC3C7', '#95A5A6']
    ax_a.barh(range(len(species_counts)), species_counts.values, color=colors_sp[:len(species_counts)], edgecolor='white')
    ax_a.set_yticks(range(len(species_counts)))
    ax_a.set_yticklabels([s.split()[-1] for s in species_counts.index], fontsize=8, fontstyle='italic')
    ax_a.set_xlabel('Count', fontsize=9)
    ax_a.set_title('A. Host Species', fontsize=11, fontweight='bold', loc='left')
    ax_a.invert_yaxis()
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)
    
    # Panel B: Geography (top-center)
    ax_b = fig.add_subplot(gs[0, 1])
    geo = df_dedup[df_dedup['geographic_origin'] != 'Not specified']['geographic_origin'].value_counts()
    geo_colors = ['#E74C3C', '#3498DB', '#2ECC71', '#9B59B6', '#F1C40F', '#E67E22', '#1ABC9C', '#34495E', '#E91E63', '#00BCD4', '#FF5722']
    ax_b.barh(range(len(geo)), geo.values, color=geo_colors[:len(geo)], edgecolor='white')
    ax_b.set_yticks(range(len(geo)))
    ax_b.set_yticklabels(geo.index, fontsize=8)
    ax_b.set_xlabel('Count', fontsize=9)
    ax_b.set_title('B. Geographic Origin', fontsize=11, fontweight='bold', loc='left')
    ax_b.invert_yaxis()
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)
    
    # Panel C: Source type pie (top-right)
    ax_c = fig.add_subplot(gs[0, 2])
    source_counts = df_dedup['source_type'].value_counts()
    wedges, texts, autotexts = ax_c.pie(
        source_counts.values, labels=source_counts.index,
        colors=['#C0392B', '#2980B9', '#3498DB'],
        autopct='%1.1f%%', startangle=90, textprops={'fontsize': 9},
        wedgeprops={'edgecolor': 'white', 'linewidth': 1.5}
    )
    ax_c.set_title('C. Genomic Location', fontsize=11, fontweight='bold', loc='left')
    
    # Panel D: Identity histogram (bottom-left)
    ax_d = fig.add_subplot(gs[1, 0])
    ax_d.hist(df_dedup['percent_identity'], bins=40, color='#C0392B', edgecolor='white', alpha=0.85)
    ax_d.set_xlabel('Identity (%)', fontsize=9)
    ax_d.set_ylabel('Count', fontsize=9)
    ax_d.set_title('D. Sequence Identity', fontsize=11, fontweight='bold', loc='left')
    ax_d.spines['top'].set_visible(False)
    ax_d.spines['right'].set_visible(False)
    
    # Panel E: Host organism (bottom-center)
    ax_e = fig.add_subplot(gs[1, 1])
    hosts = df_dedup['host_organism'].value_counts()
    hosts_clean = hosts[hosts.index != 'Not specified'].head(8)
    if len(hosts_clean) > 0:
        ax_e.barh(range(len(hosts_clean)), hosts_clean.values, color='#2980B9', edgecolor='white')
        ax_e.set_yticks(range(len(hosts_clean)))
        ax_e.set_yticklabels(hosts_clean.index, fontsize=8)
        ax_e.invert_yaxis()
    ax_e.set_xlabel('Count', fontsize=9)
    ax_e.set_title('E. Host Organism', fontsize=11, fontweight='bold', loc='left')
    ax_e.spines['top'].set_visible(False)
    ax_e.spines['right'].set_visible(False)
    
    # Panel F: Species by geography heatmap (bottom-right)
    ax_f = fig.add_subplot(gs[1, 2])
    top_sp = df_dedup['species'].value_counts().head(6).index
    top_geo = df_dedup[df_dedup['geographic_origin'] != 'Not specified']['geographic_origin'].value_counts().head(6).index
    cross = pd.crosstab(
        df_dedup[df_dedup['species'].isin(top_sp) & df_dedup['geographic_origin'].isin(top_geo)]['species'],
        df_dedup[df_dedup['species'].isin(top_sp) & df_dedup['geographic_origin'].isin(top_geo)]['geographic_origin']
    )
    if len(cross) > 0:
        sns.heatmap(cross, annot=True, fmt='d', cmap='YlOrRd', ax=ax_f, 
                    linewidths=0.5, linecolor='white', cbar_kws={'shrink': 0.8})
        ax_f.set_yticklabels([s.split()[-1] for s in cross.index], fontsize=8, fontstyle='italic', rotation=0)
        ax_f.set_xticklabels(cross.columns, fontsize=8, rotation=45, ha='right')
    ax_f.set_title('F. Species × Geography', fontsize=11, fontweight='bold', loc='left')
    
    fig.suptitle('Global Dissemination of Tn4401-like $\\it{bla}$$_{KPC-2}$ Cassette\n(BLASTn analysis: 492 hits, ≥90% identity, ≥80% coverage)', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    save_fig(fig, 'Fig_Summary_Panel_6x')
    plt.close()

# ============================================================
# FIGURE 6: Plasmid Names Word Cloud / Top Plasmids
# ============================================================
def fig_plasmid_names(df_dedup):
    print("\nGenerating Figure 6: Top Plasmid Names...")
    
    plasmid_df = df_dedup[df_dedup['source_type'] == 'Plasmid'].copy()
    
    # Filter out generic unnamed plasmids for a separate count
    named = plasmid_df[~plasmid_df['plasmid_name'].str.contains('unnamed|Not specified', case=False, na=False)]
    unnamed = plasmid_df[plasmid_df['plasmid_name'].str.contains('unnamed|Not specified', case=False, na=False)]
    
    named_counts = named['plasmid_name'].value_counts().head(30)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8), gridspec_kw={'width_ratios': [2, 1]})
    
    # Left: Top named plasmids
    colors_bar = plt.cm.YlOrRd(np.linspace(0.3, 0.9, len(named_counts)))
    ax1.barh(range(len(named_counts)), named_counts.values, color=colors_bar, edgecolor='white', linewidth=0.5)
    ax1.set_yticks(range(len(named_counts)))
    ax1.set_yticklabels(named_counts.index, fontsize=7, fontfamily='monospace')
    ax1.set_xlabel('Number of Accessions', fontsize=10)
    ax1.set_title(f'Top 30 Named Plasmids Carrying\nTn4401-like $\\it{{bla}}$$_{{KPC-2}}$ Cassette', fontsize=12, fontweight='bold')
    ax1.invert_yaxis()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    for i, v in enumerate(named_counts.values):
        ax1.text(v + 0.1, i, str(v), va='center', fontsize=7)
    
    # Right: Named vs unnamed summary
    summary_data = {
        'Named plasmids': len(named),
        'Unnamed plasmids': len(unnamed),
        'Chromosomal': len(df_dedup[df_dedup['source_type'].str.contains('Chromosome')]),
    }
    
    colors_pie = ['#C0392B', '#E74C3C', '#2980B9']
    wedges, texts, autotexts = ax2.pie(
        summary_data.values(), labels=summary_data.keys(),
        colors=colors_pie, autopct='%1.1f%%', startangle=90,
        textprops={'fontsize': 9},
        wedgeprops={'edgecolor': 'white', 'linewidth': 1.5}
    )
    ax2.set_title('Plasmid Naming Status', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Plasmid_Names')
    plt.close()

# ============================================================
# FIGURE 7: Collection Date Timeline
# ============================================================
def fig_timeline(df_dedup):
    print("\nGenerating Figure 7: Collection Date Timeline...")
    
    dates = df_dedup[df_dedup['collection_date'] != 'Not specified']['collection_date'].copy()
    
    # Extract year
    years = []
    for d in dates:
        d = str(d)
        if len(d) >= 4:
            try:
                year = int(d[:4])
                if 2000 <= year <= 2026:
                    years.append(year)
            except ValueError:
                pass
    
    if len(years) < 5:
        print("  Insufficient date data for timeline, skipping...")
        return
    
    year_counts = pd.Series(years).value_counts().sort_index()
    
    fig, ax = plt.subplots(figsize=(12, 5))
    
    bars = ax.bar(year_counts.index, year_counts.values, color='#C0392B', edgecolor='white', linewidth=0.5, alpha=0.85)
    
    # Add cumulative line
    ax2 = ax.twinx()
    cumulative = year_counts.cumsum()
    ax2.plot(cumulative.index, cumulative.values, color='#2C3E50', linewidth=2, marker='o', markersize=4)
    ax2.set_ylabel('Cumulative Count', fontsize=11, color='#2C3E50')
    ax2.tick_params(axis='y', labelcolor='#2C3E50')
    
    ax.set_xlabel('Year of Collection', fontsize=11)
    ax.set_ylabel('Number of Isolates', fontsize=11)
    ax.set_title('Temporal Distribution of Tn4401-like $\\it{bla}$$_{KPC-2}$ Cassette Detections', 
                 fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    
    for bar, v in zip(bars, year_counts.values):
        if v >= 3:
            ax.text(bar.get_x() + bar.get_width()/2, v + 0.5, str(v), ha='center', fontsize=7)
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Collection_Timeline')
    plt.close()

# ============================================================
# FIGURE 8: Species × Country Heatmap (detailed)
# ============================================================
def fig_species_country_heatmap(df_dedup):
    print("\nGenerating Figure 8: Species × Country Heatmap...")
    
    df_geo = df_dedup[df_dedup['geographic_origin'] != 'Not specified'].copy()
    
    top_sp = df_geo['species'].value_counts().head(10).index
    top_geo = df_geo['geographic_origin'].value_counts().head(8).index
    
    cross = pd.crosstab(
        df_geo[df_geo['species'].isin(top_sp)]['species'],
        df_geo[df_geo['geographic_origin'].isin(top_geo)]['geographic_origin']
    )
    
    # Reorder
    cross = cross.loc[[s for s in top_sp if s in cross.index]]
    cross = cross[[c for c in top_geo if c in cross.columns]]
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    sns.heatmap(cross, annot=True, fmt='d', cmap='YlOrRd', ax=ax,
                linewidths=1, linecolor='white', cbar_kws={'label': 'Number of Accessions', 'shrink': 0.8},
                square=False)
    
    ylabels_hm = []
    for s in cross.index:
        sp = s.replace(' ', '\\ ')
        ylabels_hm.append(f'$\\it{{{sp}}}$')
    ax.set_yticklabels(ylabels_hm, fontsize=9, rotation=0)
    ax.set_xticklabels(cross.columns, fontsize=9, rotation=45, ha='right')
    ax.set_title('Species × Geographic Origin Distribution\nof Tn4401-like $\\it{bla}$$_{KPC-2}$ Cassette', 
                 fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Species_Country_Heatmap')
    plt.close()

# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print("Generating Publication-Quality Visualizations")
    print("Global Dissemination of Tn4401-like blaKPC-2 Cassette")
    print("=" * 70)
    
    df, df_dedup = load_data()
    print(f"\nLoaded {len(df)} total records, {len(df_dedup)} unique accessions")
    
    fig_species_distribution(df_dedup)
    fig_geographic_distribution(df_dedup)
    fig_source_type(df_dedup)
    fig_identity_distribution(df_dedup)
    fig_summary_panel(df_dedup)
    fig_plasmid_names(df_dedup)
    fig_timeline(df_dedup)
    fig_species_country_heatmap(df_dedup)
    
    print("\n" + "=" * 70)
    print(f"All figures saved to: {OUTDIR}/")
    print(f"Files generated: {len(os.listdir(OUTDIR))}")
    print("=" * 70)

if __name__ == '__main__':
    main()
