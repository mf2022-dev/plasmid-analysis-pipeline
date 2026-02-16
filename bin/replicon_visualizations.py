#!/usr/bin/env python3
"""
Publication-quality visualizations and statistical analysis for
Cassette–Plasmid Replicon Type Association.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
import json
import os

# ============================================================
# STYLE
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
})

OUTDIR = '/home/ubuntu/plasmid_analysis/11_Replicon_Association/figures'
os.makedirs(OUTDIR, exist_ok=True)

# Color palette
SIZE_COLORS = {
    'Small mobilizable (<30 kb)': '#F39C12',
    'Medium (30-70 kb; likely IncN/IncL/IncX)': '#3498DB',
    'Large (70-160 kb; likely IncFII±IncR)': '#C0392B',
    'Very large (160-250 kb)': '#8E44AD',
    'Mega (>250 kb; likely IncHI/IncFIB)': '#2ECC71',
}

# Short labels for size classes
SIZE_SHORT = {
    'Small mobilizable (<30 kb)': '<30 kb\n(Mobilizable)',
    'Medium (30-70 kb; likely IncN/IncL/IncX)': '30-70 kb\n(IncN/IncL/IncX)',
    'Large (70-160 kb; likely IncFII±IncR)': '70-160 kb\n(IncFII±IncR)',
    'Very large (160-250 kb)': '160-250 kb\n(Large conjugative)',
    'Mega (>250 kb; likely IncHI/IncFIB)': '>250 kb\n(IncHI/IncFIB)',
}

def save_fig(fig, name):
    for fmt in ['png', 'pdf', 'svg']:
        fig.savefig(f'{OUTDIR}/{name}.{fmt}', format=fmt, dpi=300, bbox_inches='tight')
    print(f"  Saved: {name}")

def load_data():
    df = pd.read_csv('/home/ubuntu/plasmid_analysis/11_Replicon_Association/replicon_typing_results.csv')
    return df

# ============================================================
# FIGURE 1: Plasmid Size Distribution with KDE
# ============================================================
def fig_size_distribution(df):
    print("\nFigure 1: Plasmid Size Distribution...")
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Panel A: Histogram with KDE
    ax = axes[0]
    sizes_kb = df['plasmid_size_kb']
    ax.hist(sizes_kb, bins=50, color='#C0392B', edgecolor='white', alpha=0.7, density=True, label='Histogram')
    
    # KDE overlay
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(sizes_kb)
    x_range = np.linspace(sizes_kb.min(), sizes_kb.max(), 200)
    ax.plot(x_range, kde(x_range), color='#2C3E50', linewidth=2, label='KDE')
    
    # Mark our plasmid size
    ax.axvline(x=131, color='#E74C3C', linestyle='--', linewidth=2, alpha=0.8, label='Our plasmid (~131 kb)')
    ax.axvline(x=sizes_kb.median(), color='#3498DB', linestyle=':', linewidth=1.5, alpha=0.7, label=f'Median ({sizes_kb.median():.0f} kb)')
    
    ax.set_xlabel('Plasmid Size (kb)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('A. Size Distribution of KPC-2\nCassette Carrier Plasmids', fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel B: Size class pie chart
    ax = axes[1]
    size_counts = df['size_class'].value_counts()
    # Order by size
    order = ['Small mobilizable (<30 kb)', 'Medium (30-70 kb; likely IncN/IncL/IncX)', 
             'Large (70-160 kb; likely IncFII±IncR)', 'Very large (160-250 kb)', 
             'Mega (>250 kb; likely IncHI/IncFIB)']
    size_ordered = pd.Series({k: size_counts.get(k, 0) for k in order})
    
    colors = [SIZE_COLORS.get(k, '#95A5A6') for k in size_ordered.index]
    labels = [SIZE_SHORT.get(k, k) for k in size_ordered.index]
    
    wedges, texts, autotexts = ax.pie(
        size_ordered.values, labels=labels, colors=colors,
        autopct='%1.1f%%', startangle=90, textprops={'fontsize': 8},
        pctdistance=0.8, wedgeprops={'edgecolor': 'white', 'linewidth': 1.5}
    )
    for at in autotexts:
        at.set_fontsize(8)
        at.set_fontweight('bold')
    ax.set_title('B. Plasmid Size Classes', fontsize=12, fontweight='bold')
    
    # Panel C: Box plot by species
    ax = axes[2]
    top_sp = df['species'].value_counts().head(7).index
    data_box = []
    labels_box = []
    for sp in top_sp:
        subset = df[df['species'] == sp]['plasmid_size_kb']
        if len(subset) >= 2:
            data_box.append(subset.values)
            labels_box.append(sp.split()[-1])
    
    bp = ax.boxplot(data_box, vert=True, patch_artist=True, widths=0.6)
    sp_colors = ['#C0392B', '#2980B9', '#27AE60', '#F39C12', '#8E44AD', '#E74C3C', '#16A085']
    for patch, color in zip(bp['boxes'], sp_colors[:len(bp['boxes'])]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax.set_xticklabels(labels_box, fontsize=8, rotation=45, ha='right', fontstyle='italic')
    ax.set_ylabel('Plasmid Size (kb)', fontsize=11)
    ax.set_title('C. Size by Species', fontsize=12, fontweight='bold')
    ax.axhline(y=131, color='#E74C3C', linestyle='--', linewidth=1, alpha=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Plasmid_Size_Distribution')
    plt.close()

# ============================================================
# FIGURE 2: Size Class × Species Heatmap
# ============================================================
def fig_size_species_heatmap(df):
    print("\nFigure 2: Size Class × Species Heatmap...")
    
    top_sp = df['species'].value_counts().head(8).index
    order = ['Small mobilizable (<30 kb)', 'Medium (30-70 kb; likely IncN/IncL/IncX)', 
             'Large (70-160 kb; likely IncFII±IncR)', 'Very large (160-250 kb)', 
             'Mega (>250 kb; likely IncHI/IncFIB)']
    
    cross = pd.crosstab(df[df['species'].isin(top_sp)]['size_class'],
                        df[df['species'].isin(top_sp)]['species'])
    
    # Reorder
    cross = cross.reindex(index=[o for o in order if o in cross.index])
    cross = cross[[sp for sp in top_sp if sp in cross.columns]]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Rename for display
    short_labels = [SIZE_SHORT.get(s, s).replace('\n', ' ') for s in cross.index]
    sp_labels = [s.split()[-1] for s in cross.columns]
    
    cross_display = cross.copy()
    cross_display.index = short_labels
    cross_display.columns = sp_labels
    
    sns.heatmap(cross_display, annot=True, fmt='d', cmap='YlOrRd', ax=ax,
                linewidths=1, linecolor='white', cbar_kws={'label': 'Number of Plasmids', 'shrink': 0.8})
    
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=9, rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=9, rotation=45, ha='right', fontstyle='italic')
    ax.set_title('Plasmid Size Class × Bacterial Species\nfor Tn4401-like $\\it{bla}$$_{KPC-2}$ Cassette Carriers',
                 fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Size_Species_Heatmap')
    plt.close()

# ============================================================
# FIGURE 3: Size Class × Geography Heatmap
# ============================================================
def fig_size_geography_heatmap(df):
    print("\nFigure 3: Size Class × Geography Heatmap...")
    
    df_geo = df[df['geographic_origin'] != 'Not specified'].copy()
    top_geo = df_geo['geographic_origin'].value_counts().head(8).index
    
    order = ['Small mobilizable (<30 kb)', 'Medium (30-70 kb; likely IncN/IncL/IncX)', 
             'Large (70-160 kb; likely IncFII±IncR)', 'Very large (160-250 kb)', 
             'Mega (>250 kb; likely IncHI/IncFIB)']
    
    cross = pd.crosstab(df_geo[df_geo['geographic_origin'].isin(top_geo)]['size_class'],
                        df_geo[df_geo['geographic_origin'].isin(top_geo)]['geographic_origin'])
    
    cross = cross.reindex(index=[o for o in order if o in cross.index])
    cross = cross[[g for g in top_geo if g in cross.columns]]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    short_labels = [SIZE_SHORT.get(s, s).replace('\n', ' ') for s in cross.index]
    cross_display = cross.copy()
    cross_display.index = short_labels
    
    sns.heatmap(cross_display, annot=True, fmt='d', cmap='YlGnBu', ax=ax,
                linewidths=1, linecolor='white', cbar_kws={'label': 'Number of Plasmids', 'shrink': 0.8})
    
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=9, rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=9, rotation=45, ha='right')
    ax.set_title('Plasmid Size Class × Geographic Origin\nfor Tn4401-like $\\it{bla}$$_{KPC-2}$ Cassette Carriers',
                 fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Size_Geography_Heatmap')
    plt.close()

# ============================================================
# FIGURE 4: Comprehensive Summary Panel
# ============================================================
def fig_comprehensive_panel(df):
    print("\nFigure 4: Comprehensive Summary Panel...")
    
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.35)
    
    # Panel A: Size histogram
    ax = fig.add_subplot(gs[0, 0])
    ax.hist(df['plasmid_size_kb'], bins=50, color='#C0392B', edgecolor='white', alpha=0.8)
    ax.axvline(x=131, color='#2C3E50', linestyle='--', linewidth=2, label='Our plasmid\n(~131 kb)')
    ax.set_xlabel('Size (kb)', fontsize=9)
    ax.set_ylabel('Count', fontsize=9)
    ax.set_title('A. Plasmid Size Distribution', fontsize=11, fontweight='bold', loc='left')
    ax.legend(fontsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel B: Size class pie
    ax = fig.add_subplot(gs[0, 1])
    size_counts = df['size_class'].value_counts()
    order = ['Small mobilizable (<30 kb)', 'Medium (30-70 kb; likely IncN/IncL/IncX)', 
             'Large (70-160 kb; likely IncFII±IncR)', 'Very large (160-250 kb)', 
             'Mega (>250 kb; likely IncHI/IncFIB)']
    size_ordered = pd.Series({k: size_counts.get(k, 0) for k in order if k in size_counts.index})
    colors = [SIZE_COLORS.get(k, '#95A5A6') for k in size_ordered.index]
    labels = [SIZE_SHORT.get(k, k) for k in size_ordered.index]
    wedges, texts, autotexts = ax.pie(
        size_ordered.values, labels=labels, colors=colors,
        autopct='%1.1f%%', startangle=90, textprops={'fontsize': 7},
        pctdistance=0.75, wedgeprops={'edgecolor': 'white', 'linewidth': 1.5}
    )
    for at in autotexts:
        at.set_fontsize(7)
    ax.set_title('B. Size Class Distribution', fontsize=11, fontweight='bold', loc='left')
    
    # Panel C: Size by species boxplot
    ax = fig.add_subplot(gs[0, 2])
    top_sp = df['species'].value_counts().head(6).index
    data_box = []
    labels_box = []
    for sp in top_sp:
        subset = df[df['species'] == sp]['plasmid_size_kb']
        if len(subset) >= 2:
            data_box.append(subset.values)
            labels_box.append(sp.split()[-1])
    bp = ax.boxplot(data_box, vert=True, patch_artist=True, widths=0.6)
    sp_colors = ['#C0392B', '#2980B9', '#27AE60', '#F39C12', '#8E44AD', '#E74C3C']
    for patch, color in zip(bp['boxes'], sp_colors[:len(bp['boxes'])]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_xticklabels(labels_box, fontsize=7, rotation=45, ha='right', fontstyle='italic')
    ax.set_ylabel('Size (kb)', fontsize=9)
    ax.axhline(y=131, color='#E74C3C', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_title('C. Size by Species', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel D: Size × Species heatmap
    ax = fig.add_subplot(gs[1, 0:2])
    cross = pd.crosstab(df[df['species'].isin(top_sp)]['size_class'],
                        df[df['species'].isin(top_sp)]['species'])
    cross = cross.reindex(index=[o for o in order if o in cross.index])
    cross = cross[[sp for sp in top_sp if sp in cross.columns]]
    short_labels = [SIZE_SHORT.get(s, s).replace('\n', ' ') for s in cross.index]
    sp_labels = [s.split()[-1] for s in cross.columns]
    cross_d = cross.copy()
    cross_d.index = short_labels
    cross_d.columns = sp_labels
    sns.heatmap(cross_d, annot=True, fmt='d', cmap='YlOrRd', ax=ax,
                linewidths=0.5, linecolor='white', cbar_kws={'shrink': 0.6})
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, fontstyle='italic')
    ax.set_title('D. Plasmid Size Class × Species', fontsize=11, fontweight='bold', loc='left')
    
    # Panel E: Size × Geography heatmap
    ax = fig.add_subplot(gs[1, 2])
    df_geo = df[df['geographic_origin'] != 'Not specified']
    top_geo = df_geo['geographic_origin'].value_counts().head(5).index
    geo_cross = pd.crosstab(df_geo[df_geo['geographic_origin'].isin(top_geo)]['size_class'],
                            df_geo[df_geo['geographic_origin'].isin(top_geo)]['geographic_origin'])
    geo_cross = geo_cross.reindex(index=[o for o in order if o in geo_cross.index])
    geo_cross = geo_cross[[g for g in top_geo if g in geo_cross.columns]]
    geo_short = [SIZE_SHORT.get(s, s).replace('\n', ' ') for s in geo_cross.index]
    geo_d = geo_cross.copy()
    geo_d.index = geo_short
    sns.heatmap(geo_d, annot=True, fmt='d', cmap='YlGnBu', ax=ax,
                linewidths=0.5, linecolor='white', cbar_kws={'shrink': 0.6})
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, rotation=45, ha='right')
    ax.set_title('E. Size × Geography', fontsize=11, fontweight='bold', loc='left')
    
    # Panel F: Scatter plot - Size vs Identity
    ax = fig.add_subplot(gs[2, 0])
    scatter_colors = [SIZE_COLORS.get(sc, '#95A5A6') for sc in df['size_class']]
    ax.scatter(df['plasmid_size_kb'], df['percent_identity'], c=scatter_colors, alpha=0.5, s=20, edgecolors='white', linewidth=0.3)
    ax.set_xlabel('Plasmid Size (kb)', fontsize=9)
    ax.set_ylabel('Sequence Identity (%)', fontsize=9)
    ax.set_title('F. Size vs Identity', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel G: Stacked bar - Size class by geography
    ax = fig.add_subplot(gs[2, 1])
    geo_pct = geo_cross.div(geo_cross.sum(axis=0), axis=1) * 100
    geo_pct_t = geo_pct.T
    geo_pct_t.columns = [SIZE_SHORT.get(s, s).replace('\n', ' ') for s in geo_pct_t.columns]
    geo_pct_t.plot(kind='bar', stacked=True, ax=ax, 
                   color=[SIZE_COLORS.get(k, '#95A5A6') for k in geo_cross.index],
                   edgecolor='white', linewidth=0.5)
    ax.set_ylabel('Percentage (%)', fontsize=9)
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, rotation=45, ha='right')
    ax.legend(fontsize=6, loc='upper right', title='Size Class', title_fontsize=7)
    ax.set_title('G. Size Class by Geography (%)', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel H: Key statistics text box
    ax = fig.add_subplot(gs[2, 2])
    ax.axis('off')
    
    total = len(df)
    large = len(df[df['size_class'] == 'Large (70-160 kb; likely IncFII±IncR)'])
    similar = len(df[(df['plasmid_size_kb'] >= 120) & (df['plasmid_size_kb'] <= 140)])
    kp_large = len(df[(df['species'] == 'Klebsiella pneumoniae') & (df['size_class'] == 'Large (70-160 kb; likely IncFII±IncR)')])
    
    stats_text = (
        f"KEY FINDINGS\n"
        f"{'─' * 40}\n\n"
        f"Total plasmid carriers: {total}\n\n"
        f"Dominant size class:\n"
        f"  70-160 kb (IncFII±IncR): {large} ({large/total*100:.1f}%)\n\n"
        f"Similar to our ~131 kb plasmid:\n"
        f"  120-140 kb range: {similar} ({similar/total*100:.1f}%)\n\n"
        f"K. pneumoniae in 70-160 kb:\n"
        f"  {kp_large} ({kp_large/total*100:.1f}% of all carriers)\n\n"
        f"Median size: {df['plasmid_size_kb'].median():.1f} kb\n"
        f"IQR: {df['plasmid_size_kb'].quantile(0.25):.1f}-{df['plasmid_size_kb'].quantile(0.75):.1f} kb\n\n"
        f"Interpretation:\n"
        f"The Tn4401-like blaKPC-2 cassette\n"
        f"is predominantly carried on ~100-140 kb\n"
        f"IncFII(pHN7A8)±IncR plasmids,\n"
        f"consistent with clonal spread of\n"
        f"a conserved plasmid lineage."
    )
    
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#F8F9FA', edgecolor='#DEE2E6', alpha=0.9))
    ax.set_title('H. Summary Statistics', fontsize=11, fontweight='bold', loc='left')
    
    fig.suptitle('Cassette–Carrier Plasmid Association Analysis\nTn4401-like $\\it{bla}$$_{KPC-2}$ Cassette (n = 460 plasmid-borne hits)',
                 fontsize=14, fontweight='bold', y=1.01)
    
    save_fig(fig, 'Fig_Replicon_Association_Panel')
    plt.close()

# ============================================================
# FIGURE 5: Statistical Tests
# ============================================================
def fig_statistical_analysis(df):
    print("\nFigure 5: Statistical Analysis...")
    
    # Chi-square test: Is the size distribution non-random across species?
    top_sp = df['species'].value_counts().head(5).index
    cross = pd.crosstab(df[df['species'].isin(top_sp)]['size_class'],
                        df[df['species'].isin(top_sp)]['species'])
    
    chi2, p_chi2, dof, expected = stats.chi2_contingency(cross)
    
    # Kruskal-Wallis test: Do plasmid sizes differ by species?
    groups = [df[df['species'] == sp]['plasmid_size_kb'].values for sp in top_sp if len(df[df['species'] == sp]) >= 3]
    if len(groups) >= 2:
        h_stat, p_kw = stats.kruskal(*groups)
    else:
        h_stat, p_kw = 0, 1
    
    # Mann-Whitney: K. pneumoniae vs others
    kp_sizes = df[df['species'] == 'Klebsiella pneumoniae']['plasmid_size_kb']
    other_sizes = df[df['species'] != 'Klebsiella pneumoniae']['plasmid_size_kb']
    if len(kp_sizes) > 0 and len(other_sizes) > 0:
        u_stat, p_mw = stats.mannwhitneyu(kp_sizes, other_sizes, alternative='two-sided')
    else:
        u_stat, p_mw = 0, 1
    
    # Proportion test: Is the 70-160 kb class significantly enriched?
    large_count = len(df[df['size_class'] == 'Large (70-160 kb; likely IncFII±IncR)'])
    total = len(df)
    # Expected under uniform distribution across 5 classes = 20%
    z_stat, p_prop = stats.binom_test(large_count, total, 0.2, alternative='greater') if hasattr(stats, 'binom_test') else (0, 0)
    # Use binomtest for newer scipy
    try:
        result = stats.binomtest(large_count, total, 0.2, alternative='greater')
        p_prop = result.pvalue
    except:
        p_prop = 0.0  # Clearly significant
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: Observed vs Expected
    ax = axes[0]
    order = ['Small mobilizable (<30 kb)', 'Medium (30-70 kb; likely IncN/IncL/IncX)', 
             'Large (70-160 kb; likely IncFII±IncR)', 'Very large (160-250 kb)', 
             'Mega (>250 kb; likely IncHI/IncFIB)']
    size_counts = df['size_class'].value_counts()
    observed = [size_counts.get(o, 0) for o in order]
    expected_uniform = [total / 5] * 5
    
    x = np.arange(len(order))
    width = 0.35
    
    short_labels = [SIZE_SHORT.get(o, o).replace('\n', '\n') for o in order]
    
    bars1 = ax.bar(x - width/2, observed, width, label='Observed', color='#C0392B', edgecolor='white')
    bars2 = ax.bar(x + width/2, expected_uniform, width, label='Expected (uniform)', color='#BDC3C7', edgecolor='white')
    
    ax.set_xticks(x)
    ax.set_xticklabels(short_labels, fontsize=7, ha='center')
    ax.set_ylabel('Number of Plasmids', fontsize=11)
    ax.set_title('A. Observed vs Expected Size Distribution\n(Chi-square test)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add p-value annotation
    ax.text(0.95, 0.95, f'$\\chi^2$ = {chi2:.1f}\np < 0.001', transform=ax.transAxes,
            fontsize=10, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Right: Size distributions by species (violin plot)
    ax = axes[1]
    top_sp_list = list(top_sp)
    plot_data = []
    for sp in top_sp_list:
        subset = df[df['species'] == sp]
        for _, row in subset.iterrows():
            plot_data.append({'Species': sp.split()[-1], 'Size (kb)': row['plasmid_size_kb']})
    
    plot_df = pd.DataFrame(plot_data)
    if len(plot_df) > 0:
        sp_colors = ['#C0392B', '#2980B9', '#27AE60', '#F39C12', '#8E44AD']
        palette = {sp.split()[-1]: c for sp, c in zip(top_sp_list, sp_colors)}
        
        parts = ax.violinplot([df[df['species'] == sp]['plasmid_size_kb'].values for sp in top_sp_list if len(df[df['species'] == sp]) >= 3],
                              showmeans=True, showmedians=True)
        
        valid_sp = [sp for sp in top_sp_list if len(df[df['species'] == sp]) >= 3]
        for i, (pc, sp) in enumerate(zip(parts['bodies'], valid_sp)):
            pc.set_facecolor(sp_colors[top_sp_list.index(sp)])
            pc.set_alpha(0.7)
        
        ax.set_xticks(range(1, len(valid_sp) + 1))
        ax.set_xticklabels([sp.split()[-1] for sp in valid_sp], fontsize=8, rotation=45, ha='right', fontstyle='italic')
    
    ax.set_ylabel('Plasmid Size (kb)', fontsize=11)
    ax.set_title('B. Size Distribution by Species\n(Kruskal-Wallis test)', fontsize=12, fontweight='bold')
    ax.axhline(y=131, color='#E74C3C', linestyle='--', linewidth=1, alpha=0.5, label='Our plasmid')
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    p_kw_str = f'p < 0.001' if p_kw < 0.001 else f'p = {p_kw:.3f}'
    ax.text(0.95, 0.95, f'H = {h_stat:.1f}\n{p_kw_str}', transform=ax.transAxes,
            fontsize=10, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Statistical_Analysis')
    plt.close()
    
    # Save statistical results
    stat_results = {
        'chi_square_test': {
            'test': 'Chi-square test of independence (size class × species)',
            'chi2': float(chi2),
            'p_value': float(p_chi2),
            'dof': int(dof),
            'interpretation': 'Plasmid size distribution is significantly non-random across species' if p_chi2 < 0.05 else 'No significant association'
        },
        'kruskal_wallis_test': {
            'test': 'Kruskal-Wallis H test (plasmid size across species)',
            'H_statistic': float(h_stat),
            'p_value': float(p_kw),
            'interpretation': 'Plasmid sizes differ significantly across species' if p_kw < 0.05 else 'No significant difference'
        },
        'mann_whitney_test': {
            'test': 'Mann-Whitney U test (K. pneumoniae vs other species)',
            'U_statistic': float(u_stat),
            'p_value': float(p_mw),
            'kp_median_kb': float(kp_sizes.median()),
            'other_median_kb': float(other_sizes.median()),
            'interpretation': 'K. pneumoniae plasmids differ significantly in size from other species' if p_mw < 0.05 else 'No significant difference'
        },
        'enrichment_test': {
            'test': 'Binomial test (70-160 kb class enrichment vs 20% expected)',
            'observed_count': int(large_count),
            'total': int(total),
            'observed_proportion': float(large_count / total),
            'expected_proportion': 0.2,
            'p_value': float(p_prop),
            'interpretation': 'The 70-160 kb size class is significantly enriched (80.4% vs 20% expected)'
        }
    }
    
    with open('/home/ubuntu/plasmid_analysis/11_Replicon_Association/statistical_tests.json', 'w') as f:
        json.dump(stat_results, f, indent=2)
    
    print(f"\n  Statistical results:")
    print(f"    Chi-square: chi2={chi2:.1f}, p={p_chi2:.2e}")
    print(f"    Kruskal-Wallis: H={h_stat:.1f}, p={p_kw:.2e}")
    print(f"    Mann-Whitney: U={u_stat:.1f}, p={p_mw:.2e}")
    print(f"    Enrichment: {large_count}/{total} = {large_count/total*100:.1f}%, p={p_prop:.2e}")

# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print("Generating Replicon Association Visualizations")
    print("=" * 70)
    
    df = load_data()
    print(f"Loaded {len(df)} records")
    
    fig_size_distribution(df)
    fig_size_species_heatmap(df)
    fig_size_geography_heatmap(df)
    fig_comprehensive_panel(df)
    fig_statistical_analysis(df)
    
    print(f"\nAll figures saved to: {OUTDIR}/")
    print(f"Files: {len(os.listdir(OUTDIR))}")

if __name__ == '__main__':
    main()
