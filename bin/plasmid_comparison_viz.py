#!/usr/bin/env python3
"""
Publication-quality visualizations for whole-plasmid comparison:
1. Pairwise identity heatmap (SNP matrix)
2. Aligned coverage heatmap
3. Gene content comparison matrix
4. Comprehensive summary panel
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import seaborn as sns
import json
import os

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

WORKDIR = "/home/ubuntu/plasmid_analysis/12_Plasmid_Comparison"
FIGDIR = os.path.join(WORKDIR, "figures")
os.makedirs(FIGDIR, exist_ok=True)

def save_fig(fig, name):
    for fmt in ['png', 'pdf', 'svg']:
        fig.savefig(f'{FIGDIR}/{name}.{fmt}', format=fmt, dpi=300, bbox_inches='tight')
    print(f"  Saved: {name}")

# ============================================================
# FIGURE 1: Identity Heatmap
# ============================================================
def fig_identity_heatmap():
    print("\nFigure 1: Pairwise Identity Heatmap...")
    
    id_df = pd.read_csv(f'{WORKDIR}/identity_matrix.csv', index_col=0)
    
    # Rename labels for display
    rename = {}
    for col in id_df.columns:
        if col.startswith('SA_'):
            rename[col] = col.replace('SA_', '') + ' (SA)'
        else:
            rename[col] = col.replace('Ref_', '')
    
    id_display = id_df.rename(index=rename, columns=rename)
    
    # Reorder: our plasmids first, then references
    our_labels = [l for l in id_display.index if '(SA)' in l]
    ref_labels = [l for l in id_display.index if '(SA)' not in l]
    order = our_labels + ref_labels
    id_display = id_display.loc[order, order]
    
    fig, ax = plt.subplots(figsize=(14, 12))
    
    # Mask diagonal
    mask = np.eye(len(id_display), dtype=bool)
    
    sns.heatmap(id_display, annot=True, fmt='.2f', cmap='RdYlGn', ax=ax,
                mask=mask, vmin=99.0, vmax=100.0,
                linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'Average Nucleotide Identity (%)', 'shrink': 0.7},
                annot_kws={'fontsize': 7})
    
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, rotation=0)
    
    # Add separator between our and reference plasmids
    n_our = len(our_labels)
    ax.axhline(y=n_our, color='black', linewidth=2)
    ax.axvline(x=n_our, color='black', linewidth=2)
    
    # Add group labels
    ax.text(-0.5, n_our/2, 'This Study\n(Saudi Arabia)', ha='right', va='center',
            fontsize=10, fontweight='bold', color='#C0392B', rotation=90)
    ax.text(-0.5, n_our + len(ref_labels)/2, 'Global\nReferences', ha='right', va='center',
            fontsize=10, fontweight='bold', color='#2980B9', rotation=90)
    
    ax.set_title('Pairwise Average Nucleotide Identity (%)\nKPC-2 Carrying IncFII+IncR Plasmids',
                 fontsize=13, fontweight='bold', pad=15)
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Identity_Heatmap')
    plt.close()

# ============================================================
# FIGURE 2: SNP Matrix Heatmap
# ============================================================
def fig_snp_heatmap():
    print("\nFigure 2: SNP Matrix Heatmap...")
    
    snp_df = pd.read_csv(f'{WORKDIR}/snp_matrix.csv', index_col=0)
    
    rename = {}
    for col in snp_df.columns:
        if col.startswith('SA_'):
            rename[col] = col.replace('SA_', '') + ' (SA)'
        else:
            rename[col] = col.replace('Ref_', '')
    
    snp_display = snp_df.rename(index=rename, columns=rename)
    
    our_labels = [l for l in snp_display.index if '(SA)' in l]
    ref_labels = [l for l in snp_display.index if '(SA)' not in l]
    order = our_labels + ref_labels
    snp_display = snp_display.loc[order, order]
    
    fig, ax = plt.subplots(figsize=(14, 12))
    
    mask = np.eye(len(snp_display), dtype=bool)
    
    sns.heatmap(snp_display, annot=True, fmt='d', cmap='YlOrRd_r', ax=ax,
                mask=mask,
                linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'Number of SNPs', 'shrink': 0.7},
                annot_kws={'fontsize': 7})
    
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, rotation=0)
    
    n_our = len(our_labels)
    ax.axhline(y=n_our, color='black', linewidth=2)
    ax.axvline(x=n_our, color='black', linewidth=2)
    
    ax.set_title('Pairwise SNP Distances\nKPC-2 Carrying IncFII+IncR Plasmids',
                 fontsize=13, fontweight='bold', pad=15)
    
    plt.tight_layout()
    save_fig(fig, 'Fig_SNP_Heatmap')
    plt.close()

# ============================================================
# FIGURE 3: Gene Content Comparison
# ============================================================
def fig_gene_content():
    print("\nFigure 3: Gene Content Comparison...")
    
    # AMR gene presence/absence data
    isolates = ['JKPB244', 'JKPB284', 'JKPB381', 'JKPR282', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']
    genes = ['blaKPC-2', 'blaCTX-M-65', 'blaSHV-12', 'blaTEM-1', 'fosA3', 'rmtB1', 'catA2', 'merP',
             'erm(B)', 'tet(X5)', 'sul1', 'aadA12']
    
    # Presence matrix
    data = {
        'JKPB244':  [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
        'JKPB284':  [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
        'JKPB381':  [1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        'JKPR282':  [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
        'KP21802':  [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
        'KP21915':  [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
        'MDRKP121': [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
        'MDRKP122': [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
    }
    
    gene_df = pd.DataFrame(data, index=genes).T
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Custom colormap: 0=white, 1=red
    cmap = matplotlib.colors.ListedColormap(['#F5F5F5', '#C0392B'])
    
    sns.heatmap(gene_df, annot=False, cmap=cmap, ax=ax,
                linewidths=1, linecolor='white', cbar=False,
                vmin=0, vmax=1)
    
    # Add + and - symbols
    for i in range(len(gene_df)):
        for j in range(len(gene_df.columns)):
            val = gene_df.iloc[i, j]
            symbol = '+' if val == 1 else '-'
            color = 'white' if val == 1 else '#BDC3C7'
            ax.text(j + 0.5, i + 0.5, symbol, ha='center', va='center',
                    fontsize=12, fontweight='bold', color=color)
    
    ax.set_xticklabels(genes, fontsize=9, rotation=45, ha='right', fontstyle='italic')
    ax.set_yticklabels(isolates, fontsize=9, rotation=0)
    
    # Add plasmid size annotations
    sizes = {'JKPB244': '130.7 kb', 'JKPB284': '130.7 kb', 'JKPB381': '105.4 kb',
             'JKPR282': '130.6 kb*', 'KP21802': '130.7 kb', 'KP21915': '130.7 kb',
             'MDRKP121': '130.6 kb', 'MDRKP122': '130.6 kb'}
    
    for i, iso in enumerate(isolates):
        ax.text(len(genes) + 0.3, i + 0.5, sizes[iso], ha='left', va='center', fontsize=8)
    
    ax.set_title('AMR Gene Content of KPC-2 Carrying Plasmids\nacross 8 ST11 $\\it{K. pneumoniae}$ Isolates (Saudi Arabia)',
                 fontsize=12, fontweight='bold', pad=15)
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#C0392B', edgecolor='white', label='Present'),
        mpatches.Patch(facecolor='#F5F5F5', edgecolor='#BDC3C7', label='Absent'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9, framealpha=0.9)
    
    # Add footnote
    fig.text(0.5, -0.02, '* JKPR282 carries a different plasmid backbone with erm(B)/tet(X5)/sul1/aadA12',
             ha='center', fontsize=8, fontstyle='italic', color='#7F8C8D')
    
    plt.tight_layout()
    save_fig(fig, 'Fig_Gene_Content_Comparison')
    plt.close()

# ============================================================
# FIGURE 4: Comprehensive Summary Panel
# ============================================================
def fig_comprehensive_panel():
    print("\nFigure 4: Comprehensive Summary Panel...")
    
    fig = plt.figure(figsize=(20, 16))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
    
    # Load data
    id_df = pd.read_csv(f'{WORKDIR}/identity_matrix.csv', index_col=0)
    snp_df = pd.read_csv(f'{WORKDIR}/snp_matrix.csv', index_col=0)
    results_df = pd.read_csv(f'{WORKDIR}/pairwise_comparison_results.csv')
    
    rename = {}
    for col in id_df.columns:
        if col.startswith('SA_'):
            rename[col] = col.replace('SA_', '') + '\n(SA)'
        else:
            rename[col] = col.replace('Ref_', '')
    
    # Panel A: Identity heatmap (our plasmids only)
    ax = fig.add_subplot(gs[0, 0])
    our_cols = [c for c in id_df.columns if c.startswith('SA_')]
    our_id = id_df.loc[our_cols, our_cols]
    our_rename = {c: c.replace('SA_', '') for c in our_cols}
    our_id_d = our_id.rename(index=our_rename, columns=our_rename)
    mask = np.eye(len(our_id_d), dtype=bool)
    sns.heatmap(our_id_d, annot=True, fmt='.2f', cmap='RdYlGn', ax=ax,
                mask=mask, vmin=99.8, vmax=100.0,
                linewidths=0.5, linecolor='white', cbar_kws={'shrink': 0.6},
                annot_kws={'fontsize': 6})
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=7, rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
    ax.set_title('A. Intra-study ANI (%)', fontsize=11, fontweight='bold', loc='left')
    
    # Panel B: SNP heatmap (our plasmids only)
    ax = fig.add_subplot(gs[0, 1])
    our_snp = snp_df.loc[our_cols, our_cols]
    our_snp_d = our_snp.rename(index=our_rename, columns=our_rename)
    sns.heatmap(our_snp_d, annot=True, fmt='.0f', cmap='YlOrRd_r', ax=ax,
                mask=mask, linewidths=0.5, linecolor='white',
                cbar_kws={'shrink': 0.6}, annot_kws={'fontsize': 7})
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=7, rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
    ax.set_title('B. Intra-study SNPs', fontsize=11, fontweight='bold', loc='left')
    
    # Panel C: SNP distribution comparison
    ax = fig.add_subplot(gs[0, 2])
    intra = results_df[(results_df['group_1'] == 'This study') & (results_df['group_2'] == 'This study')]['total_snps']
    cross = results_df[
        ((results_df['group_1'] == 'This study') & (results_df['group_2'] == 'Reference')) |
        ((results_df['group_1'] == 'Reference') & (results_df['group_2'] == 'This study'))
    ]['total_snps']
    inter = results_df[(results_df['group_1'] == 'Reference') & (results_df['group_2'] == 'Reference')]['total_snps']
    
    bp = ax.boxplot([intra.values, cross.values, inter.values],
                    labels=['Intra-study\n(SA vs SA)', 'Study vs\nGlobal Refs', 'Inter-reference\n(Global vs Global)'],
                    patch_artist=True, widths=0.6)
    colors = ['#C0392B', '#8E44AD', '#2980B9']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_ylabel('Number of SNPs', fontsize=9)
    ax.set_title('C. SNP Distribution by Comparison Group', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel D: Identity heatmap (our vs references)
    ax = fig.add_subplot(gs[1, 0:2])
    ref_cols = [c for c in id_df.columns if c.startswith('Ref_')]
    cross_id = id_df.loc[our_cols, ref_cols]
    cross_rename_r = {c: c.replace('SA_', '') for c in our_cols}
    cross_rename_c = {c: c.replace('Ref_', '') for c in ref_cols}
    cross_id_d = cross_id.rename(index=cross_rename_r, columns=cross_rename_c)
    sns.heatmap(cross_id_d, annot=True, fmt='.2f', cmap='RdYlGn', ax=ax,
                vmin=99.0, vmax=100.0, linewidths=0.5, linecolor='white',
                cbar_kws={'shrink': 0.5, 'label': 'ANI (%)'}, annot_kws={'fontsize': 7})
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=7, rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
    ax.set_title('D. Study vs Global Reference ANI (%)', fontsize=11, fontweight='bold', loc='left')
    
    # Panel E: Gene content
    ax = fig.add_subplot(gs[1, 2])
    isolates = ['JKPB244', 'JKPB284', 'JKPB381', 'JKPR282', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']
    genes = ['blaKPC-2', 'blaCTX-M-65', 'blaSHV-12', 'blaTEM-1', 'fosA3', 'rmtB1', 'catA2']
    data = {
        'JKPB244':  [1, 1, 1, 1, 1, 1, 1],
        'JKPB284':  [1, 1, 1, 1, 1, 1, 1],
        'JKPB381':  [1, 1, 1, 0, 0, 0, 1],
        'JKPR282':  [0, 0, 0, 0, 0, 0, 0],
        'KP21802':  [1, 1, 1, 1, 1, 1, 1],
        'KP21915':  [1, 1, 1, 1, 1, 1, 1],
        'MDRKP121': [1, 1, 1, 1, 1, 1, 1],
        'MDRKP122': [1, 1, 1, 1, 1, 1, 1],
    }
    gene_df = pd.DataFrame(data, index=genes).T
    cmap = matplotlib.colors.ListedColormap(['#F5F5F5', '#C0392B'])
    sns.heatmap(gene_df, annot=False, cmap=cmap, ax=ax,
                linewidths=0.5, linecolor='white', cbar=False, vmin=0, vmax=1)
    for i in range(len(gene_df)):
        for j in range(len(gene_df.columns)):
            val = gene_df.iloc[i, j]
            symbol = '+' if val == 1 else '-'
            color = 'white' if val == 1 else '#BDC3C7'
            ax.text(j + 0.5, i + 0.5, symbol, ha='center', va='center',
                    fontsize=9, fontweight='bold', color=color)
    ax.set_xticklabels(genes, fontsize=6, rotation=45, ha='right', fontstyle='italic')
    ax.set_yticklabels(isolates, fontsize=7, rotation=0)
    ax.set_title('E. AMR Gene Content', fontsize=11, fontweight='bold', loc='left')
    
    # Panel F: Plasmid size comparison
    ax = fig.add_subplot(gs[2, 0])
    our_sizes = [130704, 130703, 105383, 130556, 130704, 130704, 130585, 130586]
    ref_sizes = [128563, 129252, 127240, 125619, 135837, 130677, 129466, 121348, 129869, 127668]
    our_labels_short = ['JKPB244', 'JKPB284', 'JKPB381', 'JKPR282', 'KP21802', 'KP21915', 'MDRKP121', 'MDRKP122']
    
    all_sizes = our_sizes + ref_sizes
    all_labels_short = our_labels_short + ['CP025458', 'CP045264', 'CP047161', 'CP073303', 'CP082020', 'CP097693', 'CP192284', 'KX236178', 'MT648513', 'MZ709016']
    all_colors = ['#C0392B'] * 8 + ['#2980B9'] * 10
    
    bars = ax.barh(range(len(all_sizes)), [s/1000 for s in all_sizes], color=all_colors, edgecolor='white', alpha=0.8)
    ax.set_yticks(range(len(all_sizes)))
    ax.set_yticklabels(all_labels_short, fontsize=7)
    ax.set_xlabel('Plasmid Size (kb)', fontsize=9)
    ax.set_title('F. Plasmid Sizes', fontsize=11, fontweight='bold', loc='left')
    ax.axvline(x=130.7, color='#E74C3C', linestyle='--', linewidth=1, alpha=0.5)
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    legend_elements = [
        mpatches.Patch(facecolor='#C0392B', alpha=0.8, label='This study (Saudi Arabia)'),
        mpatches.Patch(facecolor='#2980B9', alpha=0.8, label='Global references'),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc='lower right')
    
    # Panel G: Aligned coverage heatmap (our vs refs)
    ax = fig.add_subplot(gs[2, 1])
    aligned_df = pd.read_csv(f'{WORKDIR}/aligned_coverage_matrix.csv', index_col=0)
    cross_aligned = aligned_df.loc[our_cols, ref_cols]
    cross_aligned_d = cross_aligned.rename(index=cross_rename_r, columns=cross_rename_c)
    sns.heatmap(cross_aligned_d, annot=True, fmt='.1f', cmap='YlGnBu', ax=ax,
                vmin=30, vmax=100, linewidths=0.5, linecolor='white',
                cbar_kws={'shrink': 0.6, 'label': 'Aligned (%)'}, annot_kws={'fontsize': 6})
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=6, rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
    ax.set_title('G. Aligned Coverage (Study vs Refs)', fontsize=11, fontweight='bold', loc='left')
    
    # Panel H: Key findings text
    ax = fig.add_subplot(gs[2, 2])
    ax.axis('off')
    
    findings = (
        "KEY FINDINGS\n"
        f"{'─' * 45}\n\n"
        "INTRA-STUDY COMPARISON:\n"
        "  • 6/8 plasmids nearly identical\n"
        "    (0-14 SNPs, >99.9% ANI)\n"
        "  • JKPB381: truncated variant (105 kb)\n"
        "  • JKPR282: different backbone\n\n"
        "STUDY vs GLOBAL REFERENCES:\n"
        "  • ANI: 99.06-99.99%\n"
        "  • Mean SNPs: 34.9\n"
        "  • Highest similarity to CP097693\n"
        "    (China, 130.7 kb, same size)\n\n"
        "CORE AMR MODULE (6/8 isolates):\n"
        "  blaKPC-2 + blaCTX-M-65 + blaSHV-12\n"
        "  + blaTEM-1 + fosA3 + rmtB1 + catA2\n\n"
        "CONCLUSION:\n"
        "  The IncFII(pHN7A8)+IncR KPC-2\n"
        "  plasmid in Saudi Arabia is part of\n"
        "  a globally conserved plasmid lineage\n"
        "  circulating in ST11 K. pneumoniae,\n"
        "  with highest similarity to Chinese\n"
        "  isolates, suggesting intercontinental\n"
        "  plasmid transfer."
    )
    
    ax.text(0.05, 0.95, findings, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#F8F9FA', edgecolor='#DEE2E6', alpha=0.9))
    ax.set_title('H. Summary', fontsize=11, fontweight='bold', loc='left')
    
    fig.suptitle('Whole-Plasmid Comparison: KPC-2 Carrying IncFII+IncR Plasmids\n'
                 '8 ST11 $\\it{K. pneumoniae}$ (Saudi Arabia) vs 10 Global References',
                 fontsize=14, fontweight='bold', y=1.01)
    
    save_fig(fig, 'Fig_Plasmid_Comparison_Panel')
    plt.close()

# ============================================================
# MAIN
# ============================================================
def main():
    print("=" * 70)
    print("Generating Plasmid Comparison Visualizations")
    print("=" * 70)
    
    fig_identity_heatmap()
    fig_snp_heatmap()
    fig_gene_content()
    fig_comprehensive_panel()
    
    print(f"\nAll figures saved to: {FIGDIR}/")
    print(f"Files: {len(os.listdir(FIGDIR))}")

if __name__ == '__main__':
    main()
