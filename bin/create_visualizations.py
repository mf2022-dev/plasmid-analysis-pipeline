#!/usr/bin/env python3
"""
Create publication-ready visualizations for ST11 K. pneumoniae plasmid analysis.
Figures: AMR heatmap, plasmid architecture, carbapenemase distribution, MGE summary.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import numpy as np
import json
import os

OUTPUT_DIR = "/home/ubuntu/plasmid_analysis/07_Visualizations"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Apply style first, then set font
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'figure.dpi': 300,
})

# ============================================================
# Load data
# ============================================================
amr_df = pd.read_csv("/home/ubuntu/plasmid_analysis/04_AMR_Plasmid_Association/amr_genes_all_samples.csv")
mge_df = pd.read_csv("/home/ubuntu/plasmid_analysis/04_AMR_Plasmid_Association/mge_elements_all_samples.csv")

samples = ['JKPB244', 'JKPB284', 'JKPB381', 'JKPR282', 'KP21802']

# ============================================================
# Figure 1: Carbapenemase Distribution Bar Chart
# ============================================================
print("Creating Figure 1: Carbapenemase Distribution...")

carb_data = {
    'Isolate': samples,
    'blaKPC-2': [1, 1, 1, 1, 1],
    'blaNDM-1': [1, 1, 0, 0, 1],
    'blaOXA-232': [0, 0, 1, 0, 0],
    'blaCTX-M-65': [1, 1, 1, 1, 1],
    'blaSHV-12': [1, 1, 1, 1, 1],
    'blaTEM-1': [1, 1, 0, 1, 1],
}
carb_df = pd.DataFrame(carb_data)

fig, ax = plt.subplots(figsize=(10, 5))
genes = ['blaKPC-2', 'blaNDM-1', 'blaOXA-232', 'blaCTX-M-65', 'blaSHV-12', 'blaTEM-1']
colors = ['#d62728', '#ff7f0e', '#9467bd', '#2ca02c', '#1f77b4', '#8c564b']

x = np.arange(len(samples))
width = 0.12

for i, (gene, color) in enumerate(zip(genes, colors)):
    vals = carb_df[gene].values
    ax.bar(x + i * width - width * 2.5, vals, width, label=gene, color=color, edgecolor='white')

ax.set_xlabel('Isolate')
ax.set_ylabel('Presence (1=Yes, 0=No)')
ax.set_title('Beta-Lactamase Gene Distribution in ST11 K. pneumoniae Isolates')
ax.set_xticks(x)
ax.set_xticklabels(samples, rotation=45, ha='right')
ax.set_ylim(0, 1.3)
ax.legend(loc='upper right', fontsize=8, ncol=2)
ax.set_yticks([0, 1])
ax.set_yticklabels(['Absent', 'Present'])
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/Figure1_Carbapenemase_Distribution.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUTPUT_DIR}/Figure1_Carbapenemase_Distribution.pdf", bbox_inches='tight')
plt.close()
print("  Saved Figure 1")

# ============================================================
# Figure 2: AMR Gene Heatmap (Plasmid-borne genes)
# ============================================================
print("Creating Figure 2: AMR Heatmap...")

# Get plasmid-borne AMR genes
plasmid_amr = amr_df[amr_df['contig_type'] == 'plasmid'].copy()

# Key AMR genes to show
key_genes = [
    'blaKPC-2', 'blaNDM-1', 'blaOXA-232', 'blaCTX-M-65', 'blaSHV-12', 'blaTEM-1',
    'aph(3\')-Ia', 'aph(3\')-VI', 'armA', 'rmtB1', 'aadA12',
    'fosA3', 'catA2', 'ble',
    'sul1', 'sul2', 'dfrA5',
    'qnrS1', 'mph(A)', 'mph(E)', 'msr(E)',
    'iucA', 'iucB', 'iucC', 'iutA', 'rmpA', 'rmpC',
    'terC', 'terD', 'merA', 'merR',
]

# Build presence/absence matrix
heatmap_data = pd.DataFrame(0, index=key_genes, columns=samples)
for _, row in plasmid_amr.iterrows():
    gene = row['gene']
    sample = row['sample']
    if gene in key_genes and sample in samples:
        heatmap_data.loc[gene, sample] = 1

# Also check chromosome for some genes
chromo_amr = amr_df[amr_df['contig_type'] == 'chromosome']
for _, row in chromo_amr.iterrows():
    gene = row['gene']
    sample = row['sample']
    if gene in key_genes and sample in samples:
        if heatmap_data.loc[gene, sample] == 0:
            heatmap_data.loc[gene, sample] = 0.5  # chromosome only

# Remove genes with all zeros
heatmap_data = heatmap_data.loc[heatmap_data.sum(axis=1) > 0]

fig, ax = plt.subplots(figsize=(8, 12))
cmap = sns.color_palette(["#f0f0f0", "#ffd700", "#d62728"])
sns.heatmap(heatmap_data, annot=False, cmap=cmap, linewidths=0.5, linecolor='white',
            cbar_kws={'label': 'Location', 'ticks': [0, 0.5, 1]},
            vmin=0, vmax=1, ax=ax)
ax.set_title('Plasmid-Borne AMR Gene Distribution\nin ST11 K. pneumoniae', fontsize=13, fontweight='bold')
ax.set_xlabel('Isolate')
ax.set_ylabel('AMR Gene')

# Custom legend
legend_elements = [
    mpatches.Patch(facecolor='#f0f0f0', edgecolor='gray', label='Absent'),
    mpatches.Patch(facecolor='#ffd700', edgecolor='gray', label='Chromosome only'),
    mpatches.Patch(facecolor='#d62728', edgecolor='gray', label='Plasmid-borne'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/Figure2_AMR_Heatmap.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUTPUT_DIR}/Figure2_AMR_Heatmap.pdf", bbox_inches='tight')
plt.close()
print("  Saved Figure 2")

# ============================================================
# Figure 3: Plasmid Architecture Comparison
# ============================================================
print("Creating Figure 3: Plasmid Architecture...")

plasmid_data = {
    'JKPB244': [('Chromosome', 5454042), ('pNDM-Mar ~363kb', 363129), ('pKPC ~131kb', 130704), ('ColRNAI ~10kb', 10060), ('Col ~6kb', 5596)],
    'JKPB284': [('Chromosome', 5453998), ('pNDM-Mar ~362kb', 361821), ('pKPC ~131kb', 130703), ('ColRNAI ~10kb', 10060), ('Col ~6kb', 5596)],
    'JKPB381': [('Chromosome', 5327974), ('Large plasmid ~353kb', 353936), ('pKPC ~105kb', 105362), ('ColRNAI ~10kb', 10060), ('Col ~6kb', 5596)],
    'JKPR282': [('Chromosome', 5355006), ('Large plasmid ~353kb', 353936), ('pKPC ~131kb', 130704), ('ColRNAI ~12kb', 11970), ('Col ~6kb', 5596)],
    'KP21802': [('Chromosome', 5449982), ('pNDM-Mar ~360kb', 360000), ('pKPC ~131kb', 130704), ('ColRNAI ~12kb', 11970), ('Col ~6kb', 5596)],
}

fig, axes = plt.subplots(1, 5, figsize=(16, 6), sharey=True)

colors_plasmid = {
    'Chromosome': '#4a90d9',
    'pNDM-Mar': '#d62728',
    'pKPC': '#ff7f0e',
    'Large plasmid': '#9467bd',
    'ColRNAI': '#2ca02c',
    'Col': '#8c564b',
}

for idx, (sample, contigs) in enumerate(plasmid_data.items()):
    ax = axes[idx]
    sizes = [c[1]/1000 for c in contigs[1:]]  # Skip chromosome, convert to kb
    labels = [c[0] for c in contigs[1:]]
    
    bar_colors = []
    for label in labels:
        matched = False
        for key, color in colors_plasmid.items():
            if key in label:
                bar_colors.append(color)
                matched = True
                break
        if not matched:
            bar_colors.append('#999999')
    
    bars = ax.barh(range(len(sizes)), sizes, color=bar_colors, edgecolor='white', height=0.6)
    ax.set_yticks(range(len(sizes)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('Size (kb)')
    ax.set_title(sample, fontsize=10, fontweight='bold')
    ax.invert_yaxis()
    
    # Add size labels
    for bar, size in zip(bars, sizes):
        ax.text(bar.get_width() + 5, bar.get_y() + bar.get_height()/2, 
                f'{size:.0f} kb', va='center', fontsize=7)

fig.suptitle('Plasmid Architecture Comparison in ST11 K. pneumoniae', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/Figure3_Plasmid_Architecture.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUTPUT_DIR}/Figure3_Plasmid_Architecture.pdf", bbox_inches='tight')
plt.close()
print("  Saved Figure 3")

# ============================================================
# Figure 4: MGE Element Distribution
# ============================================================
print("Creating Figure 4: MGE Distribution...")

mge_summary = mge_df.groupby(['sample', 'contig_type']).size().unstack(fill_value=0)

fig, ax = plt.subplots(figsize=(10, 5))
mge_summary.plot(kind='bar', stacked=True, ax=ax, color=['#4a90d9', '#d62728'], edgecolor='white')
ax.set_xlabel('Isolate')
ax.set_ylabel('Number of MGE Elements')
ax.set_title('Mobile Genetic Element Distribution\nin ST11 K. pneumoniae', fontweight='bold')
ax.legend(title='Location', labels=['Chromosome', 'Plasmid'])
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/Figure4_MGE_Distribution.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUTPUT_DIR}/Figure4_MGE_Distribution.pdf", bbox_inches='tight')
plt.close()
print("  Saved Figure 4")

# ============================================================
# Figure 5: AMR Class Distribution (Pie Chart)
# ============================================================
print("Creating Figure 5: AMR Class Distribution...")

amr_classes = amr_df[amr_df['contig_type'] == 'plasmid']['amr_class'].value_counts()

fig, ax = plt.subplots(figsize=(8, 8))
colors_pie = plt.cm.Set3(np.linspace(0, 1, len(amr_classes)))
wedges, texts, autotexts = ax.pie(amr_classes.values, labels=amr_classes.index, 
                                    autopct='%1.1f%%', colors=colors_pie,
                                    pctdistance=0.85, startangle=90)
for text in texts:
    text.set_fontsize(9)
for autotext in autotexts:
    autotext.set_fontsize(8)
ax.set_title('Plasmid-Borne AMR Gene Classes\nin ST11 K. pneumoniae', fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/Figure5_AMR_Class_Distribution.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUTPUT_DIR}/Figure5_AMR_Class_Distribution.pdf", bbox_inches='tight')
plt.close()
print("  Saved Figure 5")

print(f"\n{'='*60}")
print(f"All 5 figures saved to: {OUTPUT_DIR}")
print(f"{'='*60}")
