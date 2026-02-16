#!/usr/bin/env python3
"""
Extract the Tn4401-like blaKPC-2 cassette sequence from all 8 ST11 K. pneumoniae genomes.
Fixes:
- JKPB284: file prefix is "2_JKPB284" not "3_JKPB284"  
- JKPB381: KPC is on a different contig (contig_3 = 105kb, not ~131kb). Need to check which contig.
- JKPR282: KPC is on contig_2, not contig_3
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

base_dir = "/home/ubuntu/plasmid_analysis"
lr_dir = os.path.join(base_dir, "LR_assemblies")
bakta_dir = os.path.join(base_dir, "04_AMR_Plasmid_Association", "Bakta_Annotations")
output_dir = os.path.join(base_dir, "09_Cassette_Sequences")
os.makedirs(output_dir, exist_ok=True)

# Correct file prefix mapping
isolates = {
    "JKPB244": "1_JKPB244",
    "JKPB284": "2_JKPB284",
    "JKPB381": "4_JKPB381",
    "JKPR282": "5_JKPR282",
    "KP21802": "6_KP21802",
    "KP21915": "7_KP21915",
    "MDRKP121": "8_MDRKP121",
    "MDRKP122": "3_MDRKP122",
}

def find_kpc_cassette(tsv_file):
    """Find the KPC-2 cassette boundaries from Bakta TSV - search ALL contigs."""
    all_genes = {}
    with open(tsv_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                contig = parts[0]
                if contig not in all_genes:
                    all_genes[contig] = []
                all_genes[contig].append({
                    'start': int(parts[2]),
                    'stop': int(parts[3]),
                    'strand': parts[4],
                    'gene': parts[6] if parts[6] else '',
                    'product': parts[7]
                })
    
    # Find which contig has KPC
    kpc_contig = None
    kpc_gene = None
    for contig, genes in all_genes.items():
        for g in genes:
            if 'KPC' in g['gene'].upper() or 'KPC' in g['product'].upper():
                kpc_contig = contig
                kpc_gene = g
                break
        if kpc_contig:
            break
    
    if kpc_contig is None:
        return None, None, None, []
    
    genes = sorted(all_genes[kpc_contig], key=lambda x: x['start'])
    
    # Find KPC index
    kpc_idx = None
    for i, g in enumerate(genes):
        if 'KPC' in g['gene'].upper() or 'KPC' in g['product'].upper():
            kpc_idx = i
            break
    
    # Find flanking IS15DIV
    cassette_start = genes[kpc_idx]['start']
    cassette_end = genes[kpc_idx]['stop']
    
    for i in range(kpc_idx, max(kpc_idx - 15, -1), -1):
        g = genes[i]
        if 'IS15DIV' in g['product'] or 'IS6 family' in g['product']:
            cassette_start = g['start']
            break
    
    for i in range(kpc_idx, min(kpc_idx + 15, len(genes))):
        g = genes[i]
        if i > kpc_idx and ('IS15DIV' in g['product'] or 'IS6 family' in g['product']):
            cassette_end = g['stop']
            break
    
    cassette_genes = [g for g in genes if g['start'] >= cassette_start and g['stop'] <= cassette_end]
    
    return kpc_contig, cassette_start, cassette_end, cassette_genes


def find_contig_in_fasta(fasta_file, contig_name, contig_number):
    """Find the matching contig in FASTA by name or number."""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Try exact match
    for r in records:
        if contig_name == r.id or contig_name in r.description:
            return r
    
    # Try by number (contig_3 = 3rd record, or record with id "3")
    num = contig_name.replace("contig_", "")
    for r in records:
        if r.id == num:
            return r
    
    # Try by index (contig_3 = index 2)
    try:
        idx = int(num) - 1
        if 0 <= idx < len(records):
            return records[idx]
    except:
        pass
    
    return None


all_cassettes = []
print("=" * 80)
print("EXTRACTING Tn4401-like blaKPC-2 CASSETTE FROM ALL 8 GENOMES")
print("=" * 80)

for isolate, prefix in sorted(isolates.items()):
    print(f"\n--- {isolate} ---")
    
    # Find TSV file
    tsv_file = os.path.join(bakta_dir, isolate, f"{isolate}.tsv")
    if not os.path.exists(tsv_file):
        alt_tsv = os.path.join(base_dir, "04_AMR_Plasmid_Association", f"{isolate}_bakta.tsv")
        if os.path.exists(alt_tsv):
            tsv_file = alt_tsv
        else:
            print(f"  WARNING: No TSV file found")
            continue
    
    # Find cassette
    kpc_contig, start, end, cassette_genes = find_kpc_cassette(tsv_file)
    if kpc_contig is None:
        print(f"  WARNING: blaKPC-2 not found")
        continue
    
    print(f"  KPC contig: {kpc_contig}")
    print(f"  Cassette: {start:,}-{end:,} bp ({end-start+1:,} bp)")
    print(f"  Genes: {len(cassette_genes)}")
    for g in cassette_genes:
        name = g['gene'] if g['gene'] else g['product'][:50]
        is_tag = " [IS]" if any(x in g['product'].lower() for x in ['transpos', 'is6', 'is15', 'iskpn', 'is481', 'is4', 'is110', 'is5075']) else ""
        amr_tag = " [AMR]" if any(x in (g['gene'] + g['product']).lower() for x in ['bla', 'kpc']) else ""
        print(f"    {g['start']:>8}-{g['stop']:>8} {g['strand']} {name}{is_tag}{amr_tag}")
    
    # Find FASTA
    fna_file = os.path.join(bakta_dir, isolate, f"{isolate}.fna")
    if not os.path.exists(fna_file):
        fna_file = os.path.join(lr_dir, f"{prefix}_LR.fasta")
    
    if not os.path.exists(fna_file):
        print(f"  WARNING: No FASTA file found")
        continue
    
    # Extract contig
    contig_num = kpc_contig.replace("contig_", "")
    contig = find_contig_in_fasta(fna_file, kpc_contig, contig_num)
    if contig is None:
        print(f"  WARNING: {kpc_contig} not found in FASTA")
        continue
    
    print(f"  Source: {contig.id} ({len(contig.seq):,} bp)")
    
    # Extract cassette sequence
    cassette_seq = contig.seq[start-1:end]
    
    gene_list = []
    for g in cassette_genes:
        name = g['gene'] if g['gene'] else g['product'].split()[0]
        gene_list.append(name)
    
    record = SeqRecord(
        cassette_seq,
        id=f"{isolate}_Tn4401_KPC2_cassette",
        name=f"{isolate}_KPC2",
        description=f"Tn4401-like blaKPC-2 cassette from {isolate} {kpc_contig}:{start}-{end} ({end-start+1}bp) genes:[{'|'.join(gene_list)}]"
    )
    
    all_cassettes.append(record)
    
    # Save individual FASTA
    out_file = os.path.join(output_dir, f"{isolate}_Tn4401_KPC2_cassette.fasta")
    SeqIO.write([record], out_file, "fasta")
    print(f"  Saved: {out_file}")

# Save combined
combined_file = os.path.join(output_dir, "All_Tn4401_KPC2_cassettes.fasta")
SeqIO.write(all_cassettes, combined_file, "fasta")

print(f"\n{'=' * 80}")
print(f"COMBINED FILE: {combined_file}")
print(f"Total cassettes extracted: {len(all_cassettes)}/8")
print(f"\nSUMMARY:")
for rec in all_cassettes:
    print(f"  {rec.id}: {len(rec.seq):,} bp")
print(f"{'=' * 80}")
