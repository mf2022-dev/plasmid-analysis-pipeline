#!/usr/bin/env python3
"""
In silico replicon typing of BLAST hit plasmids.

Strategy:
1. Download PlasmidFinder Enterobacteriaceae replicon sequences
2. For each BLAST hit plasmid accession, search for replicon sequences using NCBI BLAST
3. Alternatively: use plasmid size + name clustering as a proxy

Since running 460 individual BLAST searches is impractical, we use a hybrid approach:
- Parse plasmid names for known naming conventions that indicate replicon type
- Use plasmid size clustering to group similar plasmids
- Cross-reference with published literature on KPC-carrying plasmid types
- Use NCBI Assembly/BioSample links to find associated PlasmidFinder results
"""

import pandas as pd
import numpy as np
import re
import json
import os
import time
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from collections import Counter

def fetch_url(url, max_retries=3, timeout=120):
    for attempt in range(max_retries):
        try:
            req = Request(url, headers={'User-Agent': 'Python/3.11'})
            with urlopen(req, timeout=timeout) as response:
                return response.read().decode('utf-8')
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(2 * (attempt + 1))
    return None

def classify_plasmid_by_name(name):
    """Classify plasmid replicon type based on naming conventions in literature."""
    name_lower = str(name).lower()
    
    # Direct Inc mentions
    if re.search(r'incfii.*incr|incr.*incfii', name_lower):
        return 'IncFII+IncR'
    if re.search(r'incfii', name_lower):
        return 'IncFII'
    if re.search(r'incfib.*inchi1b|inchi1b.*incfib', name_lower):
        return 'IncFIB+IncHI1B'
    if re.search(r'incn', name_lower):
        return 'IncN'
    if re.search(r'incl|incm', name_lower):
        return 'IncL/M'
    if re.search(r'incx', name_lower):
        return 'IncX'
    if re.search(r'inchi', name_lower):
        return 'IncHI'
    if re.search(r'incr\b', name_lower):
        return 'IncR'
    if re.search(r'incfib', name_lower):
        return 'IncFIB'
    
    return None

def classify_by_size_and_context(size, species, name):
    """
    Classify plasmid replicon type based on size and known KPC plasmid epidemiology.
    
    Key literature-based size ranges for KPC-carrying plasmids:
    - IncFII(pHN7A8)+IncR: ~80-140 kb (most common in ST11 K. pneumoniae in China)
    - IncFII: ~80-140 kb
    - IncN: ~40-70 kb  
    - IncL/M: ~30-60 kb
    - IncX3: ~45-55 kb
    - IncP: ~50-70 kb
    - IncFIB+IncHI1B: ~200-400 kb (large conjugative, often NDM)
    - ColKP3: ~10-20 kb (small, non-conjugative)
    """
    size_kb = size / 1000
    
    # Very small plasmids
    if size_kb < 30:
        return 'Small mobilizable (<30 kb)'
    
    # Small-medium
    if 30 <= size_kb < 70:
        return 'Medium (30-70 kb; likely IncN/IncL/IncX)'
    
    # The dominant KPC plasmid size range
    if 70 <= size_kb < 160:
        return 'Large (70-160 kb; likely IncFII±IncR)'
    
    # Large plasmids
    if 160 <= size_kb < 250:
        return 'Very large (160-250 kb)'
    
    # Mega plasmids
    if size_kb >= 250:
        return 'Mega (>250 kb; likely IncHI/IncFIB)'
    
    return 'Unclassified'

def fetch_biosample_replicon_data(accessions, batch_size=100):
    """
    Try to get replicon typing from linked BioSample/Assembly records.
    Many genome projects include PlasmidFinder results in their metadata.
    """
    results = {}
    total = len(accessions)
    
    for i in range(0, total, batch_size):
        batch = accessions[i:i+batch_size]
        ids_str = ','.join(batch)
        batch_num = i // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size
        
        print(f"  Batch {batch_num}/{total_batches}...")
        
        # Use elink to find linked BioSample records
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={ids_str}&retmode=json&email=analysis@example.com"
        
        text = fetch_url(url)
        if text:
            try:
                data = json.loads(text)
                if 'result' in data:
                    for uid in data['result'].get('uids', []):
                        rec = data['result'].get(uid, {})
                        acc = rec.get('caption', uid)
                        
                        # Check subtype for plasmid info
                        subtype = str(rec.get('subtype', ''))
                        subname = str(rec.get('subname', ''))
                        title = str(rec.get('title', ''))
                        
                        # Extract replicon info from title
                        inc_matches = re.findall(r'(Inc[A-Z][A-Za-z0-9_\(\)\-]+)', title)
                        if inc_matches:
                            results[acc] = inc_matches
                        
                        # Also check subname (plasmid name field)
                        types = subtype.split('|')
                        names = subname.split('|')
                        for t, n in zip(types, names):
                            if t.strip() == 'plasmid':
                                name_inc = re.findall(r'(Inc[A-Z][A-Za-z0-9_\(\)\-]+)', n)
                                if name_inc and acc not in results:
                                    results[acc] = name_inc
            except json.JSONDecodeError:
                pass
        
        time.sleep(0.4)
    
    return results

def main():
    print("=" * 70)
    print("Comprehensive Replicon Typing of KPC-2 Cassette Carrier Plasmids")
    print("=" * 70)
    
    # Load data
    df = pd.read_csv('/home/ubuntu/plasmid_analysis/10_BLAST_Results/blast_kpc_enriched_results.csv')
    plasmid_df = df[df['source_type'] == 'Plasmid'].drop_duplicates('accession').copy()
    print(f"\nPlasmid-borne hits: {len(plasmid_df)} unique accessions")
    
    accessions = plasmid_df['accession'].tolist()
    
    # Step 1: Try to get replicon types from NCBI metadata
    print("\nStep 1: Querying NCBI for replicon metadata...")
    ncbi_replicons = fetch_biosample_replicon_data(accessions, batch_size=100)
    print(f"  Found replicon info for {len(ncbi_replicons)} accessions from NCBI")
    
    # Step 2: Classify by plasmid name
    print("\nStep 2: Classifying by plasmid name patterns...")
    name_classified = 0
    
    # Step 3: Classify by size
    print("Step 3: Classifying by plasmid size...")
    
    # Build comprehensive classification
    results = []
    for idx, row in plasmid_df.iterrows():
        acc = row['accession']
        name = row['plasmid_name']
        size = row['hit_length']
        species = row['species']
        
        # Priority 1: NCBI metadata
        replicon = None
        source = None
        
        if acc in ncbi_replicons:
            replicon = '+'.join(sorted(set(ncbi_replicons[acc])))
            source = 'NCBI metadata'
        
        # Priority 2: Plasmid name
        if not replicon:
            name_rep = classify_plasmid_by_name(name)
            if name_rep:
                replicon = name_rep
                source = 'Plasmid name'
                name_classified += 1
        
        # Priority 3: Size-based classification
        size_class = classify_by_size_and_context(size, species, name)
        
        if not replicon:
            replicon = size_class
            source = 'Size-based inference'
        
        results.append({
            'accession': acc,
            'species': species,
            'plasmid_name': name,
            'plasmid_size_bp': size,
            'plasmid_size_kb': round(size / 1000, 1),
            'replicon_type': replicon,
            'classification_source': source,
            'size_class': size_class,
            'geographic_origin': row['geographic_origin'],
            'percent_identity': row['percent_identity'],
            'host_organism': row.get('host_organism', 'Not specified'),
            'collection_date': row.get('collection_date', 'Not specified'),
        })
    
    rep_df = pd.DataFrame(results)
    
    print(f"\n  Name-classified: {name_classified}")
    print(f"  NCBI-classified: {len(ncbi_replicons)}")
    print(f"  Size-classified: {(rep_df['classification_source'] == 'Size-based inference').sum()}")
    
    # ============================================================
    # ANALYSIS
    # ============================================================
    
    print("\n" + "=" * 70)
    print("REPLICON TYPE DISTRIBUTION")
    print("=" * 70)
    
    rep_counts = rep_df['replicon_type'].value_counts()
    for t, n in rep_counts.items():
        pct = n / len(rep_df) * 100
        print(f"  {t}: {n} ({pct:.1f}%)")
    
    print("\n" + "=" * 70)
    print("SIZE CLASS DISTRIBUTION")
    print("=" * 70)
    
    size_counts = rep_df['size_class'].value_counts()
    for t, n in size_counts.items():
        pct = n / len(rep_df) * 100
        print(f"  {t}: {n} ({pct:.1f}%)")
    
    print("\n" + "=" * 70)
    print("CLASSIFICATION SOURCE")
    print("=" * 70)
    
    src_counts = rep_df['classification_source'].value_counts()
    for t, n in src_counts.items():
        pct = n / len(rep_df) * 100
        print(f"  {t}: {n} ({pct:.1f}%)")
    
    print("\n" + "=" * 70)
    print("SIZE CLASS × SPECIES")
    print("=" * 70)
    
    top_sp = rep_df['species'].value_counts().head(8).index
    cross = pd.crosstab(rep_df[rep_df['species'].isin(top_sp)]['size_class'], 
                        rep_df[rep_df['species'].isin(top_sp)]['species'])
    print(cross.to_string())
    
    print("\n" + "=" * 70)
    print("SIZE CLASS × GEOGRAPHY")
    print("=" * 70)
    
    top_geo = rep_df[rep_df['geographic_origin'] != 'Not specified']['geographic_origin'].value_counts().head(6).index
    geo_cross = pd.crosstab(rep_df[rep_df['geographic_origin'].isin(top_geo)]['size_class'],
                            rep_df[rep_df['geographic_origin'].isin(top_geo)]['geographic_origin'])
    print(geo_cross.to_string())
    
    print("\n" + "=" * 70)
    print("PLASMID SIZE STATISTICS")
    print("=" * 70)
    print(f"  Overall median: {rep_df['plasmid_size_kb'].median():.1f} kb")
    print(f"  Overall mean: {rep_df['plasmid_size_kb'].mean():.1f} kb")
    print(f"  Range: {rep_df['plasmid_size_kb'].min():.1f} - {rep_df['plasmid_size_kb'].max():.1f} kb")
    print(f"  IQR: {rep_df['plasmid_size_kb'].quantile(0.25):.1f} - {rep_df['plasmid_size_kb'].quantile(0.75):.1f} kb")
    
    # Key finding: what % are in the ~100-140 kb range (typical IncFII+IncR)?
    incfii_range = rep_df[(rep_df['plasmid_size_kb'] >= 90) & (rep_df['plasmid_size_kb'] <= 145)]
    print(f"\n  Plasmids in 90-145 kb range (typical IncFII+IncR): {len(incfii_range)} ({len(incfii_range)/len(rep_df)*100:.1f}%)")
    
    # Our study plasmid is ~131 kb
    similar_to_ours = rep_df[(rep_df['plasmid_size_kb'] >= 120) & (rep_df['plasmid_size_kb'] <= 140)]
    print(f"  Plasmids in 120-140 kb range (similar to our ~131 kb): {len(similar_to_ours)} ({len(similar_to_ours)/len(rep_df)*100:.1f}%)")
    
    # Save results
    outdir = '/home/ubuntu/plasmid_analysis/11_Replicon_Association'
    os.makedirs(outdir, exist_ok=True)
    
    rep_df.to_csv(f'{outdir}/replicon_typing_results.csv', index=False)
    
    summary = {
        'total_plasmid_accessions': len(rep_df),
        'ncbi_classified': len(ncbi_replicons),
        'name_classified': name_classified,
        'size_classified': int((rep_df['classification_source'] == 'Size-based inference').sum()),
        'replicon_distribution': {k: int(v) for k, v in rep_counts.items()},
        'size_distribution': {k: int(v) for k, v in size_counts.items()},
        'median_size_kb': float(rep_df['plasmid_size_kb'].median()),
        'mean_size_kb': float(rep_df['plasmid_size_kb'].mean()),
        'pct_90_145kb': float(len(incfii_range) / len(rep_df) * 100),
        'pct_120_140kb': float(len(similar_to_ours) / len(rep_df) * 100),
    }
    
    with open(f'{outdir}/replicon_typing_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nResults saved to {outdir}/")

if __name__ == '__main__':
    main()
