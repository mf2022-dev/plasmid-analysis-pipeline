#!/usr/bin/env python3
"""
Download reference KPC-carrying plasmid sequences from NCBI for comparison.
Select diverse geographic representatives from the 100% cassette identity matches.
"""

import os
import time
import json
from urllib.request import urlopen, Request

NCBI_EMAIL = "analysis@example.com"
OUTDIR = "/home/ubuntu/plasmid_analysis/12_Plasmid_Comparison/reference_plasmids"
os.makedirs(OUTDIR, exist_ok=True)

# Select diverse references:
# - CP097693: 130.7 kb, China (closest size to our ~130.7 kb)
# - CP082020: 135.8 kb, China (pKP18-2050-KPC2, named)
# - MT648513: 129.9 kb, Egypt (geographic diversity)
# - CP192284: 129.5 kb, United Kingdom (geographic diversity)
# - CP025458: 128.6 kb, China (p69-2, well-characterized)
# - CP045264: 129.3 kb, China (p16HN-263_KPC, named)
# - MZ709016: 127.7 kb, China (pJNKPN52_kpc_fosA)
# - CP047161: 127.2 kb, China (pKP19-2029-KPC2)
# - KX236178: 121.3 kb, China (pHS091147)
# - CP073303: 125.6 kb, China (pSW1780-KPC)

references = {
    "CP097693": {"name": "unnamed1_China_130kb", "size": 130704, "origin": "China"},
    "CP082020": {"name": "pKP18-2050-KPC2_China_136kb", "size": 135835, "origin": "China"},
    "MT648513": {"name": "pEBSI036-2-KPC_Egypt_130kb", "size": 129877, "origin": "Egypt"},
    "CP192284": {"name": "p-unmaned1_UK_130kb", "size": 129488, "origin": "UK"},
    "CP025458": {"name": "p69-2_China_129kb", "size": 128619, "origin": "China"},
    "CP045264": {"name": "p16HN-263_KPC_China_129kb", "size": 129314, "origin": "China"},
    "MZ709016": {"name": "pJNKPN52_kpc_fosA_China_128kb", "size": 127743, "origin": "China"},
    "CP047161": {"name": "pKP19-2029-KPC2_China_127kb", "size": 127176, "origin": "China"},
    "KX236178": {"name": "pHS091147_China_121kb", "size": 121277, "origin": "China"},
    "CP073303": {"name": "pSW1780-KPC_China_126kb", "size": 125575, "origin": "China"},
}

def download_fasta(accession, outfile):
    """Download FASTA sequence from NCBI."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=fasta&retmode=text&email={NCBI_EMAIL}"
    try:
        req = Request(url, headers={'User-Agent': 'Python/3.11'})
        with urlopen(req, timeout=120) as response:
            data = response.read().decode('utf-8')
        
        if data.startswith('>'):
            with open(outfile, 'w') as f:
                f.write(data)
            # Count sequence length
            seq = ''.join(line.strip() for line in data.split('\n') if not line.startswith('>'))
            return len(seq)
        else:
            print(f"  ERROR: Invalid FASTA for {accession}")
            return 0
    except Exception as e:
        print(f"  ERROR downloading {accession}: {e}")
        return 0

print("Downloading reference plasmid sequences from NCBI...")
print("=" * 60)

downloaded = {}
for acc, info in references.items():
    outfile = os.path.join(OUTDIR, f"{acc}_{info['name']}.fasta")
    print(f"  {acc} ({info['name']})...", end=" ", flush=True)
    
    size = download_fasta(acc, outfile)
    if size > 0:
        print(f"OK ({size:,} bp)")
        downloaded[acc] = {
            'file': outfile,
            'size': size,
            'name': info['name'],
            'origin': info['origin']
        }
    else:
        print("FAILED")
    
    time.sleep(0.5)

print(f"\nDownloaded {len(downloaded)}/{len(references)} reference plasmids")

# Save metadata
with open(os.path.join(OUTDIR, 'reference_metadata.json'), 'w') as f:
    json.dump(downloaded, f, indent=2)

print("Done!")
