#!/usr/bin/env python3
"""
Fast NCBI metadata enrichment using esummary + BioSample approach.
Extracts country from BLAST hit descriptions and NCBI esummary.
"""

import pandas as pd
import re
import json
import time
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

NCBI_EMAIL = "analysis@example.com"

def fetch_url(url, max_retries=3):
    for attempt in range(max_retries):
        try:
            req = Request(url, headers={'User-Agent': 'Python/3.11'})
            with urlopen(req, timeout=120) as response:
                return response.read().decode('utf-8')
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(2 * (attempt + 1))
    return None

def batch_esummary_nucleotide(accessions, batch_size=100):
    """Use esummary which is much faster than efetch for basic metadata."""
    results = {}
    total = len(accessions)
    
    for i in range(0, total, batch_size):
        batch = accessions[i:i+batch_size]
        ids_str = ','.join(batch)
        
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={ids_str}&retmode=json&email={NCBI_EMAIL}"
        
        batch_num = i // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size
        print(f"  esummary batch {batch_num}/{total_batches}...")
        
        text = fetch_url(url)
        if text:
            try:
                data = json.loads(text)
                if 'result' in data:
                    for uid in data['result'].get('uids', []):
                        rec = data['result'].get(uid, {})
                        acc = rec.get('caption', uid)
                        results[acc] = {
                            'title': rec.get('title', ''),
                            'organism': rec.get('organism', ''),
                            'strain': rec.get('strain', ''),
                            'biosample': rec.get('biosample', ''),
                            'subtype': rec.get('subtype', ''),
                            'subname': rec.get('subname', ''),
                            'extra': rec.get('extra', ''),
                        }
            except json.JSONDecodeError:
                print(f"    JSON parse error for batch {batch_num}")
        
        time.sleep(0.4)
    
    return results

def extract_country_from_description(desc):
    """Try to extract country from NCBI description/title."""
    # Common patterns in plasmid descriptions
    country_patterns = {
        'China': r'(?:China|Chinese|Beijing|Shanghai|Hangzhou|Guangzhou|Shenzhen|Wuhan|Chongqing|Nanjing|Sichuan|Zhejiang|Jiangsu|Guangdong|Henan|Hubei|Hunan|Fujian|Shandong|Anhui|Hebei|Liaoning|Jilin|Heilongjiang|Yunnan|Guizhou|Shanxi|Shaanxi|Gansu|Qinghai|Hainan|Taiwan|Hong Kong|Tianjin|Chengdu|Changsha|Zhengzhou|Jinan|Kunming|Nanning|Harbin|Dalian|Xiamen|Qingdao|Suzhou|Wenzhou|Ningbo|Fuzhou|Hefei|Wuxi|Dongguan|Foshan|Changchun|Urumqi|Guiyang|Nanchang|Lanzhou|Hohhot|Yinchuan|Xining|Lhasa|HN|GD|ZJ|JS|SD|AH|HB|FJ|SC|CQ|BJ|SH|TJ|GZ|SZ|NJ|WH|CD|CS|ZZ|JN|KM|NN|HRB|DL|XM|QD|SZ|WZ|NB|FZ|HF|WX)',
        'USA': r'(?:United States|USA|U\.S\.|American|New York|California|Texas|Pennsylvania|Illinois|Ohio|Georgia|North Carolina|Michigan|New Jersey|Virginia|Washington|Arizona|Massachusetts|Tennessee|Indiana|Maryland|Wisconsin|Colorado|Minnesota|South Carolina|Alabama|Louisiana|Kentucky|Oregon|Oklahoma|Connecticut|Utah|Iowa|Nevada|Arkansas|Mississippi|Kansas|New Mexico|Nebraska|Idaho|West Virginia|Hawaii|New Hampshire|Maine|Montana|Rhode Island|Delaware|South Dakota|North Dakota|Alaska|Vermont|Wyoming)',
        'India': r'(?:India|Indian|Mumbai|Delhi|Bangalore|Hyderabad|Chennai|Kolkata|Pune|Ahmedabad|Jaipur)',
        'Brazil': r'(?:Brazil|Brazilian|São Paulo|Rio de Janeiro|Brasilia)',
        'South Korea': r'(?:South Korea|Korean|Seoul|Busan|Korea)',
        'Japan': r'(?:Japan|Japanese|Tokyo|Osaka|Kyoto)',
        'Italy': r'(?:Italy|Italian|Rome|Milan)',
        'Germany': r'(?:Germany|German|Berlin|Munich)',
        'France': r'(?:France|French|Paris)',
        'UK': r'(?:United Kingdom|UK|England|London|British)',
        'Spain': r'(?:Spain|Spanish|Madrid|Barcelona)',
        'Australia': r'(?:Australia|Australian|Sydney|Melbourne)',
        'Canada': r'(?:Canada|Canadian|Toronto|Vancouver)',
        'Thailand': r'(?:Thailand|Thai|Bangkok)',
        'Vietnam': r'(?:Vietnam|Vietnamese|Hanoi|Ho Chi Minh)',
        'Pakistan': r'(?:Pakistan|Pakistani|Karachi|Lahore)',
        'Egypt': r'(?:Egypt|Egyptian|Cairo)',
        'Turkey': r'(?:Turkey|Turkish|Istanbul|Ankara)',
        'Iran': r'(?:Iran|Iranian|Tehran)',
        'Colombia': r'(?:Colombia|Colombian|Bogota)',
        'Argentina': r'(?:Argentina|Argentine|Buenos Aires)',
        'Mexico': r'(?:Mexico|Mexican|Mexico City)',
        'Russia': r'(?:Russia|Russian|Moscow)',
        'Poland': r'(?:Poland|Polish|Warsaw)',
        'Greece': r'(?:Greece|Greek|Athens)',
        'Portugal': r'(?:Portugal|Portuguese|Lisbon)',
        'Netherlands': r'(?:Netherlands|Dutch|Amsterdam)',
        'Belgium': r'(?:Belgium|Belgian|Brussels)',
        'Switzerland': r'(?:Switzerland|Swiss|Zurich)',
        'Sweden': r'(?:Sweden|Swedish|Stockholm)',
        'Norway': r'(?:Norway|Norwegian|Oslo)',
        'Denmark': r'(?:Denmark|Danish|Copenhagen)',
        'Finland': r'(?:Finland|Finnish|Helsinki)',
        'Austria': r'(?:Austria|Austrian|Vienna)',
        'Czech Republic': r'(?:Czech|Prague)',
        'Romania': r'(?:Romania|Romanian|Bucharest)',
        'Hungary': r'(?:Hungary|Hungarian|Budapest)',
        'Israel': r'(?:Israel|Israeli|Tel Aviv|Jerusalem)',
        'Saudi Arabia': r'(?:Saudi Arabia|Saudi|Riyadh|Jeddah)',
        'UAE': r'(?:UAE|United Arab Emirates|Dubai|Abu Dhabi)',
        'Singapore': r'(?:Singapore|Singaporean)',
        'Malaysia': r'(?:Malaysia|Malaysian|Kuala Lumpur)',
        'Philippines': r'(?:Philippines|Filipino|Manila)',
        'Indonesia': r'(?:Indonesia|Indonesian|Jakarta)',
        'Nigeria': r'(?:Nigeria|Nigerian|Lagos)',
        'South Africa': r'(?:South Africa|Johannesburg|Cape Town)',
        'Kenya': r'(?:Kenya|Kenyan|Nairobi)',
        'Chile': r'(?:Chile|Chilean|Santiago)',
        'Peru': r'(?:Peru|Peruvian|Lima)',
        'Ecuador': r'(?:Ecuador|Ecuadorian|Quito)',
        'Venezuela': r'(?:Venezuela|Venezuelan|Caracas)',
        'Cuba': r'(?:Cuba|Cuban|Havana)',
        'Dominican Republic': r'(?:Dominican Republic|Santo Domingo)',
        'Guatemala': r'(?:Guatemala|Guatemalan)',
        'Costa Rica': r'(?:Costa Rica|Costa Rican)',
        'Panama': r'(?:Panama|Panamanian)',
        'New Zealand': r'(?:New Zealand|Auckland|Wellington)',
        'Bangladesh': r'(?:Bangladesh|Bangladeshi|Dhaka)',
        'Sri Lanka': r'(?:Sri Lanka|Lankan|Colombo)',
        'Nepal': r'(?:Nepal|Nepalese|Kathmandu)',
        'Myanmar': r'(?:Myanmar|Burmese|Yangon)',
        'Cambodia': r'(?:Cambodia|Cambodian|Phnom Penh)',
        'Laos': r'(?:Laos|Laotian|Vientiane)',
        'Lebanon': r'(?:Lebanon|Lebanese|Beirut)',
        'Jordan': r'(?:Jordan|Jordanian|Amman)',
        'Iraq': r'(?:Iraq|Iraqi|Baghdad)',
        'Algeria': r'(?:Algeria|Algerian|Algiers)',
        'Morocco': r'(?:Morocco|Moroccan|Rabat)',
        'Tunisia': r'(?:Tunisia|Tunisian|Tunis)',
        'Libya': r'(?:Libya|Libyan|Tripoli)',
        'Ethiopia': r'(?:Ethiopia|Ethiopian|Addis Ababa)',
        'Tanzania': r'(?:Tanzania|Tanzanian|Dar es Salaam)',
        'Ghana': r'(?:Ghana|Ghanaian|Accra)',
        'Cameroon': r'(?:Cameroon|Cameroonian|Yaoundé)',
    }
    
    for country, pattern in country_patterns.items():
        if re.search(pattern, desc, re.IGNORECASE):
            return country
    return 'Not specified'

def extract_country_from_subtype(subtype, subname):
    """Extract country from esummary subtype/subname fields."""
    if not subtype or not subname:
        return None
    
    # subtype and subname are pipe-delimited parallel fields
    types = str(subtype).split('|')
    names = str(subname).split('|')
    
    for t, n in zip(types, names):
        t = t.strip()
        n = n.strip()
        if t.lower() in ['country', 'geo_loc_name']:
            # Extract just the country part (before colon)
            country = n.split(':')[0].strip()
            return country
        if t.lower() == 'isolation_source':
            pass  # Could use this too
    
    return None

def main():
    print("=" * 70)
    print("Fast NCBI Metadata Enrichment for KPC-2 Cassette BLAST Hits")
    print("=" * 70)
    
    # Load BLAST results
    df = pd.read_csv('/home/ubuntu/plasmid_analysis/10_BLAST_Results/blast_kpc_cassette_results.csv')
    hq = df[(df['percent_identity'] >= 90) & (df['query_coverage'] >= 80)].copy()
    
    print(f"High-quality hits: {len(hq)}")
    accessions = hq['accession'].unique().tolist()
    print(f"Unique accessions: {len(accessions)}")
    
    # Step 1: Try esummary for metadata
    print("\nStep 1: Fetching esummary metadata...")
    esummary_data = batch_esummary_nucleotide(accessions, batch_size=100)
    print(f"  Got esummary for {len(esummary_data)} accessions")
    
    # Step 2: Extract countries from multiple sources
    print("\nStep 2: Extracting geographic information...")
    
    country_map = {}
    host_map = {}
    collection_date_map = {}
    isolation_source_map = {}
    
    for acc in accessions:
        country = None
        
        # Try esummary subtype/subname first (most reliable)
        if acc in esummary_data:
            es = esummary_data[acc]
            country = extract_country_from_subtype(es.get('subtype', ''), es.get('subname', ''))
            
            # Also extract host and date from subtype/subname
            types = str(es.get('subtype', '')).split('|')
            names = str(es.get('subname', '')).split('|')
            for t, n in zip(types, names):
                t = t.strip().lower()
                n = n.strip()
                if t == 'host':
                    host_map[acc] = n
                if t == 'collection_date':
                    collection_date_map[acc] = n
                if t == 'isolation_source':
                    isolation_source_map[acc] = n
        
        # If no country from esummary, try description
        if not country:
            desc_rows = hq[hq['accession'] == acc]
            if len(desc_rows) > 0:
                desc = str(desc_rows.iloc[0]['description'])
                country = extract_country_from_description(desc)
                if country == 'Not specified':
                    # Try from esummary title
                    if acc in esummary_data:
                        title = esummary_data[acc].get('title', '')
                        country = extract_country_from_description(title)
        
        country_map[acc] = country if country and country != 'Not specified' else 'Not specified'
    
    # Step 3: Add metadata to dataframe
    print("\nStep 3: Merging metadata...")
    hq['geographic_origin'] = hq['accession'].map(country_map).fillna('Not specified')
    hq['host_organism'] = hq['accession'].map(host_map).fillna('Not specified')
    hq['collection_date'] = hq['accession'].map(collection_date_map).fillna('Not specified')
    hq['isolation_source_detail'] = hq['accession'].map(isolation_source_map).fillna('Not specified')
    
    # Save enriched results
    outdir = '/home/ubuntu/plasmid_analysis/10_BLAST_Results'
    hq.to_csv(f'{outdir}/blast_kpc_enriched_results.csv', index=False)
    
    # Save metadata
    metadata = {
        'countries': country_map,
        'hosts': host_map,
        'collection_dates': collection_date_map,
        'isolation_sources': isolation_source_map,
        'esummary': esummary_data
    }
    with open(f'{outdir}/ncbi_metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    
    # Print summary
    print("\n" + "=" * 70)
    print("GEOGRAPHIC DISTRIBUTION")
    print("=" * 70)
    geo = hq.drop_duplicates('accession')['geographic_origin'].value_counts()
    specified = geo[geo.index != 'Not specified']
    print(f"Countries identified: {len(specified)}")
    print(f"Accessions with country: {specified.sum()}/{len(accessions)}")
    for c, n in specified.head(25).items():
        print(f"  {c}: {n}")
    
    print("\n" + "=" * 70)
    print("SPECIES DISTRIBUTION")
    print("=" * 70)
    sp = hq.drop_duplicates('accession')['species'].value_counts()
    for s, n in sp.head(20).items():
        print(f"  {s}: {n}")
    
    print("\n" + "=" * 70)
    print("SOURCE TYPE DISTRIBUTION")
    print("=" * 70)
    st = hq.drop_duplicates('accession')['source_type'].value_counts()
    for s, n in st.items():
        print(f"  {s}: {n}")
    
    print("\n" + "=" * 70)
    print("HOST ORGANISM DISTRIBUTION")
    print("=" * 70)
    ho = pd.Series(host_map).value_counts()
    for h, n in ho.head(10).items():
        print(f"  {h}: {n}")
    
    print("\n" + "=" * 70)
    print("PLASMID NAME DISTRIBUTION (top 20)")
    print("=" * 70)
    pn = hq.drop_duplicates('accession')[hq.drop_duplicates('accession')['source_type'] == 'Plasmid']['plasmid_name'].value_counts()
    for p, n in pn.head(20).items():
        print(f"  {p}: {n}")
    
    print(f"\nTotal unique plasmid names: {len(pn)}")
    print(f"\nDone! Results saved to {outdir}/")

if __name__ == '__main__':
    main()
