#!/usr/bin/env python3
"""
Sync all ST11 Plasmid Analysis results to Dropbox.
Uploads new/updated files from the local analysis directory.
"""
import os
import json
import requests

TOKEN = open('/home/ubuntu/.dropbox_token').read().strip()
LOCAL_BASE = "/home/ubuntu/plasmid_analysis"
DBX_BASE = "/ST11_Plasmid_Analysis"

# Map local dirs to Dropbox dirs
DIR_MAP = {
    '01_plasmid_enumeration': '01_Data_Assessment',
    '02_replicon_typing': '02_Replicon_Typing',
    '04_AMR_Plasmid_Association': '04_AMR_Plasmid_Association',
    '05_MGE_Analysis': '05_MGE_Analysis',
    '07_Visualizations': '07_Visualizations',
    '08_Reports': '08_Reports',
    '09_Scripts': '09_Scripts',
}

def upload_file(local_path, dbx_path):
    """Upload a file to Dropbox (overwrite if exists)."""
    url = "https://content.dropboxapi.com/2/files/upload"
    headers = {
        "Authorization": f"Bearer {TOKEN}",
        "Dropbox-API-Arg": json.dumps({
            "path": dbx_path,
            "mode": "overwrite",
            "autorename": False,
        }),
        "Content-Type": "application/octet-stream",
    }
    try:
        with open(local_path, 'rb') as f:
            response = requests.post(url, headers=headers, data=f, timeout=120)
        if response.status_code == 200:
            return True
        else:
            print(f"  FAIL: {dbx_path} - {response.status_code}: {response.text[:100]}")
            return False
    except Exception as e:
        print(f"  ERROR: {dbx_path} - {e}")
        return False

def main():
    print("=== Syncing all results to Dropbox ===\n")
    
    success = 0
    fail = 0
    
    # Upload files from mapped directories
    for local_dir, dbx_dir in DIR_MAP.items():
        local_path = os.path.join(LOCAL_BASE, local_dir)
        if not os.path.exists(local_path):
            continue
        
        for root, dirs, files in os.walk(local_path):
            for fname in files:
                fpath = os.path.join(root, fname)
                rel = os.path.relpath(fpath, local_path)
                dbx_path = f"{DBX_BASE}/{dbx_dir}/{rel}"
                
                if upload_file(fpath, dbx_path):
                    success += 1
                    print(f"  OK: {dbx_dir}/{rel}")
                else:
                    fail += 1
    
    # Upload key root-level files
    root_files = [
        'extract_amr_from_tsv.py',
        'create_visualizations.py',
        'recovery_status.md',
        'sync_all.py',
    ]
    for fname in root_files:
        fpath = os.path.join(LOCAL_BASE, fname)
        if os.path.exists(fpath):
            dbx_path = f"{DBX_BASE}/09_Scripts/{fname}"
            if upload_file(fpath, dbx_path):
                success += 1
                print(f"  OK: 09_Scripts/{fname}")
            else:
                fail += 1
    
    print(f"\n{'='*60}")
    print(f"Sync complete: {success} uploaded, {fail} failed")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
