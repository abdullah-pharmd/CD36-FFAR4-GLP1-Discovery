import os
import urllib.request
import re

STRUCTURES_DIR = r"C:\Users\Abdullah\.gemini\antigravity\scratch\cd36_ffar4_glp1_study\structures"
os.makedirs(STRUCTURES_DIR, exist_ok=True)

PDB_IDS = ["5LGD", "8G59", "8ID6"]

def download_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    raw_path = os.path.join(STRUCTURES_DIR, f"{pdb_id}.pdb")
    print(f"Downloading {pdb_id} from {url}...")
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req) as response, open(raw_path, 'wb') as out_file:
            data = response.read()
            out_file.write(data)
        size = os.path.getsize(raw_path)
        print(f"Successfully downloaded {pdb_id}.pdb ({size:,} bytes)")
        return raw_path
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        return None

def inspect_and_clean_pdb(pdb_id, raw_path):
    print(f"\n--- Analyzing {pdb_id} ---")
    with open(raw_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    # Identify Title, Molecules, Chains, Heteroatoms
    title = ""
    compounds = []
    hetatms = set()
    chains = set()
    res_per_chain = {}
    
    for line in lines:
        if line.startswith("TITLE"):
            title += line[10:].strip() + " "
        elif line.startswith("COMPND"):
            compounds.append(line.strip())
        elif line.startswith("ATOM  "):
            ch = line[21:22].strip()
            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            chains.add(ch)
            if ch not in res_per_chain:
                res_per_chain[ch] = set()
            res_per_chain[ch].add((res_name, res_num))
        elif line.startswith("HETATM"):
            res_name = line[17:20].strip()
            ch = line[21:22].strip()
            hetatms.add((res_name, ch))

    print(f"Title: {title.strip()}")
    print(f"Chains present: {sorted(list(chains))}")
    for ch in sorted(list(chains)):
        print(f"  Chain {ch}: {len(res_per_chain[ch])} residues")
    print(f"Heteroatoms (res, chain): {sorted(list(hetatms))}")

if __name__ == "__main__":
    for pid in PDB_IDS:
        p = download_pdb(pid)
        if p:
            inspect_and_clean_pdb(pid, p)
