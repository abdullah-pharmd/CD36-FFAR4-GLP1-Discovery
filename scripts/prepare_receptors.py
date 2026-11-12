import os
import numpy as np

STRUCTURES_DIR = r"C:\Users\Abdullah\.gemini\antigravity\scratch\cd36_ffar4_glp1_study\structures"

STANDARD_AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE'  # Selenomethionine
}

def clean_pdb(raw_pdb, clean_pdb, chain_id, target_name):
    """Extract standard amino acid residues for a given chain, excluding waters and heteroatoms."""
    raw_path = os.path.join(STRUCTURES_DIR, raw_pdb)
    clean_path = os.path.join(STRUCTURES_DIR, clean_pdb)
    
    atoms = []
    residues = set()
    coords = []
    
    with open(raw_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith("ATOM  "):
                ch = line[21:22].strip()
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                if ch == chain_id and res_name in STANDARD_AMINO_ACIDS:
                    atoms.append(line)
                    residues.add((res_name, res_seq))
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
    
    with open(clean_path, 'w', encoding='utf-8') as f:
        f.write(f"REMARK   Cleaned structure: {target_name}\n")
        f.write(f"REMARK   Source: {raw_pdb} Chain {chain_id}\n")
        f.write(f"REMARK   Total residues: {len(residues)}, Total atoms: {len(atoms)}\n")
        for at in atoms:
            f.write(at)
        f.write("END\n")
    
    coords = np.array(coords)
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)
    
    print(f"[RECEPTOR] {target_name} ({clean_pdb}):")
    print(f"  Chain: {chain_id} | Residues: {len(residues)} | Atoms: {len(atoms)}")
    print(f"  Receptor geometric center: x={center[0]:.3f}, y={center[1]:.3f}, z={center[2]:.3f}")
    print(f"  Receptor dimension spans: dx={span[0]:.2f}, dy={span[1]:.2f}, dz={span[2]:.2f}")
    return clean_path, center, span

def extract_ligand(raw_pdb, out_pdb, chain_id, res_name, res_seq=None, ligand_name=""):
    """Extract specific HETATM ligand coordinates for redocking validation and pocket grid center."""
    raw_path = os.path.join(STRUCTURES_DIR, raw_pdb)
    out_path = os.path.join(STRUCTURES_DIR, out_pdb)
    
    lines = []
    coords = []
    with open(raw_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith("HETATM"):
                r_name = line[17:20].strip()
                ch = line[21:22].strip()
                r_seq = line[22:26].strip()
                if ch == chain_id and r_name == res_name:
                    if res_seq is None or r_seq == str(res_seq):
                        lines.append(line)
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
    
    if lines:
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write(f"REMARK   Extracted Ligand: {ligand_name} ({res_name})\n")
            f.write(f"REMARK   Source: {raw_pdb} Chain {chain_id} Res {res_name} {res_seq if res_seq else ''}\n")
            for l in lines:
                f.write(l)
            f.write("END\n")
        coords = np.array(coords)
        center = coords.mean(axis=0)
        span = coords.max(axis=0) - coords.min(axis=0)
        print(f"[LIGAND] {ligand_name} ({res_name}) -> {out_pdb}:")
        print(f"  Atoms: {len(lines)} | Pocket Center: x={center[0]:.3f}, y={center[1]:.3f}, z={center[2]:.3f}")
        print(f"  Recommended Grid Box: center_x={center[0]:.2f}, center_y={center[1]:.2f}, center_z={center[2]:.2f}, size ~22x22x22 A")
        return center
    else:
        print(f"[WARN] No ligand found for {res_name} in {raw_pdb}")
        return None

def find_cd36_tunnel_entrance(raw_pdb, chain_id='A'):
    """Locate CD36 Lys164 and Lys166 to define the extracellular tunnel entrance grid box."""
    raw_path = os.path.join(STRUCTURES_DIR, raw_pdb)
    lys_coords = []
    with open(raw_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith("ATOM  "):
                ch = line[21:22].strip()
                res_seq = line[22:26].strip()
                if ch == chain_id and res_seq in ['164', '166', '335']:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    lys_coords.append((x, y, z))
    if lys_coords:
        lys_coords = np.array(lys_coords)
        center = lys_coords.mean(axis=0)
        print(f"[CD36 TUNNEL ENTRANCE] Lys164/Lys166/Glu335 Pocket Center:")
        print(f"  center_x={center[0]:.3f}, center_y={center[1]:.3f}, center_z={center[2]:.3f}")
        return center
    return None

if __name__ == "__main__":
    print("=== PREPARING RECEPTOR STRUCTURES & BINDING POCKETS ===\n")
    
    # 1. CD36 Ectodomain (5LGD Chain A)
    clean_pdb("5LGD.pdb", "5LGD_CD36_clean.pdb", "A", "Human CD36 Ectodomain")
    extract_ligand("5LGD.pdb", "5LGD_PLM_A512.pdb", "A", "PLM", "512", "Palmitate Tunnel Ligand 1")
    extract_ligand("5LGD.pdb", "5LGD_PLM_A513.pdb", "A", "PLM", "513", "Palmitate Tunnel Ligand 2")
    find_cd36_tunnel_entrance("5LGD.pdb", "A")
    print()

    # 2. FFAR4 / GPR120 active Cryo-EM with TUG-891 (8G59 Chain R)
    clean_pdb("8G59.pdb", "8G59_FFAR4_clean.pdb", "R", "Human FFAR4 (GPR120) Receptor")
    extract_ligand("8G59.pdb", "8G59_TUG891.pdb", "R", "YN9", "401", "TUG-891 Agonist")
    print()

    # 3. FFAR4 / GPR120 active Cryo-EM with Oleic Acid (8ID6 Chain R)
    clean_pdb("8ID6.pdb", "8ID6_FFAR4_clean.pdb", "R", "Human FFAR4 (GPR120) Receptor")
    extract_ligand("8ID6.pdb", "8ID6_OLA.pdb", "R", "OLA", "401", "Oleic Acid Agonist")
    print("\nReceptor preparation and pocket definition complete.")
