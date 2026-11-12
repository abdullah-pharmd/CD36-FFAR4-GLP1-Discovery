"""
04_extract_atomic_contacts.py
Extracts exact heavy-atom interatomic distances (<= 4.0 A) between top-ranked
AutoDock Vina ligand poses (Model 1) and receptor binding-pocket residues in
human FFAR4 (PDB: 8G59) and human CD36 (PDB: 5LGD).
"""
import os
import math
import csv

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FFAR4_PDB = os.path.join(BASE_DIR, "structures", "8G59_FFAR4_clean.pdb")
CD36_PDB = os.path.join(BASE_DIR, "structures", "5LGD_CD36_clean.pdb")
OUT_CSV = os.path.join(BASE_DIR, "results", "docked_pose_residue_contacts.csv")


def parse_pdbqt_model1(pdbqt_path):
    atoms = []
    in_model = False
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("MODEL 1") or line.startswith("MODEL        1"):
                in_model = True
            elif line.startswith("ENDMDL") and in_model:
                break
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                if not atom_name.startswith("H"):
                    atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54]), atom_name))
    return atoms


def get_nearby_residues(pdb_path, lig_atoms, cutoff=4.0):
    res_dist = {}
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name.startswith("H"):
                    continue
                res_name = line[17:20].strip()
                res_seq = int(line[22:26].strip())
                key = f"{res_name}{res_seq}"
                rx, ry, rz = float(line[30:38]), float(line[38:46]), float(line[46:54])
                d = min(math.sqrt((lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2) for lx, ly, lz, _ in lig_atoms)
                if d <= cutoff:
                    if key not in res_dist or d < res_dist[key]:
                        res_dist[key] = round(d, 2)
    return dict(sorted(res_dist.items(), key=lambda x: x[1]))


def main():
    rows = []
    ffar4_panel = [
        ("57522038_TUG_891_FFAR4.pdbqt", "TUG-891 (Control)", "FFAR4 (8G59)", -11.73),
        ("638024_Piperine_FFAR4.pdbqt", "Piperine", "FFAR4 (8G59)", -10.34),
        ("1794427_Chlorogenic_acid_FFAR4.pdbqt", "Chlorogenic acid", "FFAR4 (8G59)", -9.62),
        ("439246_Naringenin_FFAR4.pdbqt", "Naringenin", "FFAR4 (8G59)", -9.53),
        ("5280343_Quercetin_FFAR4.pdbqt", "Quercetin", "FFAR4 (8G59)", -9.30),
        ("5318516_Pellitorine_FFAR4.pdbqt", "Pellitorine", "FFAR4 (8G59)", -7.62),
    ]
    for fn, label, target, score in ffar4_panel:
        pose_path = os.path.join(BASE_DIR, "results", "poses", "FFAR4", fn)
        atoms = parse_pdbqt_model1(pose_path)
        contacts = get_nearby_residues(FFAR4_PDB, atoms, cutoff=4.0)
        contact_str = "; ".join(f"{k} ({v:.2f} A)" for k, v in list(contacts.items())[:7])
        rows.append({
            "Ligand": label,
            "Receptor": target,
            "Vina_Affinity_kcal_mol": score,
            "Closest_Residue_Contacts_le_4A": contact_str,
        })

    cd36_panel = [
        ("46782948_Sulfo_N_succinimidyl_oleate_CD36.pdbqt", "SSO (Control)", "CD36 (5LGD)", -5.80),
        ("86656470_Nitro_oleic_acid_CD36.pdbqt", "Nitro-oleic acid (Control)", "CD36 (5LGD)", -6.51),
        ("10114_Glycyrrhetinic_acid_CD36.pdbqt", "Glycyrrhetinic acid", "CD36 (5LGD)", -6.51),
        ("6918774_Corosolic_acid_CD36.pdbqt", "Corosolic acid", "CD36 (5LGD)", -6.34),
        ("439246_Naringenin_CD36.pdbqt", "Naringenin", "CD36 (5LGD)", -6.15),
        ("638024_Piperine_CD36.pdbqt", "Piperine", "CD36 (5LGD)", -5.46),
    ]
    for fn, label, target, score in cd36_panel:
        pose_path = os.path.join(BASE_DIR, "results", "poses", "CD36", fn)
        atoms = parse_pdbqt_model1(pose_path)
        contacts = get_nearby_residues(CD36_PDB, atoms, cutoff=4.0)
        contact_str = "; ".join(f"{k} ({v:.2f} A)" for k, v in list(contacts.items())[:7])
        rows.append({
            "Ligand": label,
            "Receptor": target,
            "Vina_Affinity_kcal_mol": score,
            "Closest_Residue_Contacts_le_4A": contact_str,
        })

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Ligand", "Receptor", "Vina_Affinity_kcal_mol", "Closest_Residue_Contacts_le_4A"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"[SUCCESS] Wrote {len(rows)} real structural contact profiles to {OUT_CSV}")


if __name__ == "__main__":
    main()
