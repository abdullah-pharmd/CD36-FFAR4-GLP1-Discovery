"""
Shared cheminformatics, receptor preparation, deterministic AutoDock Vina
execution, symmetry-corrected RMSD calculation, and 3D polar/hydrophobic
interaction analysis utilities shared across Protocol A and Protocol B.
"""
import os
import math
import subprocess
import numpy as np
from scipy.optimize import linear_sum_assignment
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Crippen, QED
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

OBABEL_EXE = r"C:\Users\Abdullah\AppData\Local\Programs\Python\Python312\Scripts\obabel.exe"

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"
}

_pains_params = FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
PAINS_CATALOG = FilterCatalog(_pains_params)

_brenk_params = FilterCatalogParams()
_brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
BRENK_CATALOG = FilterCatalog(_brenk_params)


def clean_pdb_strict(pdb_in, pdb_out, chain_id, het_resname=None, het_resseq=None, het_out=None):
    """
    Isolates a single receptor chain, removes waters, N-glycans, and co-factors,
    and strictly excludes duplicate alternative conformations (altLoc != ' ' and 'A').
    Optionally extracts the native co-crystallized ligand for redocking validation.
    """
    atom_lines = []
    het_lines = []
    het_coords = []

    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            alt = line[16:17]
            if alt not in (" ", "A"):
                continue
            clean_line = line[:16] + " " + line[17:]
            if line.startswith("ATOM  "):
                ch = line[21:22].strip()
                res_name = line[17:20].strip()
                if ch == chain_id and res_name in STANDARD_AA:
                    atom_lines.append(clean_line)
            elif het_resname and line.startswith("HETATM"):
                ch = line[21:22].strip()
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                if ch == chain_id and res_name == het_resname and (het_resseq is None or res_seq == str(het_resseq)):
                    het_lines.append(clean_line)
                    het_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

    with open(pdb_out, "w", encoding="utf-8") as f:
        f.writelines(atom_lines)
        f.write("END\n")

    if het_out and het_lines:
        with open(het_out, "w", encoding="utf-8") as f:
            f.writelines(het_lines)
            f.write("END\n")

    center = np.mean(het_coords, axis=0) if het_coords else np.zeros(3)
    return len(atom_lines), center


def prepare_receptor_pdbqt(clean_pdb, rec_pdbqt, ph=7.4):
    """Converts cleaned PDB chain to rigid PDBQT with physiological polar HD donors."""
    cmd = [OBABEL_EXE, clean_pdb, "-O", rec_pdbqt, "-p", str(ph), "-xr"]
    subprocess.run(cmd, capture_output=True, text=True)
    return os.path.exists(rec_pdbqt) and os.path.getsize(rec_pdbqt) > 100


def prepare_ligand_pdbqt(sdf_or_pdb_in, lig_pdbqt, ph=7.4):
    """
    Re-protonates ligand 3D coordinates at target physiological pH (-d -p <ph>)
    so ionizable carboxylic acids (-COO-) and amines (-NH3+) carry accurate formal states.
    """
    cmd = [OBABEL_EXE, sdf_or_pdb_in, "-O", lig_pdbqt, "-d", "-p", str(ph)]
    subprocess.run(cmd, capture_output=True, text=True)
    return os.path.exists(lig_pdbqt) and os.path.getsize(lig_pdbqt) > 100


def shift_and_dock_vina(vina_exe, rec_pdbqt, lig_pdbqt, out_pdbqt, box, exhaustiveness=8, seed=42, cpu=4):
    """Centers ligand centroid inside docking box and runs deterministic AutoDock Vina."""
    lines, coords = [], []
    with open(lig_pdbqt, "r", encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith(("ATOM", "HETATM")):
                coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            lines.append(l)
    if not coords:
        return []
    shift = np.array([box["center_x"], box["center_y"], box["center_z"]]) - np.mean(coords, axis=0)
    prep = out_pdbqt.replace(".pdbqt", "_prep.pdbqt")
    with open(prep, "w", encoding="utf-8") as f:
        for l in lines:
            if l.startswith(("ATOM", "HETATM")):
                x = float(l[30:38]) + shift[0]
                y = float(l[38:46]) + shift[1]
                z = float(l[46:54]) + shift[2]
                f.write(l[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + l[54:])
            else:
                f.write(l)

    cmd = [
        vina_exe,
        "--receptor", rec_pdbqt,
        "--ligand", prep,
        "--center_x", f"{box['center_x']:.2f}",
        "--center_y", f"{box['center_y']:.2f}",
        "--center_z", f"{box['center_z']:.2f}",
        "--size_x", f"{box['size_x']:.2f}",
        "--size_y", f"{box['size_y']:.2f}",
        "--size_z", f"{box['size_z']:.2f}",
        "--out", out_pdbqt,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "9",
        "--seed", str(seed),
        "--cpu", str(cpu)
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if os.path.exists(prep):
        os.remove(prep)

    modes = []
    in_table = False
    for line in res.stdout.splitlines():
        if "-----+------------+----------+----------" in line:
            in_table = True
            continue
        if in_table:
            parts = line.strip().split()
            if len(parts) >= 4 and parts[0].isdigit():
                modes.append({
                    "mode": int(parts[0]),
                    "affinity": float(parts[1]),
                    "rmsd_lb": float(parts[2]),
                    "rmsd_ub": float(parts[3])
                })
            elif not parts:
                break
    return modes


def parse_heavy_atoms(path, is_pdbqt=False):
    """Parses heavy atoms (element, x, y, z, atom_name) from PDB or Model 1 of PDBQT."""
    atoms = []
    in_model = not is_pdbqt
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for l in f:
            if is_pdbqt:
                if l.startswith("MODEL 1") or l.startswith("MODEL        1"):
                    in_model = True
                elif l.startswith("ENDMDL") and in_model:
                    break
            if in_model and l.startswith(("ATOM", "HETATM")):
                aname = l[12:16].strip()
                raw_el = l[76:78].strip() if len(l) >= 78 else ""
                if not raw_el:
                    raw_el = "".join(c for c in aname if c.isalpha())[:1]
                el_map = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H", "CL": "CL", "BR": "BR", "F": "F"}
                el = el_map.get(raw_el.upper(), raw_el.upper()[:1] if raw_el.upper() not in ("CL", "BR") else raw_el.upper())
                if el == "H" or aname.startswith("H"):
                    continue
                atoms.append((el, float(l[30:38]), float(l[38:46]), float(l[46:54]), aname))
    return atoms


def hungarian_heavy_atom_rmsd(path_a, path_b, is_pdbqt_a=False, is_pdbqt_b=True):
    """
    Computes exact symmetry-corrected heavy-atom RMSD (Å) via optimal bipartite
    Hungarian assignment (scipy.optimize.linear_sum_assignment) per chemical element.
    """
    atoms_a = parse_heavy_atoms(path_a, is_pdbqt=is_pdbqt_a)
    atoms_b = parse_heavy_atoms(path_b, is_pdbqt=is_pdbqt_b)
    elements = set(e for e, *_ in atoms_a)
    sq_dists = []
    for el in elements:
        pa = np.array([[x, y, z] for e, x, y, z, _ in atoms_a if e == el])
        pb = np.array([[x, y, z] for e, x, y, z, _ in atoms_b if e == el])
        if len(pa) == 0 or len(pa) != len(pb):
            return float("nan")
        cost = np.sum((pa[:, None, :] - pb[None, :, :]) ** 2, axis=-1)
        ri, ci = linear_sum_assignment(cost)
        sq_dists.extend(cost[ri, ci].tolist())
    return round(math.sqrt(sum(sq_dists) / len(sq_dists)), 3)


def analyze_pocket_interactions_3d(clean_pdb, pose_pdbqt):
    """
    Classifies Model 1 protein-ligand contacts from 3D coordinates into:
      - Polar H-bond / Electrostatic interactions (Ligand N/O <-> Receptor N/O/S <= 3.5 Å)
      - Closest heavy-atom pocket contacts (<= 4.0 Å)
      - Direct steric clashes (< 2.2 Å between non-bonded heavy atoms)
    """
    lig_atoms = parse_heavy_atoms(pose_pdbqt, is_pdbqt=True)
    if not lig_atoms:
        return {"polar_hbonds": "None", "all_contacts_4A": "None", "steric_clashes": "0"}

    polar_els = {"N", "O", "S"}
    res_min_all = {}
    res_min_polar = {}
    clash_list = []

    with open(clean_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for l in f:
            if not l.startswith("ATOM"):
                continue
            aname = l[12:16].strip()
            if aname.startswith("H"):
                continue
            rel = "".join(c for c in aname if c.isalpha())[:1].upper()
            rname = l[17:20].strip()
            rseq = int(l[22:26].strip())
            rkey = f"{rname}{rseq}"
            rx, ry, rz = float(l[30:38]), float(l[38:46]), float(l[46:54])

            for lel, lx, ly, lz, _ in lig_atoms:
                d = math.sqrt((lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2)
                if d <= 4.0 and (rkey not in res_min_all or d < res_min_all[rkey]):
                    res_min_all[rkey] = round(d, 2)
                if d <= 3.5 and (lel in polar_els) and (rel in polar_els):
                    if rkey not in res_min_polar or d < res_min_polar[rkey]:
                        res_min_polar[rkey] = round(d, 2)
                if d < 2.2:
                    clash_list.append((rkey, round(d, 2)))

    top_all = ", ".join(f"{k} ({v:.2f} Å)" for k, v in sorted(res_min_all.items(), key=lambda x: x[1])[:7])
    top_hb = ", ".join(f"{k} ({v:.2f} Å)" for k, v in sorted(res_min_polar.items(), key=lambda x: x[1])[:5]) or "Hydrophobic/π-only"
    clash_summary = f"{len(clash_list)} (<2.2 Å: {', '.join(f'{r} {d:.2f}Å' for r, d in sorted(clash_list, key=lambda x: x[1])[:3])})" if clash_list else "0"

    return {
        "polar_hbonds": top_hb,
        "all_contacts_4A": top_all,
        "steric_clashes": clash_summary
    }


def compute_admet_descriptors(smiles):
    """Computes RDKit ADMET, Delaney ESOL, Egan/BOILED-Egg GI/BBB, PAINS & Brenk alerts."""
    mol = Chem.MolFromSmiles(smiles)
    mw = round(Descriptors.MolWt(mol), 2)
    n_heavy = int(mol.GetNumHeavyAtoms())
    clogp = round(Crippen.MolLogP(mol), 2)
    tpsa = round(rdMolDescriptors.CalcTPSA(mol), 2)
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    rotb = int(Lipinski.NumRotatableBonds(mol))
    fsp3 = round(rdMolDescriptors.CalcFractionCSP3(mol), 2)
    qed = round(QED.qed(mol), 3)

    ro5_v = sum([mw > 500, clogp > 5, hbd > 5, hba > 10])
    lipinski_pass = "Yes" if ro5_v <= 1 else "No"
    veber_pass = "Yes" if (rotb <= 10 and tpsa <= 140) else "No"

    # Egan / BOILED-Egg Passive GI & BBB Permeability Windows (Daina & Zoete 2016)
    hia_egg = "High (Inside White Egg)" if (tpsa <= 131.6 and -0.7 <= clogp <= 5.88) else "Low (Outside Egg)"
    bbb_egg = "BBB+ (Yolk)" if (tpsa < 79.0 and 0.4 <= clogp <= 6.0) else "BBB- (Non-Permeant)"

    # Delaney (2004) ESOL Aqueous Solubility LogS
    ar_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    ap = ar_atoms / n_heavy if n_heavy > 0 else 0.0
    esol_logs = round(0.16 - 0.63 * clogp - 0.0062 * mw + 0.066 * rotb - 0.74 * ap, 2)
    sol_mg_ml = round((10 ** esol_logs) * mw * 1000.0, 4)

    pains_m = PAINS_CATALOG.GetMatches(mol)
    pains_flag = "0 Alerts" if not pains_m else f"1 Alert ({pains_m[0].GetDescription()})"

    brenk_m = BRENK_CATALOG.GetMatches(mol)
    brenk_flag = "0 Alerts" if not brenk_m else f"{len(brenk_m)} ({brenk_m[0].GetDescription()})"

    return {
        "MW_Da": mw,
        "Heavy_Atoms": n_heavy,
        "LogP_WildmanCrippen": clogp,
        "TPSA_A2": tpsa,
        "HBD": hbd,
        "HBA": hba,
        "Rotatable_Bonds": rotb,
        "Fraction_CSP3": fsp3,
        "QED_Score": qed,
        "Lipinski_Violations": ro5_v,
        "Lipinski_Pass": lipinski_pass,
        "Veber_Pass": veber_pass,
        "BOILED_Egg_GI_Absorption": hia_egg,
        "BOILED_Egg_BBB": bbb_egg,
        "ESOL_LogS_mol_L": esol_logs,
        "Estimated_Solubility_mg_mL": sol_mg_ml,
        "PAINS_Structural_Alert": pains_flag,
        "Brenk_Structural_Alert": brenk_flag
    }
