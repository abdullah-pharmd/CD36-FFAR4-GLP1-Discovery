"""
Shared cheminformatics, receptor preparation, deterministic AutoDock Vina
execution, symmetry-corrected RMSD calculation, and 3D polar/hydrophobic
interaction analysis utilities shared across Protocol A and Protocol B.
"""
import os
import math
import shutil
import subprocess
import numpy as np
from scipy.optimize import linear_sum_assignment
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Crippen, QED
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

OBABEL_EXE = (
    os.environ.get("OBABEL_EXE")
    or shutil.which("obabel")
    or r"C:\Users\Abdullah\AppData\Local\Programs\Python\Python312\Scripts\obabel.exe"
)

# Exact published Daina & Zoete (ChemMedChem 2016, 11:1117-1121; PMC5089604) BOILED-Egg ellipses:
# R plotrix::draw.ellipse(x0, y0, a, b, angle_deg, deg=TRUE) in (TPSA, WLOGP) space.
BOILED_EGG_HIA = {"x0": 71.051, "y0": 2.292, "a": 142.081 / 2.0, "b": 8.740 / 2.0, "angle_deg": -1.031325}
BOILED_EGG_BBB = {"x0": 38.117, "y0": 3.177, "a": 82.061 / 2.0, "b": 5.557 / 2.0, "angle_deg": -0.171887}


def point_in_boiled_egg_ellipse(tpsa, wlogp, params):
    """Exact rotated-ellipse membership test matching Daina & Zoete (2016) ChemMedChem."""
    th = math.radians(params["angle_deg"])
    dx = tpsa - params["x0"]
    dy = wlogp - params["y0"]
    xr = dx * math.cos(th) + dy * math.sin(th)
    yr = -dx * math.sin(th) + dy * math.cos(th)
    return ((xr / params["a"]) ** 2 + (yr / params["b"]) ** 2) <= 1.0

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


def clean_pdb_strict(pdb_in, pdb_out, chain_id, het_resname=None, het_resseq=None, het_out=None, keep_altloc="A"):
    """
    Isolates a single receptor chain, removes waters, N-glycans, and co-factors,
    and retains a single crystallographic alternate-conformation state (`keep_altloc`,
    default 'A'; pass 'B' for alternate-rotamer sensitivity cross-checks).
    Optionally extracts the native co-crystallized ligand for redocking validation.
    """
    atom_lines = []
    het_lines = []
    het_coords = []

    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            alt = line[16:17]
            if alt not in (" ", keep_altloc):
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


def compute_tunnel_radial_error(path_a, path_b, is_pdbqt_a=False, is_pdbqt_b=True):
    """
    Computes 3D center-of-mass shift, XY radial shift, and perpendicular deviation
    relative to the crystal reference's principal longitudinal axis (SVD axis 0).
    """
    pa = np.array([[x, y, z] for _, x, y, z, _ in parse_heavy_atoms(path_a, is_pdbqt=is_pdbqt_a)])
    pb = np.array([[x, y, z] for _, x, y, z, _ in parse_heavy_atoms(path_b, is_pdbqt=is_pdbqt_b)])
    if len(pa) == 0 or len(pb) == 0:
        return {"axis_perp_err_A": float("nan"), "com_3d_err_A": float("nan"), "xy_radial_err_A": float("nan")}
    ca, cb = pa.mean(axis=0), pb.mean(axis=0)
    d_com = cb - ca
    axis = np.linalg.svd(pa - ca)[2][0]
    perp = np.linalg.norm(d_com - np.dot(d_com, axis) * axis)
    com_3d = np.linalg.norm(d_com)
    xy_rad = np.linalg.norm(d_com[:2])
    return {
        "axis_perp_err_A": round(float(perp), 3),
        "com_3d_err_A": round(float(com_3d), 3),
        "xy_radial_err_A": round(float(xy_rad), 3),
    }


def parse_ligand_pdbqt_with_hd(pose_pdbqt):
    """Parses Model 1 ligand heavy atoms from PDBQT along with AutoDock type and bonded HD status."""
    all_atoms = []
    in_model = False
    with open(pose_pdbqt, "r", encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("MODEL 1") or l.startswith("MODEL        1"):
                in_model = True
            elif l.startswith("ENDMDL") and in_model:
                break
            elif in_model and l.startswith(("ATOM", "HETATM")):
                aname = l[12:16].strip()
                adtype = (l[77:79].strip() if len(l) >= 78 else "").upper()
                if not adtype:
                    adtype = "".join(c for c in aname if c.isalpha())[:1].upper()
                all_atoms.append((aname, adtype, float(l[30:38]), float(l[38:46]), float(l[46:54])))

    hd_coords = [(x, y, z) for _, ad, x, y, z in all_atoms if ad == "HD"]
    heavy = []
    for aname, adtype, x, y, z in all_atoms:
        if adtype in ("H", "HD") or aname.startswith("H"):
            continue
        has_hd = any(((x - hx) ** 2 + (y - hy) ** 2 + (z - hz) ** 2) <= 1.25 ** 2 for hx, hy, hz in hd_coords)
        heavy.append((aname, adtype, x, y, z, has_hd))
    return heavy


def _classify_receptor_polar_role(rname, aname):
    """Returns (role, tag) for receptor polar N/O/S atoms: 'donor', 'acceptor', 'both', or None."""
    if aname == "O":
        return "acceptor", "bb-O"
    if aname == "N":
        return (None, None) if rname == "PRO" else ("donor", "bb-NH")
    sc_roles = {
        ("ARG", "NE"): "donor", ("ARG", "NH1"): "donor", ("ARG", "NH2"): "donor",
        ("LYS", "NZ"): "donor", ("TRP", "NE1"): "donor",
        ("ASP", "OD1"): "acceptor", ("ASP", "OD2"): "acceptor",
        ("GLU", "OE1"): "acceptor", ("GLU", "OE2"): "acceptor",
        ("ASN", "OD1"): "acceptor", ("ASN", "ND2"): "donor",
        ("GLN", "OE1"): "acceptor", ("GLN", "NE2"): "donor",
        ("SER", "OG"): "both", ("THR", "OG1"): "both", ("TYR", "OH"): "both",
        ("HIS", "ND1"): "both", ("HIS", "NE2"): "both", ("CYS", "SG"): "both",
    }
    role = sc_roles.get((rname, aname))
    return (role, "sc") if role else (None, None)


def analyze_pocket_interactions_3d(clean_pdb, pose_pdbqt):
    """
    Classifies Model 1 protein-ligand contacts from 3D coordinates with chemical donor/acceptor validation:
      - Chemically valid H-bonds (<= 3.5 Å; verifying ligand HD donor & receptor acceptor, or ligand OA/NA/SA acceptor & receptor donor, plus backbone vs side-chain tags)
      - Additional non-H-bonded heavy-atom pocket contacts (<= 4.0 Å; strictly excluding residues already reported as H-bonds)
      - Direct steric clashes (< 2.2 Å between non-bonded heavy atoms)
    """
    lig_atoms = parse_ligand_pdbqt_with_hd(pose_pdbqt)
    if not lig_atoms:
        return {"polar_hbonds": "None", "all_contacts_4A": "None", "steric_clashes": "0"}

    res_min_all = {}
    res_min_hbond = {}
    clash_list = []

    with open(clean_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for l in f:
            if not l.startswith("ATOM"):
                continue
            aname = l[12:16].strip()
            if aname.startswith("H"):
                continue
            rname = l[17:20].strip()
            rseq = int(l[22:26].strip())
            rkey = f"{rname}{rseq}"
            rx, ry, rz = float(l[30:38]), float(l[38:46]), float(l[46:54])
            rec_role, rec_tag = _classify_receptor_polar_role(rname, aname)

            for _, adtype, lx, ly, lz, has_hd in lig_atoms:
                d = math.sqrt((lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2)
                if d <= 4.0 and (rkey not in res_min_all or d < res_min_all[rkey]):
                    res_min_all[rkey] = round(d, 2)
                if d <= 3.5 and rec_role is not None:
                    lig_is_donor = has_hd and adtype in ("N", "NA", "O", "OA", "S", "SA")
                    lig_is_acceptor = adtype in ("NA", "OA", "SA")
                    valid_hb = (
                        (lig_is_donor and rec_role in ("acceptor", "both"))
                        or (lig_is_acceptor and rec_role in ("donor", "both"))
                    )
                    if valid_hb:
                        d_r = round(d, 2)
                        if rkey not in res_min_hbond or d_r < res_min_hbond[rkey][0]:
                            res_min_hbond[rkey] = (d_r, rec_tag)
                if d < 2.2:
                    clash_list.append((rkey, round(d, 2)))

    sorted_hb = sorted(res_min_hbond.items(), key=lambda x: x[1][0])[:5]
    hb_resids = {k for k, _ in sorted_hb}
    non_hb_contacts = [(k, v) for k, v in sorted(res_min_all.items(), key=lambda x: x[1]) if k not in hb_resids][:7]
    top_all = ", ".join(f"{k} ({v:.2f} Å)" for k, v in non_hb_contacts) or "None additional"
    top_hb = ", ".join(f"{k} ({v[0]:.2f} Å, {v[1]})" for k, v in sorted_hb) or "Hydrophobic/π-stacking only"
    clash_summary = f"{len(clash_list)} (<2.2 Å: {', '.join(f'{r} {d:.2f}Å' for r, d in sorted(clash_list, key=lambda x: x[1])[:3])})" if clash_list else "0"

    return {
        "polar_hbonds": top_hb,
        "all_contacts_4A": top_all,
        "steric_clashes": clash_summary
    }


def extract_isomeric_smiles(sdf_path, fallback_smiles=""):
    """Reads 3D SDF and returns canonical isomeric SMILES preserving chiral @/@@ centers."""
    try:
        sup = Chem.SDMolSupplier(sdf_path, removeHs=True)
        mol = sup[0] if sup and len(sup) > 0 else None
        if mol is not None:
            return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        pass
    return fallback_smiles


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

    # Exact Daina & Zoete (2016) BOILED-Egg Elliptical GI & BBB Permeability Classification
    has_quat_n = any(a.GetSymbol() == "N" and a.GetFormalCharge() > 0 and a.GetTotalNumHs() == 0 for a in mol.GetAtoms())
    if has_quat_n:
        hia_egg = "Low (Quat. Cation / P-gp)"
        bbb_egg = "BBB- (Quat. Cation)"
    else:
        in_hia = point_in_boiled_egg_ellipse(tpsa, clogp, BOILED_EGG_HIA)
        in_bbb = point_in_boiled_egg_ellipse(tpsa, clogp, BOILED_EGG_BBB)
        hia_egg = "High (Inside White Egg)" if in_hia else "Low (Outside Egg)"
        bbb_egg = "BBB+ (Yolk)" if in_bbb else "BBB- (Non-Permeant)"

    # Delaney (2004) ESOL Aqueous Solubility LogS (mol/L -> g/L=mg/mL and ug/mL)
    ar_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    ap = ar_atoms / n_heavy if n_heavy > 0 else 0.0
    esol_logs = round(0.16 - 0.63 * clogp - 0.0062 * mw + 0.066 * rotb - 0.74 * ap, 2)
    sol_molar = 10.0 ** esol_logs
    sol_mg_ml = round(sol_molar * mw, 4)          # 1 g/L == 1 mg/mL
    sol_ug_ml = round(sol_molar * mw * 1000.0, 4) # 1 mg/L == 1 ug/mL
    assert abs(sol_ug_ml - sol_mg_ml * 1000.0) <= 0.15, f"Solubility unit mismatch: {sol_mg_ml} mg/mL vs {sol_ug_ml} ug/mL"

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
        "Estimated_Solubility_ug_mL": sol_ug_ml,
        "PAINS_Structural_Alert": pains_flag,
        "Brenk_Structural_Alert": brenk_flag
    }
