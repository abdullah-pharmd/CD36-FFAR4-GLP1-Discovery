"""
07_audit_8g59_phe115_and_fix_admet.py
Protocol A (FFAR4 8G59 + CD36 5LGD, salivary pH 6.8) integrity repair & structural sensitivity audit:
  1. Repairs the 1,000x solubility unit bug across all 22 rows of
     results/admet_and_physicochemical_properties.csv by recomputing from exact molar ESOL LogS
     into separate `Estimated_Solubility_mg_mL` and `Estimated_Solubility_ug_mL` columns
     (with strict unit-conversion assertion), preserving isomeric SMILES + full InChIKey from SDFs,
     and updating BOILED-Egg classification to exact published Daina & Zoete (2016) ellipses.
  2. Updates results/molecular_docking_results.csv to retire pseudo-thermodynamic `pKd`/`LLE`
     labels in favor of `aLLE` (docking-derived apparent lipophilic efficiency) and remove the
     unselective `"Efficiency-Validated"` tier wording (given 0% specificity at LE >= 0.30 on
     the CHEMBL5339 panel).
  3. Runs a complete 8G59 Phe115 alternate-conformation (altLoc A vs altLoc B, 0.50/0.50 occupancy)
     sensitivity audit:
       - Fixed-pose closest heavy-atom distances to Phe115 A vs B AND rest of receptor (kept distinct!)
       - Fresh AutoDock Vina v1.2.7 redocking (exhaustiveness=8, seed=42) into 8G59 Conformer B
         across all 4 triterpenoids (Gymnemagenin, Oleanolic acid, Ursolic acid, Corosolic acid)
         and key Tier 1 leads/controls (TUG-891, Naringenin, Phloretin, Resveratrol, Spilanthol).
     Saves results/audit_8g59_phe115_conformer_ab.csv.
"""
import os
import sys
import glob
import numpy as np
import pandas as pd
from rdkit import Chem

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))

from structural_prep_utils import (
    clean_pdb_strict, prepare_receptor_pdbqt, shift_and_dock_vina,
    hungarian_heavy_atom_rmsd, compute_admet_descriptors
)

STRUCT_DIR = os.path.join(BASE_DIR, "structures")
LIG_DIR = os.path.join(BASE_DIR, "ligands")
SDF_DIR = os.path.join(LIG_DIR, "sdf")
PDBQT_DIR = os.path.join(LIG_DIR, "pdbqt_pH68")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
POSES_FFAR4 = os.path.join(RESULTS_DIR, "poses", "FFAR4")
VINA_EXE = os.path.join(BASE_DIR, "scripts", "vina.exe")


def parse_pdb_atoms_excl_h(pdb_path, chain, alt_keep):
    """Returns list of (resname, resseq, coords_xyz) for heavy atoms of a chain retaining alt_keep."""
    res_atoms = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("ATOM  ") and l[21] == chain and l[16] in (" ", alt_keep):
                el = l[76:78].strip()
                if el != "H":
                    res_atoms.append((l[17:20].strip(), int(l[22:26].strip()), [float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    return res_atoms


def parse_pdbqt_model1_heavy(pdbqt_path):
    pts = []
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("ENDMDL"):
                break
            if l.startswith(("ATOM  ", "HETATM")):
                at = l[77:79].strip()
                if not at.startswith("H"):
                    pts.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(pts, dtype=float)


def min_dist_to_residue(lig_pts, res_atoms, target_resseq=None):
    if target_resseq is not None:
        pts = np.array([x[2] for x in res_atoms if x[1] == target_resseq], dtype=float)
        if not len(pts):
            return np.nan, ""
        d = np.min(np.linalg.norm(lig_pts[:, None, :] - pts[None, :, :], axis=-1))
        return float(d), f"PHE{target_resseq}"
    # Global minimum + list of all residues < 2.2 A
    clashes = []
    min_d = 999.0
    min_lbl = ""
    by_res = {}
    for rn, rs, xyz in res_atoms:
        by_res.setdefault((rn, rs), []).append(xyz)
    for (rn, rs), xyzs in by_res.items():
        pts = np.array(xyzs, dtype=float)
        d = float(np.min(np.linalg.norm(lig_pts[:, None, :] - pts[None, :, :], axis=-1)))
        if d < min_d:
            min_d = d
            min_lbl = f"{rn}{rs}"
        if d < 2.2:
            clashes.append(f"{rn}{rs} ({d:.3f} Å)")
    return min_d, ("; ".join(sorted(clashes)) if clashes else "None")


def repair_protocol_a_admet_and_docking_csv():
    admet_csv = os.path.join(RESULTS_DIR, "admet_and_physicochemical_properties.csv")
    df = pd.read_csv(admet_csv)
    new_rows = []
    for _, row in df.iterrows():
        r = dict(row)
        sdf_name = str(r.get("SDF_File", ""))
        sdf_path = os.path.join(SDF_DIR, sdf_name)
        if not os.path.exists(sdf_path):
            sdf_path = os.path.join(LIG_DIR, sdf_name)
        mol = None
        if os.path.exists(sdf_path):
            suppl = Chem.SDMolSupplier(sdf_path, removeHs=True)
            if suppl and len(suppl) > 0 and suppl[0] is not None:
                mol = suppl[0]
        if mol is not None:
            iso_smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            ik_full = Chem.MolToInchiKey(mol)
        else:
            iso_smi = str(r["SMILES"])
            m2 = Chem.MolFromSmiles(iso_smi)
            ik_full = Chem.MolToInchiKey(m2) if m2 is not None else ""

        ad = compute_admet_descriptors(iso_smi)
        r["SMILES"] = iso_smi
        r["InChIKey_Full"] = ik_full
        r["BOILED_Egg_GI_Absorption"] = ad["BOILED_Egg_GI_Absorption"]
        r["BOILED_Egg_BBB"] = ad["BOILED_Egg_BBB"]
        r["ESOL_LogS_mol_L"] = ad["ESOL_LogS_mol_L"]
        r["Estimated_Solubility_mg_mL"] = ad["Estimated_Solubility_mg_mL"]
        r["Estimated_Solubility_ug_mL"] = ad["Estimated_Solubility_ug_mL"]
        assert abs(r["Estimated_Solubility_ug_mL"] - r["Estimated_Solubility_mg_mL"] * 1000.0) <= 0.15
        new_rows.append(r)

    df_fixed = pd.DataFrame(new_rows)
    # Place Estimated_Solubility_ug_mL right after Estimated_Solubility_mg_mL and InChIKey_Full after SMILES
    cols = list(df_fixed.columns)
    for col_move, after_col in [("InChIKey_Full", "SMILES"), ("Estimated_Solubility_ug_mL", "Estimated_Solubility_mg_mL")]:
        cols.remove(col_move)
        cols.insert(cols.index(after_col) + 1, col_move)
    df_fixed = df_fixed[cols]
    df_fixed.to_csv(admet_csv, index=False)

    nar = df_fixed[df_fixed["Compound"] == "Naringenin"].iloc[0]
    print(f"[PROTO A ADMET FIXED] Naringenin: LogS={nar['ESOL_LogS_mol_L']}, mg/mL={nar['Estimated_Solubility_mg_mL']}, ug/mL={nar['Estimated_Solubility_ug_mL']}, SMILES={nar['SMILES']}")

    # Sync molecular_docking_results.csv: retire pseudo-pKd wording -> apparent aLLE & honest tier names
    dock_csv = os.path.join(RESULTS_DIR, "molecular_docking_results.csv")
    if os.path.exists(dock_csv):
        ddf = pd.read_csv(dock_csv)
        ad_map = df_fixed.set_index("Compound")
        for idx, r in ddf.iterrows():
            nm = r["Compound"]
            if nm in ad_map.index:
                ddf.at[idx, "BOILED_Egg_GI"] = ad_map.loc[nm, "BOILED_Egg_GI_Absorption"]
            tier_old = str(r["Lead_Qualification_Tier"])
            ddf.at[idx, "Lead_Qualification_Tier"] = tier_old.replace("Efficiency-Validated ", "PAINS-Free Prioritized ").replace("PAINS-Free PAINS-Free", "PAINS-Free")
        # Rename pKd/LLE columns transparently to Apparent_Docking_pKd / Apparent_aLLE
        rename_map = {
            "FFAR4_pKd": "FFAR4_Apparent_pKd_Heuristic",
            "FFAR4_LLE": "FFAR4_aLLE_Docking",
            "CD36_pKd": "CD36_Apparent_pKd_Heuristic",
            "CD36_LLE": "CD36_aLLE_Docking",
        }
        ddf = ddf.rename(columns={k: v for k, v in rename_map.items() if k in ddf.columns})
        ddf.to_csv(dock_csv, index=False)
        print(f"[PROTO A DOCKING CSV UPDATED] Transparent aLLE + BOILED-Egg ellipse + honest tier labels saved.")
    return df_fixed


def audit_8g59_phe115_ab(df_admet):
    pdb_raw = os.path.join(STRUCT_DIR, "8G59.pdb")
    clean_b = os.path.join(STRUCT_DIR, "8G59_FFAR4_clean_altB.pdb")
    rec_a = os.path.join(STRUCT_DIR, "8G59_FFAR4_clean.pdbqt")
    rec_b = os.path.join(STRUCT_DIR, "8G59_FFAR4_pH68_altB.pdbqt")

    n_b, c_tug = clean_pdb_strict(pdb_raw, clean_b, chain_id="R", het_resname="YN9", het_resseq="401", keep_altloc="B")
    ok_b = prepare_receptor_pdbqt(clean_b, rec_b, ph=6.8)
    assert ok_b, "Failed to build 8G59 altLoc B receptor PDBQT"

    res_a = parse_pdb_atoms_excl_h(pdb_raw, "R", "A")
    res_b = parse_pdb_atoms_excl_h(pdb_raw, "R", "B")

    box_8g59 = {"center_x": float(c_tug[0]), "center_y": float(c_tug[1]), "center_z": float(c_tug[2]), "size_x": 22.0, "size_y": 22.0, "size_z": 22.0}
    out_dir = os.path.join(RESULTS_DIR, "poses", "FFAR4_altB")
    os.makedirs(out_dir, exist_ok=True)

    panel = [
        "TUG-891", "Naringenin", "Resveratrol", "Phloretin", "Spilanthol",
        "Gymnemagenin", "Oleanolic acid", "Ursolic acid", "Corosolic acid"
    ]
    rows = []
    for name in panel:
        r = df_admet[df_admet["Compound"] == name].iloc[0]
        cid = int(r["PubChem_CID"])
        safe = name.replace(" ", "_").replace("-", "_")
        # Locate existing pose in Conformer A
        pose_a_cands = glob.glob(os.path.join(POSES_FFAR4, f"{cid}_*_FFAR4_out.pdbqt"))
        pose_a = pose_a_cands[0] if pose_a_cands else ""
        # Locate ligand pH 6.8 pdbqt
        lig_cands = glob.glob(os.path.join(PDBQT_DIR, f"{cid}_*.pdbqt")) + glob.glob(os.path.join(LIG_DIR, f"{cid}_*.pdbqt"))
        lig_qt = lig_cands[0] if lig_cands else ""

        # Fixed-pose analysis of existing Conformer A pose against A vs B
        dg_a_saved = np.nan
        fp_phe_a = fp_phe_b = np.nan
        fp_clashes_a = fp_clashes_b = ""
        if pose_a and os.path.exists(pose_a):
            with open(pose_a, "r", encoding="utf-8", errors="ignore") as fp:
                for line in fp:
                    if "REMARK VINA RESULT:" in line:
                        dg_a_saved = float(line.split()[3])
                        break
            pts_a = parse_pdbqt_model1_heavy(pose_a)
            fp_phe_a, _ = min_dist_to_residue(pts_a, res_a, 115)
            fp_phe_b, _ = min_dist_to_residue(pts_a, res_b, 115)
            _, fp_clashes_a = min_dist_to_residue(pts_a, res_a, None)
            _, fp_clashes_b = min_dist_to_residue(pts_a, res_b, None)

        # Fresh redocking into 8G59 Conformer B (seed=42, exhaustiveness=8)
        out_b = os.path.join(out_dir, f"{cid}_{safe}_FFAR4_altB_out.pdbqt")
        dg_b = np.nan
        rd_phe_b = np.nan
        rd_clashes_b = ""
        rmsd_ab = np.nan
        if lig_qt and os.path.exists(lig_qt):
            modes_b = shift_and_dock_vina(VINA_EXE, rec_b, lig_qt, out_b, box_8g59, exhaustiveness=8, seed=42, cpu=4)
            dg_b = modes_b[0]["affinity"] if modes_b else np.nan
            pts_b = parse_pdbqt_model1_heavy(out_b)
            rd_phe_b, _ = min_dist_to_residue(pts_b, res_b, 115)
            _, rd_clashes_b = min_dist_to_residue(pts_b, res_b, None)
            if pose_a and os.path.exists(pose_a):
                rmsd_ab = hungarian_heavy_atom_rmsd(pose_a, out_b, is_pdbqt_a=True, is_pdbqt_b=True)

        rows.append({
            "Compound": name,
            "ConformerA_dG_kcal_mol": round(dg_a_saved, 3),
            "FixedPose_Phe115_AltA_A": round(fp_phe_a, 3),
            "FixedPose_Phe115_AltB_A": round(fp_phe_b, 3),
            "FixedPose_AllClashes_AltA": fp_clashes_a,
            "FixedPose_AllClashes_AltB": fp_clashes_b,
            "Redocked_ConformerB_dG_kcal_mol": round(dg_b, 3),
            "Delta_dG_AltB_minus_AltA": round(dg_b - dg_a_saved, 3) if not np.isnan(dg_b) else np.nan,
            "Redocked_AltB_Phe115_Dist_A": round(rd_phe_b, 3),
            "Redocked_AltB_AllClashes": rd_clashes_b,
            "Pose_RMSD_AltA_vs_AltB_A": round(rmsd_ab, 3) if not np.isnan(rmsd_ab) else np.nan,
        })
        print(f"  {name:16s} | AltA dG={dg_a_saved:6.2f} (Phe115 fixed A={fp_phe_a:.2f}/B={fp_phe_b:.2f} Å) | Redocked AltB dG={dg_b:6.2f} (Δ={dg_b-dg_a_saved:+.2f}, Phe115={rd_phe_b:.2f} Å, ClashesB={rd_clashes_b})", flush=True)

    df_ab = pd.DataFrame(rows)
    out_csv = os.path.join(RESULTS_DIR, "audit_8g59_phe115_conformer_ab.csv")
    df_ab.to_csv(out_csv, index=False)
    print(f"[PROTO A 8G59 CONFORMER A/B AUDIT SAVED] {out_csv}")
    return df_ab


if __name__ == "__main__":
    df_fixed = repair_protocol_a_admet_and_docking_csv()
    audit_8g59_phe115_ab(df_fixed)
