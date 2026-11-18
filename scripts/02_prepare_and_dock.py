"""
02_prepare_and_dock.py
Protocol A: altLoc-clean Chain Surgery + Physiological Salivary pH 6.8 Receptor/Ligand
Protonation (HD polar donors enabled) + Uniform Deterministic AutoDock Vina (--seed 42,
exhaustiveness = 8) + Multi-Exhaustiveness (exh=12) Convergence RMSD Verification.
"""
import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structural_prep_utils import (
    clean_pdb_strict, prepare_receptor_pdbqt, prepare_ligand_pdbqt,
    shift_and_dock_vina, hungarian_heavy_atom_rmsd
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT_DIR = os.path.join(BASE_DIR, "structures")
LIG_DIR = os.path.join(BASE_DIR, "ligands")
SDF_DIR = os.path.join(LIG_DIR, "sdf")
PDBQT_DIR = os.path.join(LIG_DIR, "pdbqt_pH68")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
POSES_DIR = os.path.join(RESULTS_DIR, "poses")
VINA_EXE = os.path.join(BASE_DIR, "scripts", "vina.exe")

os.makedirs(PDBQT_DIR, exist_ok=True)
os.makedirs(os.path.join(POSES_DIR, "FFAR4"), exist_ok=True)
os.makedirs(os.path.join(POSES_DIR, "CD36"), exist_ok=True)


def main():
    print("=" * 80)
    print("STEP 2 (PROTOCOL A): ALTLOC-CLEAN + SALIVARY pH 6.8 PROTONATION + VINA (--seed 42)")
    print("=" * 80, flush=True)

    # 1. Clean 8G59 (Chain R only; strip Gq chains A,B,Y,S + altLoc B PHE115) at salivary pH 6.8
    ffar4_clean = os.path.join(STRUCT_DIR, "8G59_FFAR4_clean.pdb")
    ffar4_pdbqt = os.path.join(STRUCT_DIR, "8G59_FFAR4_clean.pdbqt")
    tug_crystal = os.path.join(STRUCT_DIR, "TUG891_crystal.pdb")
    n_ffar4, c_ffar4 = clean_pdb_strict(os.path.join(STRUCT_DIR, "8G59.pdb"), ffar4_clean, "R", "YN9", "401", tug_crystal)
    prepare_receptor_pdbqt(ffar4_clean, ffar4_pdbqt, ph=6.8)
    ffar4_box = {"center_x": float(c_ffar4[0]), "center_y": float(c_ffar4[1]), "center_z": float(c_ffar4[2]), "size_x": 22.0, "size_y": 22.0, "size_z": 22.0}
    print(f"[8G59 FFAR4] altLoc-clean Chain R atoms: {n_ffar4} | TUG-891 center: ({c_ffar4[0]:.2f}, {c_ffar4[1]:.2f}, {c_ffar4[2]:.2f})", flush=True)

    # 2. Clean 5LGD (Chain A only; strip PfEMP1 Chain B + glycans + altLoc B) at salivary pH 6.8
    cd36_clean = os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdb")
    cd36_pdbqt = os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdbqt")
    plm_crystal = os.path.join(STRUCT_DIR, "PLM512_crystal.pdb")
    n_cd36, c_cd36 = clean_pdb_strict(os.path.join(STRUCT_DIR, "5LGD.pdb"), cd36_clean, "A", "PLM", "512", plm_crystal)
    prepare_receptor_pdbqt(cd36_clean, cd36_pdbqt, ph=6.8)
    cd36_box = {"center_x": float(c_cd36[0]), "center_y": float(c_cd36[1]), "center_z": float(c_cd36[2]), "size_x": 22.0, "size_y": 22.0, "size_z": 24.0}
    print(f"[5LGD CD36]  altLoc-clean Chain A atoms: {n_cd36} | PLM-512 center: ({c_cd36[0]:.2f}, {c_cd36[1]:.2f}, {c_cd36[2]:.2f})", flush=True)

    admet_df = pd.read_csv(os.path.join(RESULTS_DIR, "admet_and_physicochemical_properties.csv"))
    docking_rows = []
    conv_check_targets = {"TUG-891", "Naringenin", "Resveratrol", "Phloretin", "Spilanthol"}
    conv_rows = []

    for idx, row in admet_df.iterrows():
        name = row["Compound"]
        n_heavy = int(row["Heavy_Atoms"])
        clogp = float(row["LogP_WildmanCrippen"])
        esol_logs = float(row["ESOL_LogS_mol_L"])
        sdf_file = os.path.join(SDF_DIR, row["SDF_File"])
        base_name = row["SDF_File"].replace(".sdf", "")
        lig_pdbqt = os.path.join(PDBQT_DIR, f"{base_name}_pH68.pdbqt")

        prepare_ligand_pdbqt(sdf_file, lig_pdbqt, ph=6.8)

        out_ffar4 = os.path.join(POSES_DIR, "FFAR4", f"{base_name}_FFAR4_out.pdbqt")
        out_cd36 = os.path.join(POSES_DIR, "CD36", f"{base_name}_CD36_out.pdbqt")

        print(f"[{idx+1:02d}/22] Docking {name} (N_HA={n_heavy}, RotB={row['Rotatable_Bonds']}) at pH 6.8 [exh=8, seed=42]...", flush=True)
        m_ffar4 = shift_and_dock_vina(VINA_EXE, ffar4_pdbqt, lig_pdbqt, out_ffar4, ffar4_box, exhaustiveness=8, seed=42, cpu=4)
        m_cd36 = shift_and_dock_vina(VINA_EXE, cd36_pdbqt, lig_pdbqt, out_cd36, cd36_box, exhaustiveness=8, seed=42, cpu=4)

        aff_ffar4 = m_ffar4[0]["affinity"] if m_ffar4 else np.nan
        aff_cd36 = m_cd36[0]["affinity"] if m_cd36 else np.nan
        le_ffar4 = round(-aff_ffar4 / n_heavy, 3)
        le_cd36 = round(-aff_cd36 / n_heavy, 3)
        pkd_ffar4 = round(-aff_ffar4 / 1.364, 2)
        pkd_cd36 = round(-aff_cd36 / 1.364, 2)
        lle_ffar4 = round(pkd_ffar4 - clogp, 2)
        lle_cd36 = round(pkd_cd36 - clogp, 2)
        dual_dg = round(0.60 * aff_ffar4 + 0.40 * aff_cd36, 2)

        le_pass = "Yes" if (le_ffar4 >= 0.30 and le_cd36 >= 0.25) else "No"
        pains_clean = (row["PAINS_Structural_Alert"] == "0 Alerts")

        if "Positive Control" in row["Category"]:
            tier = "Reference Control"
        elif n_heavy >= 30 and (le_ffar4 < 0.30 or esol_logs < -6.0):
            tier = "Disqualified False-Positive (Steric Clash / High-MW Brick-Dust; LE < 0.30)"
        elif le_ffar4 >= 0.30 and esol_logs >= -5.0 and aff_ffar4 <= -7.0 and pains_clean:
            tier = "Tier 1A: PAINS-Free Efficiency-Validated Lingual Lead (LE >= 0.30)"
        elif le_ffar4 >= 0.30 and esol_logs >= -5.0 and aff_ffar4 <= -7.0:
            tier = "Tier 1B: Catechol-PAINS Polyphenolic Hit (Requires Redox Control)"
        else:
            tier = "Tier 2: Low-MW Fragment / Flexible Analog"

        print(f"    -> FFAR4: {aff_ffar4:.2f} (LE={le_ffar4:.3f}, LLE={lle_ffar4:.2f}) | CD36: {aff_cd36:.2f} (LE={le_cd36:.3f}) | {tier[:42]}", flush=True)

        if name in conv_check_targets:
            out_exh12 = os.path.join(POSES_DIR, "FFAR4", f"{base_name}_FFAR4_exh12.pdbqt")
            m12 = shift_and_dock_vina(VINA_EXE, ffar4_pdbqt, lig_pdbqt, out_exh12, ffar4_box, exhaustiveness=12, seed=42, cpu=4)
            aff12 = m12[0]["affinity"] if m12 else np.nan
            pose_rmsd = hungarian_heavy_atom_rmsd(out_ffar4, out_exh12, is_pdbqt_a=True, is_pdbqt_b=True)
            conv_rows.append({"Compound": name, "Exh8_dG": aff_ffar4, "Exh12_dG": aff12, "Delta_dG": round(abs(aff12 - aff_ffar4), 3), "Pose_RMSD_A": pose_rmsd})
            print(f"    -> [Convergence Check exh=12] dG={aff12:.2f} | Pose RMSD vs exh=8: {pose_rmsd:.3f} Å", flush=True)

        docking_rows.append({
            "Compound": name,
            "Category": row["Category"],
            "PubChem_CID": row["PubChem_CID"],
            "MW_Da": row["MW_Da"],
            "Heavy_Atoms": n_heavy,
            "LogP": clogp,
            "TPSA_A2": row["TPSA_A2"],
            "ESOL_LogS": esol_logs,
            "BOILED_Egg_GI": row["BOILED_Egg_GI_Absorption"],
            "PAINS_Alert": row["PAINS_Structural_Alert"],
            "Brenk_Alert": row["Brenk_Structural_Alert"],
            "FFAR4_8G59_Affinity_kcal_mol": aff_ffar4,
            "FFAR4_LE_kcal_mol_HA": le_ffar4,
            "FFAR4_pKd": pkd_ffar4,
            "FFAR4_LLE": lle_ffar4,
            "CD36_5LGD_Affinity_kcal_mol": aff_cd36,
            "CD36_LE_kcal_mol_HA": le_cd36,
            "CD36_pKd": pkd_cd36,
            "CD36_LLE": lle_cd36,
            "Weighted_Dual_Score_kcal_mol": dual_dg,
            "LE_Threshold_Pass_030": le_pass,
            "Lead_Qualification_Tier": tier
        })

    dock_df = pd.DataFrame(docking_rows)
    tier_order = {
        "Tier 1A: PAINS-Free Efficiency-Validated Lingual Lead (LE >= 0.30)": 1,
        "Tier 1B: Catechol-PAINS Polyphenolic Hit (Requires Redox Control)": 2,
        "Tier 2: Low-MW Fragment / Flexible Analog": 3,
        "Reference Control": 4,
        "Disqualified False-Positive (Steric Clash / High-MW Brick-Dust; LE < 0.30)": 5
    }
    dock_df["_s"] = dock_df["Lead_Qualification_Tier"].map(tier_order)
    dock_df = dock_df.sort_values(["_s", "FFAR4_8G59_Affinity_kcal_mol"]).drop(columns=["_s"]).reset_index(drop=True)
    dock_df["Efficiency_Adjusted_Rank"] = dock_df.index + 1

    dock_df.to_csv(os.path.join(RESULTS_DIR, "molecular_docking_results.csv"), index=False)
    pd.DataFrame(conv_rows).to_csv(os.path.join(RESULTS_DIR, "convergence_exh8_vs_exh12.csv"), index=False)
    print(f"\n[SUCCESS] Protocol A docking + convergence verification saved!", flush=True)


if __name__ == "__main__":
    main()
