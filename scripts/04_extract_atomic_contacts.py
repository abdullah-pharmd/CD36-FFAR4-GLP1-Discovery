"""
04_extract_atomic_contacts.py
Extracts true 3D polar hydrogen-bond contacts (<= 3.5 A), all heavy-atom
pocket contacts (<= 4.0 A), steric clashes (< 2.2 A), and exact symmetry-corrected
Hungarian heavy-atom RMSDs for Protocol A (8G59 FFAR4 & 5LGD CD36).
"""
import os
import sys
import json
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structural_prep_utils import hungarian_heavy_atom_rmsd, analyze_pocket_interactions_3d

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FFAR4_PDB = os.path.join(BASE_DIR, "structures", "8G59_FFAR4_clean.pdb")
CD36_PDB = os.path.join(BASE_DIR, "structures", "5LGD_CD36_clean.pdb")
TUG_CRYSTAL = os.path.join(BASE_DIR, "structures", "TUG891_crystal.pdb")
PLM_CRYSTAL = os.path.join(BASE_DIR, "structures", "PLM512_crystal.pdb")
OUT_CSV = os.path.join(BASE_DIR, "results", "docked_pose_residue_contacts.csv")


def main():
    admet_df = pd.read_csv(os.path.join(BASE_DIR, "results", "admet_and_physicochemical_properties.csv"))
    dock_df = pd.read_csv(os.path.join(BASE_DIR, "results", "molecular_docking_results.csv"))
    sdf_map = dict(zip(admet_df["Compound"], admet_df["SDF_File"]))

    rows = []
    redock_metrics = {}
    for _, row in dock_df.iterrows():
        cname = row["Compound"]
        base_name = sdf_map[cname].replace(".sdf", "")

        ffar4_pose = os.path.join(BASE_DIR, "results", "poses", "FFAR4", f"{base_name}_FFAR4_out.pdbqt")
        cd36_pose = os.path.join(BASE_DIR, "results", "poses", "CD36", f"{base_name}_CD36_out.pdbqt")

        f_int = analyze_pocket_interactions_3d(FFAR4_PDB, ffar4_pose)
        c_int = analyze_pocket_interactions_3d(CD36_PDB, cd36_pose)

        rows.append({
            "Compound": cname,
            "Category": row["Category"],
            "Lead_Qualification_Tier": row["Lead_Qualification_Tier"],
            "FFAR4_Affinity_kcal_mol": row["FFAR4_8G59_Affinity_kcal_mol"],
            "FFAR4_LE": row["FFAR4_LE_kcal_mol_HA"],
            "FFAR4_Polar_HBonds_35A": f_int["polar_hbonds"],
            "FFAR4_Closest_Contacts_le_4A": f_int["all_contacts_4A"],
            "FFAR4_Steric_Clashes": f_int["steric_clashes"],
            "CD36_Affinity_kcal_mol": row["CD36_5LGD_Affinity_kcal_mol"],
            "CD36_LE": row["CD36_LE_kcal_mol_HA"],
            "CD36_Polar_HBonds_35A": c_int["polar_hbonds"],
            "CD36_Closest_Contacts_le_4A": c_int["all_contacts_4A"],
            "CD36_Steric_Clashes": c_int["steric_clashes"]
        })

        if cname == "TUG-891":
            rmsd_tug = hungarian_heavy_atom_rmsd(TUG_CRYSTAL, ffar4_pose, is_pdbqt_a=False, is_pdbqt_b=True)
            redock_metrics["TUG891_8G59_Hungarian_RMSD_A"] = rmsd_tug
            print(f"[Redocking Validation] 8G59 TUG-891 symmetry-corrected Hungarian RMSD = {rmsd_tug:.3f} Å")
        if cname == "Palmitic acid":
            rmsd_plm = hungarian_heavy_atom_rmsd(PLM_CRYSTAL, cd36_pose, is_pdbqt_a=False, is_pdbqt_b=True)
            redock_metrics["PLM512_5LGD_Hungarian_RMSD_A"] = rmsd_plm
            print(f"[Redocking Validation] 5LGD PLM-512 symmetry-corrected Hungarian RMSD = {rmsd_plm:.3f} Å")

    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    with open(os.path.join(BASE_DIR, "results", "redocking_metrics.json"), "w") as fj:
        json.dump(redock_metrics, fj, indent=2)
    print(f"[SUCCESS] Wrote {len(rows)} 3D interaction & clash profiles to: {OUT_CSV}")


if __name__ == "__main__":
    main()
