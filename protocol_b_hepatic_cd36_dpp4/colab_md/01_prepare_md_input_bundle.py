"""
01_prepare_md_input_bundle.py
Prepares bond-order-verified, pH 7.4 protonated 3D SDF complexes directly from Model 1
docked Vina PDBQTs (0.000 Å heavy-atom deviation from docked poses) and packs all inputs
into protocol_b_md_inputs.zip for 1-click 100 ns explicit-water OpenMM MD + MM-GBSA
execution on Google Colab Pro (L4 / A100 GPU).
"""
import os
import sys
import glob
import json
import zipfile
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

COLAB_DIR = os.path.dirname(os.path.abspath(__file__))
PROT_B_DIR = os.path.dirname(COLAB_DIR)
ROOT_DIR = os.path.dirname(PROT_B_DIR)
STRUCT_DIR = os.path.join(ROOT_DIR, "structures")
SDF_DIR = os.path.join(PROT_B_DIR, "ligands", "sdf")
POSES_DPP4 = os.path.join(PROT_B_DIR, "results", "poses", "DPP4_4A5S")
POSES_CD36 = os.path.join(PROT_B_DIR, "results", "poses", "CD36_5LGD")
OBABEL_EXE = r"C:\Users\Abdullah\AppData\Local\Programs\Python\Python312\Scripts\obabel.exe"

LIG_OUT_DIR = os.path.join(COLAB_DIR, "docked_ligands_pH74_sdf")
os.makedirs(LIG_OUT_DIR, exist_ok=True)


def find_file(folder, pattern):
    matches = glob.glob(os.path.join(folder, pattern))
    if not matches:
        raise FileNotFoundError(f"Missing {os.path.join(folder, pattern)}")
    return matches[0]


def extract_model1_heavy_pdb(pdbqt_path, out_pdb):
    lines = []
    in_m1, has_m = False, False
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("MODEL"):
                has_m = True
                if l.strip().split()[-1] == "1":
                    in_m1 = True
                else:
                    break
            elif l.startswith("ENDMDL") and in_m1:
                break
            elif (in_m1 or not has_m) and l.startswith(("ATOM  ", "HETATM")):
                aname = l[12:16].strip()
                adtype = l[77:79].strip() if len(l) >= 78 else ""
                if adtype in ("H", "HD") or aname.startswith("H"):
                    continue
                el_map = {"A": "C", "C": "C", "N": "N", "NA": "N", "O": "O", "OA": "O", "S": "S", "SA": "S", "F": "F", "CL": "CL", "BR": "BR"}
                el = el_map.get(adtype.upper(), "".join(c for c in aname if c.isalpha())[:2].upper() or "C")
                lines.append(f"HETATM{l[6:12]}{l[12:16]} LIG L   1    {l[30:54]}  1.00 20.00          {el:>2s}\n")
    with open(out_pdb, "w", encoding="utf-8") as fo:
        fo.writelines(lines)
        fo.write("END\n")
    return out_pdb


def build_verified_docked_sdf(name, template_src, docked_pdbqt, out_sdf):
    # 1. Protonate template source at pH 7.4 via OpenBabel to obtain exact formal charges & valence
    tmp_tpl_sdf = out_sdf.replace(".sdf", "_tpl.sdf")
    tmp_m1_pdb = out_sdf.replace(".sdf", "_m1.pdb")
    subprocess.run([OBABEL_EXE, template_src, "-O", tmp_tpl_sdf, "-d", "-p", "7.4"], capture_output=True, check=True)
    tpl_h = Chem.SDMolSupplier(tmp_tpl_sdf, removeHs=False)[0]
    assert tpl_h is not None, f"Failed to parse pH 7.4 template for {name}"
    tpl_no_h = Chem.RemoveHs(tpl_h)

    # 2. Load Model 1 docked heavy-atom geometry and transfer exact bond orders & formal charges from pH 7.4 template
    extract_model1_heavy_pdb(docked_pdbqt, tmp_m1_pdb)
    dock_raw = Chem.MolFromPDBFile(tmp_m1_pdb, sanitize=False, removeHs=True)
    assert dock_raw is not None, f"Failed to load Model 1 PDB for {name}"
    dock_bo = AllChem.AssignBondOrdersFromTemplate(tpl_no_h, dock_raw)
    Chem.SanitizeMol(dock_bo)
    Chem.AssignStereochemistryFrom3D(dock_bo)
    dock_h = Chem.AddHs(dock_bo, addCoords=True)
    dock_h.SetProp("_Name", name)

    # 3. Verify 0.000 Å heavy-atom coordinate preservation relative to docked PDBQT
    orig_xyz = np.array([[dock_raw.GetConformer().GetAtomPosition(i).x,
                          dock_raw.GetConformer().GetAtomPosition(i).y,
                          dock_raw.GetConformer().GetAtomPosition(i).z]
                         for i in range(dock_raw.GetNumAtoms())])
    new_xyz = np.array([[dock_h.GetConformer().GetAtomPosition(i).x,
                         dock_h.GetConformer().GetAtomPosition(i).y,
                         dock_h.GetConformer().GetAtomPosition(i).z]
                        for i in range(dock_h.GetNumAtoms()) if dock_h.GetAtomWithIdx(i).GetAtomicNum() > 1])
    rmsd_check = float(np.sqrt(np.mean(np.sum((orig_xyz - new_xyz) ** 2, axis=1))))
    assert rmsd_check < 1e-3, f"Coordinate mismatch on {name}: {rmsd_check}"

    w = Chem.SDWriter(out_sdf)
    w.write(dock_h)
    w.close()

    for tp in (tmp_tpl_sdf, tmp_m1_pdb):
        if os.path.exists(tp):
            os.remove(tp)

    q = Chem.GetFormalCharge(dock_h)
    print(f"  [OK] {name:22s} | heavy={dock_bo.GetNumHeavyAtoms():2d} | total={dock_h.GetNumAtoms():2d} | formal_q={q:+d} | pose_RMSD={rmsd_check:.5f} Å")
    return {"name": name, "sdf": os.path.basename(out_sdf), "heavy_atoms": dock_bo.GetNumHeavyAtoms(), "total_atoms": dock_h.GetNumAtoms(), "formal_charge": int(q)}


def main():
    jobs = [
        ("4A5S_N7F901", "4A5S_DPP4_clean.pdb", os.path.join(STRUCT_DIR, "N7F901_crystal.pdb"), os.path.join(POSES_DPP4, "N7F901_redocked_4A5S.pdbqt")),
        ("4A5S_Teneligliptin", "4A5S_DPP4_clean.pdb", find_file(SDF_DIR, "*_Teneligliptin.sdf"), find_file(POSES_DPP4, "*_Teneligliptin_DPP4.pdbqt")),
        ("4A5S_Naringenin", "4A5S_DPP4_clean.pdb", find_file(SDF_DIR, "*_Naringenin.sdf"), find_file(POSES_DPP4, "*_Naringenin_DPP4.pdbqt")),
        ("4A5S_Berberine", "4A5S_DPP4_clean.pdb", find_file(SDF_DIR, "*_Berberine.sdf"), find_file(POSES_DPP4, "*_Berberine_DPP4.pdbqt")),
        ("5LGD_PLM512", "5LGD_CD36_clean.pdb", find_file(SDF_DIR, "*_Palmitic_acid.sdf"), os.path.join(POSES_CD36, "PLM512_redocked_5LGD.pdbqt")),
        ("5LGD_Linagliptin", "5LGD_CD36_clean.pdb", find_file(SDF_DIR, "*_Linagliptin.sdf"), find_file(POSES_CD36, "*_Linagliptin_CD36.pdbqt")),
        ("5LGD_Teneligliptin", "5LGD_CD36_clean.pdb", find_file(SDF_DIR, "*_Teneligliptin.sdf"), find_file(POSES_CD36, "*_Teneligliptin_CD36.pdbqt")),
        ("5LGD_Naringenin", "5LGD_CD36_clean.pdb", find_file(SDF_DIR, "*_Naringenin.sdf"), find_file(POSES_CD36, "*_Naringenin_CD36.pdbqt")),
    ]

    manifest = []
    for name, rec_pdb, tpl_src, pose_qt in jobs:
        out_sdf = os.path.join(LIG_OUT_DIR, f"{name}_pH74.sdf")
        meta = build_verified_docked_sdf(name, tpl_src, pose_qt, out_sdf)
        meta["receptor_pdb"] = rec_pdb
        manifest.append(meta)

    man_path = os.path.join(COLAB_DIR, "md_manifest.json")
    with open(man_path, "w", encoding="utf-8") as fj:
        json.dump(manifest, fj, indent=2)

    zip_path = os.path.join(COLAB_DIR, "protocol_b_md_inputs.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(man_path, "md_manifest.json")
        zf.write(os.path.join(STRUCT_DIR, "4A5S_DPP4_clean.pdb"), "receptors/4A5S_DPP4_clean.pdb")
        zf.write(os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdb"), "receptors/5LGD_CD36_clean.pdb")
        for m in manifest:
            zf.write(os.path.join(LIG_OUT_DIR, m["sdf"]), f"ligands/{m['sdf']}")

    print(f"\n[SUCCESS] Packed {len(manifest)} verified pH 7.4 MD complexes into {zip_path} ({os.path.getsize(zip_path)/1024:.1f} KB)")


if __name__ == "__main__":
    main()
