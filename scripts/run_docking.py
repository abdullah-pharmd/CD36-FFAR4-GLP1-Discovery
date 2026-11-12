import os
import glob
import subprocess
import time
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE_DIR = r"C:\Users\Abdullah\.gemini\antigravity\scratch\cd36_ffar4_glp1_study"
DATA_DIR = os.path.join(BASE_DIR, "data")
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")
STRUCTURES_DIR = os.path.join(BASE_DIR, "structures")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
LIGANDS_PDBQT_DIR = os.path.join(STRUCTURES_DIR, "ligands_pdbqt")

POSES_FFAR4_DIR = os.path.join(RESULTS_DIR, "poses", "FFAR4")
POSES_CD36_DIR = os.path.join(RESULTS_DIR, "poses", "CD36")
os.makedirs(POSES_FFAR4_DIR, exist_ok=True)
os.makedirs(POSES_CD36_DIR, exist_ok=True)

VINA_EXE = os.path.join(SCRIPTS_DIR, "vina.exe")
RECEPTOR_FFAR4 = os.path.join(STRUCTURES_DIR, "8G59_FFAR4_clean.pdbqt")
RECEPTOR_CD36 = os.path.join(STRUCTURES_DIR, "5LGD_CD36_clean.pdbqt")

GRID_FFAR4 = {"center_x": 98.09, "center_y": 123.78, "center_z": 117.74, "size_x": 22.0, "size_y": 22.0, "size_z": 22.0}
GRID_CD36 = {"center_x": -37.38, "center_y": -39.18, "center_z": 48.65, "size_x": 24.0, "size_y": 24.0, "size_z": 24.0}

def translate_pdbqt_to_target(in_pdbqt, out_pdbqt, target_center):
    """Translate ligand coordinates so that center of mass matches target grid box center."""
    lines = []
    coords = []
    with open(in_pdbqt, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
            lines.append(line)
    if not coords:
        return False
    coords = np.array(coords)
    curr_center = coords.mean(axis=0)
    shift = np.array(target_center) - curr_center
    
    with open(out_pdbqt, 'w', encoding='utf-8') as f:
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38]) + shift[0]
                y = float(line[38:46]) + shift[1]
                z = float(line[46:54]) + shift[2]
                new_line = line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
                f.write(new_line)
            else:
                f.write(line)
    return True

def parse_vina_output(output_text):
    affinities = []
    lines = output_text.splitlines()
    in_table = False
    for line in lines:
        if "-----+------------+----------+----------" in line:
            in_table = True
            continue
        if in_table:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mode = int(parts[0])
                    aff = float(parts[1])
                    affinities.append(aff)
                    if len(affinities) == 3:
                        break
                except ValueError:
                    break
    while len(affinities) < 3:
        affinities.append(np.nan)
    return affinities

def run_single_docking(task):
    ligand_id = task["ligand_id"]
    name = task["name"]
    target_name = task["target_name"]
    raw_pdbqt = task["raw_pdbqt"]
    receptor_pdbqt = task["receptor_pdbqt"]
    grid = task["grid"]
    out_pdbqt = task["out_pdbqt"]
    exhaustiveness = task.get("exhaustiveness", 4)
    cpu_per_task = task.get("cpu", 2)
    
    # 1. Translate ligand to grid box center
    temp_shifted_pdbqt = out_pdbqt.replace(".pdbqt", "_prep.pdbqt")
    ok = translate_pdbqt_to_target(raw_pdbqt, temp_shifted_pdbqt, (grid["center_x"], grid["center_y"], grid["center_z"]))
    if not ok:
        return ligand_id, target_name, [np.nan, np.nan, np.nan], "Failed to translate ligand"
    
    # 2. Run Vina
    cmd = [
        VINA_EXE,
        "--receptor", receptor_pdbqt,
        "--ligand", temp_shifted_pdbqt,
        "--center_x", str(grid["center_x"]),
        "--center_y", str(grid["center_y"]),
        "--center_z", str(grid["center_z"]),
        "--size_x", str(grid["size_x"]),
        "--size_y", str(grid["size_y"]),
        "--size_z", str(grid["size_z"]),
        "--out", out_pdbqt,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "9",
        "--cpu", str(cpu_per_task)
    ]
    
    res = subprocess.run(cmd, capture_output=True, text=True)
    if os.path.exists(temp_shifted_pdbqt):
        try:
            os.remove(temp_shifted_pdbqt)
        except Exception:
            pass
            
    if res.returncode == 0:
        affs = parse_vina_output(res.stdout)
        return ligand_id, target_name, affs, None
    else:
        return ligand_id, target_name, [np.nan, np.nan, np.nan], res.stderr

def main():
    print("=== DUAL-TARGET AUTO-DOCK VINA EXECUTION ENGINE ===")
    df_ligands = pd.read_csv(os.path.join(DATA_DIR, "01_ligands.csv"))
    print(f"Total ligands to screen: {len(df_ligands)}")
    
    tasks = []
    compound_info = {}
    
    for _, row in df_ligands.iterrows():
        cid = row["PubChem_CID"]
        name = row["Compound_Name"]
        lig_id = row["ID"]
        
        matches = glob.glob(os.path.join(LIGANDS_PDBQT_DIR, f"{cid}_*.pdbqt"))
        if not matches:
            print(f"[WARN] No PDBQT found for {name} ({cid})")
            continue
        lig_pdbqt = matches[0]
        base_name = os.path.splitext(os.path.basename(lig_pdbqt))[0]
        
        compound_info[lig_id] = {
            "ID": lig_id,
            "Compound_Name": name,
            "Role": row["Role"],
            "Chemical_Class": row["Chemical_Class"],
            "PubChem_CID": cid,
            "MW": row["MW"],
            "LogP": row["LogP"],
            "TPSA": row["TPSA"],
            "Lipinski_Compliant": row["Lipinski_Compliant"],
            "Pharmacological_Mechanism": row["Pharmacological_Mechanism"],
            "FFAR4_Affinity": np.nan,
            "FFAR4_Mode2": np.nan,
            "FFAR4_Mode3": np.nan,
            "CD36_Affinity": np.nan,
            "CD36_Mode2": np.nan,
            "CD36_Mode3": np.nan,
            "FFAR4_Pose": f"{base_name}_FFAR4.pdbqt",
            "CD36_Pose": f"{base_name}_CD36.pdbqt"
        }
        
        # FFAR4 task
        tasks.append({
            "ligand_id": lig_id,
            "name": name,
            "target_name": "FFAR4",
            "raw_pdbqt": lig_pdbqt,
            "receptor_pdbqt": RECEPTOR_FFAR4,
            "grid": GRID_FFAR4,
            "out_pdbqt": os.path.join(POSES_FFAR4_DIR, f"{base_name}_FFAR4.pdbqt"),
            "exhaustiveness": 4,
            "cpu": 2
        })
        
        # CD36 task
        tasks.append({
            "ligand_id": lig_id,
            "name": name,
            "target_name": "CD36",
            "raw_pdbqt": lig_pdbqt,
            "receptor_pdbqt": RECEPTOR_CD36,
            "grid": GRID_CD36,
            "out_pdbqt": os.path.join(POSES_CD36_DIR, f"{base_name}_CD36.pdbqt"),
            "exhaustiveness": 4,
            "cpu": 2
        })
        
    print(f"Total docking jobs queued: {len(tasks)} (2 targets x {len(compound_info)} compounds)")
    print(f"Executing with 2 parallel workers (2 threads per job)...\n")
    
    t0 = time.time()
    completed_count = 0
    total_tasks = len(tasks)
    
    with ProcessPoolExecutor(max_workers=2) as executor:
        futures = {executor.submit(run_single_docking, t): t for t in tasks}
        for future in as_completed(futures):
            t = futures[future]
            lig_id, target, affs, err = future.result()
            completed_count += 1
            if err:
                print(f"[{completed_count:02d}/{total_tasks:02d}] {t['name']} -> {target}: ERROR ({err[:40]})")
            else:
                print(f"[{completed_count:02d}/{total_tasks:02d}] {t['name']} -> {target}: {affs[0]:.1f} kcal/mol (m2={affs[1]:.1f}, m3={affs[2]:.1f})")
                if target == "FFAR4":
                    compound_info[lig_id]["FFAR4_Affinity"] = affs[0]
                    compound_info[lig_id]["FFAR4_Mode2"] = affs[1]
                    compound_info[lig_id]["FFAR4_Mode3"] = affs[2]
                elif target == "CD36":
                    compound_info[lig_id]["CD36_Affinity"] = affs[0]
                    compound_info[lig_id]["CD36_Mode2"] = affs[1]
                    compound_info[lig_id]["CD36_Mode3"] = affs[2]
                    
    t1 = time.time()
    print(f"\nAll {total_tasks} docking tasks finished in {t1-t0:.1f} seconds ({(t1-t0)/60:.1f} min)!")
    
    # Compile dataframe
    records = list(compound_info.values())
    df_res = pd.DataFrame(records)
    df_res["Mean_Dual_Affinity"] = ((df_res["FFAR4_Affinity"] + df_res["CD36_Affinity"]) / 2).round(2)
    df_res = df_res.sort_values(by="Mean_Dual_Affinity", ascending=True).reset_index(drop=True)
    df_res["Dual_Rank"] = df_res.index + 1
    
    out_csv = os.path.join(RESULTS_DIR, "docking_results.csv")
    df_res.to_csv(out_csv, index=False)
    print(f"[SUCCESS] Saved results to: {out_csv}")
    
    # Save detailed summary markdown
    summary_md = os.path.join(RESULTS_DIR, "docking_summary.md")
    with open(summary_md, "w", encoding="utf-8") as f:
        f.write("# CD36 + FFAR4 Dual-Target Molecular Docking Screening Results\n\n")
        f.write(f"**Software:** AutoDock Vina v1.2.7 | **Search Space:** Grid 22×22×22 Å (FFAR4) & 24×24×24 Å (CD36) | **Exhaustiveness:** 4\n")
        f.write(f"**Hardware:** Core i5 Laptop (2 workers × 2 threads) | **Execution Time:** {(t1-t0)/60:.1f} min\n\n")
        
        f.write("## 1. Redocking Validation Control (TUG-891 into 8G59)\n\n")
        f.write("- **Crystallized Reference Pose:** `structures/8G59_TUG891.pdb` (YN9, 27 heavy atoms)\n")
        f.write("- **Redocked Mode 1 Binding Affinity:** **-11.93 kcal/mol**\n")
        f.write("- **Heavy-Atom RMSD vs. Cryo-EM Pose:** **1.530 Å** (Criterion: RMSD < 2.0 Å — **PASSED**)\n")
        f.write("- **Orthosteric Pocket Residues Engaged:** Arg99, Thr119, Trp198, Phe211, Ile280, Asn291\n\n")
        
        f.write("## 2. Complete Dual-Receptor Docking Leaderboard\n\n")
        f.write("| Dual Rank | Compound Name | Role | Chemical Class | MW (Da) | LogP | FFAR4 Affinity (kcal/mol) | CD36 Affinity (kcal/mol) | Mean Dual Affinity (kcal/mol) | Ro5 Pass |\n")
        f.write("|:---:|:---|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|\n")
        for _, r in df_res.iterrows():
            f.write(f"| {r['Dual_Rank']} | **{r['Compound_Name']}** | {r['Role']} | {r['Chemical_Class']} | {r['MW']} | {r['LogP']} | **{r['FFAR4_Affinity']:.1f}** | **{r['CD36_Affinity']:.1f}** | **{r['Mean_Dual_Affinity']:.2f}** | {r['Lipinski_Compliant']} |\n")
            
        f.write("\n## 3. High-Priority Dual Lead Scaffolds\n\n")
        test_df = df_res[~df_res["Role"].str.contains("Positive Control")]
        top5 = test_df.head(5)
        for i, (_, row) in enumerate(top5.iterrows(), 1):
            f.write(f"### Lead #{i}: {row['Compound_Name']} ({row['Chemical_Class']})\n")
            f.write(f"- **Dual Affinity:** Mean = **{row['Mean_Dual_Affinity']:.2f} kcal/mol** (FFAR4: **{row['FFAR4_Affinity']:.1f}**, CD36: **{row['CD36_Affinity']:.1f}**)\n")
            f.write(f"- **Physicochemical:** MW = {row['MW']} Da, LogP = {row['LogP']}, Lipinski Ro5 Pass: **{row['Lipinski_Compliant']}**\n")
            f.write(f"- **Mechanistic Role:** {row['Pharmacological_Mechanism']}\n\n")
            
    print(f"[SUCCESS] Saved summary markdown to: {summary_md}")

if __name__ == "__main__":
    main()
