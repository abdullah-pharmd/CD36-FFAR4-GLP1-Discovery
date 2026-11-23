"""
03_hardening_audits_protocol_b.py
Comprehensive referee-hardening audit suite for Protocol B (CD36 5LGD + DPP-4 4A5S, pH 7.4):
  1. Updates Protocol B master ADMET CSV with exact published Daina & Zoete (2016) rotated
     BOILED-Egg ellipses (preserving all docking scores & tiers) + verifies solubility units.
  2. Executes N=60 DPP-4 PubChem BioAssay (GeneID: 1803) deep statistical & scaffold audit:
       - Stratified 20,000-resample percentile bootstrap 95% CIs (seed=42) for Vina dG,
         full-precision LE, and 5 descriptor baselines (MW, Heavy_Atoms, WLOGP, TPSA, RotB)
       - Paired bootstrap 95% CIs for AUC(Vina) - AUC(LE) and AUC(Vina) - AUC(baselines)
       - PR-AUC, full-precision Top-6 (10%) & Top-12 (20%) recovery, tie-free diagnostics
       - Full 2x2 confusion matrix at Tier 1A DPP-4 gate (dG <= -8.0 & LE >= 0.280)
       - Bemis-Murcko scaffold diversity, cross-class overlap, & 1-per-scaffold deduplicated AUC
       - Full record-level provenance CSV (CID, SMILES, Full_InChIKey, Murcko_Scaffold, units)
  3. Executes 4A5S Val546 alternate-conformation (altLoc A vs B, 0.50/0.50 occupancy)
     fixed-pose + redocking check AND 3-seed reproducibility audit (seeds 42, 137, 2026)
     across N7F-901 and flagship Tier 1 DPP-4 leads.
"""
import os
import sys
import json
import math
import glob
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind
from sklearn.metrics import average_precision_score
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

PROT_B_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(PROT_B_DIR)
sys.path.insert(0, os.path.join(ROOT_DIR, "scripts"))

from structural_prep_utils import (
    BOILED_EGG_HIA, BOILED_EGG_BBB, point_in_boiled_egg_ellipse,
    clean_pdb_strict, prepare_receptor_pdbqt, shift_and_dock_vina,
    hungarian_heavy_atom_rmsd, compute_admet_descriptors
)

RESULTS_DIR = os.path.join(PROT_B_DIR, "results")
NET_DIR = os.path.join(RESULTS_DIR, "network")
STRUCT_DIR = os.path.join(ROOT_DIR, "structures")
POSES_DPP4 = os.path.join(RESULTS_DIR, "poses", "DPP4_4A5S")
LIG_PH74_DIR = os.path.join(PROT_B_DIR, "ligands", "pdbqt_pH74")
VINA_EXE = os.path.join(ROOT_DIR, "scripts", "vina.exe")


def exact_auc_vec(pos, neg):
    """Vectorized Mann-Whitney AUC with exact 0.5 tie credit."""
    diff = pos[:, :, None] - neg[:, None, :]
    return ((diff > 0).sum(axis=(1, 2)) + 0.5 * (diff == 0).sum(axis=(1, 2))) / (pos.shape[1] * neg.shape[1])


def exact_auc_1d(y, s):
    p = s[y == 1][None, :]
    n = s[y == 0][None, :]
    return float(exact_auc_vec(p, n)[0])


def parse_pdb_res_heavy(pdb_path, chain, resnum, alt):
    pts = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("ATOM  ") and l[21] == chain and l[22:26].strip() == str(resnum):
                if l[16] in (" ", alt) and l[76:78].strip() != "H":
                    pts.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(pts, dtype=float)


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


def step1_update_protocol_b_boiled_egg():
    csv_path = os.path.join(RESULTS_DIR, "protocol_b_cd36_dpp4_complete_results.csv")
    df = pd.read_csv(csv_path)
    changed = []
    for idx, row in df.iterrows():
        smi = str(row["SMILES"])
        ad = compute_admet_descriptors(smi)
        # Verify solubility unit invariant
        assert abs(float(row["Estimated_Solubility_ug_mL"]) - float(row["Estimated_Solubility_mg_mL"]) * 1000.0) <= 0.15
        old_hia = str(row["BOILED_Egg_GI_Absorption"])
        old_bbb = str(row["BOILED_Egg_BBB"])
        if "Quat" in old_hia:
            continue
        if ad["BOILED_Egg_GI_Absorption"] != old_hia or ad["BOILED_Egg_BBB"] != old_bbb:
            changed.append({
                "Compound": row["Compound"],
                "TPSA": float(row["TPSA_A2"]),
                "WLOGP": float(row["LogP_WildmanCrippen"]),
                "Old_HIA": old_hia, "New_HIA": ad["BOILED_Egg_GI_Absorption"],
                "Old_BBB": old_bbb, "New_BBB": ad["BOILED_Egg_BBB"],
            })
            # Assert ZERO HIA changes so Tier 1A gates are 100% invariant
            assert ad["BOILED_Egg_GI_Absorption"] == old_hia, f"Unexpected HIA change on {row['Compound']}"
            df.at[idx, "BOILED_Egg_GI_Absorption"] = ad["BOILED_Egg_GI_Absorption"]
            df.at[idx, "BOILED_Egg_BBB"] = ad["BOILED_Egg_BBB"]
    df.to_csv(csv_path, index=False)
    print(f"[STEP 1] Protocol B BOILED-Egg updated to exact Daina & Zoete (2016) ellipses. BBB updates: {len(changed)} (0 HIA/Tier changes).")
    for c in changed:
        print(f"   - {c['Compound']}: BBB {c['Old_BBB']} -> {c['New_BBB']} (TPSA={c['TPSA']:.2f}, WLOGP={c['WLOGP']:.2f})")
    return changed


def step2_audit_dpp4_roc_benchmark():
    pool_path = os.path.join(NET_DIR, "dpp4_chembl_benchmark_pool.json")
    roc_csv = os.path.join(RESULTS_DIR, "roc_validation_dpp4_chembl284.csv")
    with open(pool_path, "r", encoding="utf-8") as f:
        pool = json.load(f)
    df_scores = pd.read_csv(roc_csv)

    meta_map = {}
    for lbl, grp in [(1, pool["actives"]), (0, pool["inactives"])]:
        for it in grp:
            meta_map[it["chembl_id"]] = (lbl, it)

    prov_rows = []
    for _, r in df_scores.iterrows():
        cid_str = r["ChEMBL_ID"]
        lbl, it = meta_map[cid_str]
        smi = it["smiles"]
        mol = Chem.MolFromSmiles(smi)
        ik_full = Chem.MolToInchiKey(mol)
        scaf_mol = MurckoScaffold.GetScaffoldForMol(mol)
        scaf_smi = Chem.MolToSmiles(scaf_mol, isomericSmiles=False) if scaf_mol is not None else "ACYCLIC"
        wlogp = round(Crippen.MolLogP(mol), 3)
        tpsa = round(rdMolDescriptors.CalcTPSA(mol), 2)
        rotb = int(Lipinski.NumRotatableBonds(mol))
        nha = int(mol.GetNumHeavyAtoms())
        dg = float(r["Vina_Affinity_kcal_mol"])
        le_full = -dg / nha
        prov_rows.append({
            "Target": "DPP4_4A5S_PUBCHEM1803",
            "PubChem_CID": int(it["cid"]),
            "Compound_ID": cid_str,
            "True_Active": int(r["True_Active"]),
            "Activity_Representation": "pChEMBL (-log10 M)" if int(r["True_Active"]) == 1 else "IC50/Ki_or_HTS_Threshold_nM",
            "Activity_Operator": ">=" if int(r["True_Active"]) == 1 else ">=",
            "Experimental_Val": float(r["Experimental_Val"]),
            "Source_Database": "PubChem BioAssay (Human DPP4 GeneID: 1803; covalent warheads SMARTS-filtered)",
            "Isomeric_SMILES": Chem.MolToSmiles(mol, isomericSmiles=True),
            "InChIKey_Full": ik_full,
            "Murcko_Scaffold": scaf_smi,
            "MW_Da": round(float(Descriptors.MolWt(mol)), 2),
            "Heavy_Atoms": nha,
            "WLOGP": wlogp,
            "TPSA_A2": tpsa,
            "RotB": rotb,
            "Vina_Affinity_kcal_mol": dg,
            "Ligand_Efficiency_Full": round(le_full, 6),
            "Ligand_Efficiency": round(le_full, 3),
        })

    df_prov = pd.DataFrame(prov_rows)
    out_prov_csv = os.path.join(RESULTS_DIR, "roc_validation_dpp4_pubchem1803.csv")
    df_prov.to_csv(out_prov_csv, index=False)

    y = df_prov["True_Active"].values
    scores = {
        "Vina_dG": -df_prov["Vina_Affinity_kcal_mol"].values,
        "LE_Full": df_prov["Ligand_Efficiency_Full"].values,
        "Baseline_MW": df_prov["MW_Da"].values,
        "Baseline_HeavyAtoms": df_prov["Heavy_Atoms"].values.astype(float),
        "Baseline_WLOGP": df_prov["WLOGP"].values,
        "Baseline_TPSA": df_prov["TPSA_A2"].values,
        "Baseline_RotB": df_prov["RotB"].values.astype(float),
    }

    rng = np.random.default_rng(42)
    B = 20000
    idx1 = np.where(y == 1)[0]
    idx0 = np.where(y == 0)[0]
    b1 = rng.integers(0, len(idx1), size=(B, len(idx1)))
    b0 = rng.integers(0, len(idx0), size=(B, len(idx0)))

    boot_aucs = {}
    metrics_by_score = {}
    for k, s in scores.items():
        s1, s0 = s[idx1], s[idx0]
        obs = exact_auc_1d(y, s)
        arr = exact_auc_vec(s1[b1], s0[b0])
        boot_aucs[k] = arr
        lo, hi = np.percentile(arr, [2.5, 97.5])
        pr_auc = float(average_precision_score(y, s))
        _, p_mw_u = mannwhitneyu(s1, s0, alternative="greater")
        _, p_two = mannwhitneyu(s1, s0, alternative="two-sided")
        # Full-precision tie-aware top-6 and top-12
        order = np.argsort(-s, kind="mergesort")
        top6_hits = int(y[order[:6]].sum())
        top12_hits = int(y[order[:12]].sum())
        tie_at_6 = bool(len(s) > 6 and s[order[5]] == s[order[6]])
        tie_at_12 = bool(len(s) > 12 and s[order[11]] == s[order[12]])
        metrics_by_score[k] = {
            "ROC_AUC": round(obs, 4),
            "Bootstrap_95CI": [round(float(lo), 4), round(float(hi), 4)],
            "PR_AUC": round(pr_auc, 4),
            "MannWhitney_p_one_sided": float(f"{p_mw_u:.4e}"),
            "MannWhitney_p_two_sided": float(f"{p_two:.4e}"),
            "Top6_Actives": top6_hits,
            "Top6_Boundary_Tie": tie_at_6,
            "Top12_Actives": top12_hits,
            "Top12_Boundary_Tie": tie_at_12,
            "EF10pct": round((top6_hits / 6.0) / y.mean(), 3),
            "EF20pct": round((top12_hits / 12.0) / y.mean(), 3),
        }

    # Paired bootstrap differences (Vina vs LE & vs baselines)
    paired_diffs = {}
    for comp in ["LE_Full", "Baseline_MW", "Baseline_HeavyAtoms", "Baseline_WLOGP", "Baseline_TPSA", "Baseline_RotB"]:
        diff_arr = boot_aucs["Vina_dG"] - boot_aucs[comp]
        obs_d = metrics_by_score["Vina_dG"]["ROC_AUC"] - metrics_by_score[comp]["ROC_AUC"]
        lo, hi = np.percentile(diff_arr, [2.5, 97.5])
        paired_diffs[f"Vina_minus_{comp}"] = {
            "Delta_AUC": round(float(obs_d), 4),
            "Paired_Bootstrap_95CI": [round(float(lo), 4), round(float(hi), 4)],
            "Prob_Vina_Greater": round(float(np.mean(diff_arr > 0)), 4),
        }

    # Tier 1A DPP-4 Gate Confusion Matrix (dG <= -8.0 & LE >= 0.280)
    gate_masks = {
        "Tier1A_dG_le_minus8_and_LE_ge_028": (df_prov["Vina_Affinity_kcal_mol"] <= -8.0) & (df_prov["Ligand_Efficiency_Full"] >= 0.280),
        "dG_le_minus8_only": df_prov["Vina_Affinity_kcal_mol"] <= -8.0,
        "LE_ge_028_only": df_prov["Ligand_Efficiency_Full"] >= 0.280,
        "dG_le_minus8_5_strict": df_prov["Vina_Affinity_kcal_mol"] <= -8.5,
    }
    confusion = {}
    for gname, mask in gate_masks.items():
        m = mask.values
        tp = int((m & (y == 1)).sum())
        fn = int((~m & (y == 1)).sum())
        fp = int((m & (y == 0)).sum())
        tn = int((~m & (y == 0)).sum())
        sens = tp / (tp + fn)
        spec = tn / (tn + fp)
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mcc = ((tp * tn) - (fp * fn)) / denom if denom > 0 else 0.0
        confusion[gname] = {
            "TP": tp, "FN": fn, "FP": fp, "TN": tn,
            "Sensitivity": round(sens, 4),
            "Specificity": round(spec, 4),
            "Precision": round(prec, 4),
            "MCC": round(mcc, 4),
        }

    # Bemis-Murcko Scaffold Diversity & 1-per-Scaffold Collapsed AUC
    act_scafs = set(df_prov[df_prov["True_Active"] == 1]["Murcko_Scaffold"])
    inact_scafs = set(df_prov[df_prov["True_Active"] == 0]["Murcko_Scaffold"])
    overlap_scafs = sorted(act_scafs & inact_scafs)
    scaf_collapsed = (
        df_prov.groupby(["True_Active", "Murcko_Scaffold"], as_index=False)
        .agg(Vina_Median=("Vina_Affinity_kcal_mol", "median"), LE_Median=("Ligand_Efficiency_Full", "median"), Count=("PubChem_CID", "count"))
    )
    ys = scaf_collapsed["True_Active"].values
    sa = -scaf_collapsed["Vina_Median"].values
    sl = scaf_collapsed["LE_Median"].values
    auc_scaf_vina = exact_auc_1d(ys, sa)
    auc_scaf_le = exact_auc_1d(ys, sl)
    _, p_scaf_vina = mannwhitneyu(sa[ys == 1], sa[ys == 0], alternative="greater")
    _, p_scaf_le = mannwhitneyu(sl[ys == 1], sl[ys == 0], alternative="greater")
    # Also bootstrap the scaffold-collapsed AUC
    s_i1 = np.where(ys == 1)[0]
    s_i0 = np.where(ys == 0)[0]
    sb1 = rng.integers(0, len(s_i1), size=(B, len(s_i1)))
    sb0 = rng.integers(0, len(s_i0), size=(B, len(s_i0)))
    scaf_boot_vina = exact_auc_vec(sa[s_i1][sb1], sa[s_i0][sb0])
    scaf_ci = np.percentile(scaf_boot_vina, [2.5, 97.5])

    # Descriptor balance t-tests
    act_df = df_prov[df_prov["True_Active"] == 1]
    inact_df = df_prov[df_prov["True_Active"] == 0]
    balance = {}
    for col in ["MW_Da", "Heavy_Atoms", "WLOGP", "TPSA_A2", "RotB"]:
        _, p_t = ttest_ind(act_df[col], inact_df[col])
        balance[col] = {
            "Active_Mean": round(float(act_df[col].mean()), 2),
            "Inactive_Mean": round(float(inact_df[col].mean()), 2),
            "t_test_p": round(float(p_t), 4),
        }

    audit = {
        "Benchmark": "DPP4_4A5S_PUBCHEM1803",
        "N_Total": int(len(df_prov)),
        "N_Actives": int((y == 1).sum()),
        "N_Inactives": int((y == 0).sum()),
        "Bootstrap_Resamples": B,
        "Bootstrap_Seed": 42,
        "Tie_Policy": "Full double-precision scores; exact 0.5 tie credit in Mann-Whitney AUC",
        "Descriptor_Balance": balance,
        "Metrics_By_Score": metrics_by_score,
        "Paired_Bootstrap_Diffs": paired_diffs,
        "Gate_Confusion_Matrices": confusion,
        "Bemis_Murcko_Scaffold_Audit": {
            "Unique_Active_Scaffolds": len(act_scafs),
            "Unique_Inactive_Scaffolds": len(inact_scafs),
            "Cross_Class_Shared_Scaffolds": len(overlap_scafs),
            "Scaffold_Collapsed_N": int(len(scaf_collapsed)),
            "Scaffold_Collapsed_Vina_AUC": round(auc_scaf_vina, 4),
            "Scaffold_Collapsed_Vina_95CI": [round(float(scaf_ci[0]), 4), round(float(scaf_ci[1]), 4)],
            "Scaffold_Collapsed_Vina_p": float(f"{p_scaf_vina:.4e}"),
            "Scaffold_Collapsed_LE_AUC": round(auc_scaf_le, 4),
            "Scaffold_Collapsed_LE_p": float(f"{p_scaf_le:.4e}"),
        },
    }

    out_json = os.path.join(RESULTS_DIR, "roc_dpp4_pubchem1803_audit.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(audit, f, indent=2)

    print(f"[STEP 2] DPP-4 N=60 Audit saved: {out_json}")
    print(f"   Vina AUC = {metrics_by_score['Vina_dG']['ROC_AUC']} {metrics_by_score['Vina_dG']['Bootstrap_95CI']} | PR-AUC = {metrics_by_score['Vina_dG']['PR_AUC']}")
    print(f"   LE   AUC = {metrics_by_score['LE_Full']['ROC_AUC']} {metrics_by_score['LE_Full']['Bootstrap_95CI']} | Paired Diff = {paired_diffs['Vina_minus_LE_Full']}")
    print(f"   Baselines: MW={metrics_by_score['Baseline_MW']['ROC_AUC']}, WLOGP={metrics_by_score['Baseline_WLOGP']['ROC_AUC']}, TPSA={metrics_by_score['Baseline_TPSA']['ROC_AUC']}, RotB={metrics_by_score['Baseline_RotB']['ROC_AUC']}")
    print(f"   Murcko Scaffolds: {len(act_scafs)}/30 active, {len(inact_scafs)}/30 inactive, overlap={len(overlap_scafs)} | Collapsed AUC={auc_scaf_vina:.4f} (p={p_scaf_vina:.4g}, 95% CI {np.round(scaf_ci, 3)})")
    print(f"   Tier 1A Gate Confusion: {confusion['Tier1A_dG_le_minus8_and_LE_ge_028']}")
    return audit


def step3_4a5s_altloc_and_multiseed():
    pdb_raw = os.path.join(STRUCT_DIR, "4A5S.pdb")
    clean_b = os.path.join(STRUCT_DIR, "4A5S_DPP4_clean_altB.pdb")
    rec_a = os.path.join(STRUCT_DIR, "4A5S_DPP4_pH74.pdbqt")
    rec_b = os.path.join(STRUCT_DIR, "4A5S_DPP4_pH74_altB.pdbqt")

    clean_pdb_strict(pdb_raw, clean_b, chain_id="A", keep_altloc="B")
    ok_b = prepare_receptor_pdbqt(clean_b, rec_b, ph=7.4)
    assert ok_b, "Failed to prepare 4A5S altLoc B receptor PDBQT"

    # Fixed-pose distance to Val546 A vs B across ALL 24 docked DPP-4 poses
    val_a = parse_pdb_res_heavy(pdb_raw, "A", 546, "A")
    val_b = parse_pdb_res_heavy(pdb_raw, "A", 546, "B")
    fixed_rows = []
    for p in sorted(glob.glob(os.path.join(POSES_DPP4, "*_DPP4.pdbqt"))):
        pts = parse_pdbqt_model1_heavy(p)
        if len(pts):
            da = float(np.min(np.linalg.norm(pts[:, None, :] - val_a[None, :, :], axis=-1)))
            db = float(np.min(np.linalg.norm(pts[:, None, :] - val_b[None, :, :], axis=-1)))
            fixed_rows.append({"Pose_File": os.path.basename(p), "Min_Dist_Val546_A_A": round(da, 3), "Min_Dist_Val546_B_A": round(db, 3), "Clash_in_B_lt_2p2": db < 2.2})
    any_clash_b = any(r["Clash_in_B_lt_2p2"] for r in fixed_rows)
    min_all_b = min(r["Min_Dist_Val546_B_A"] for r in fixed_rows) if fixed_rows else np.nan
    print(f"[STEP 3a] Fixed-pose Val546 A vs B check across {len(fixed_rows)} DPP-4 poses: min dist in B = {min_all_b:.3f} A | any < 2.2 A clash = {any_clash_b}")

    # Redock N7F-901 + key Tier 1 leads in 4A5S Conformer B (seed=42) + 3 seeds (42, 137, 2026) in Conformer A
    box_4a5s = {"center_x": 14.14, "center_y": 35.87, "center_z": 54.84, "size_x": 22.0, "size_y": 22.0, "size_z": 22.0}
    tmp_dir = os.path.join(RESULTS_DIR, "poses", "audit_altloc_multiseed")
    os.makedirs(tmp_dir, exist_ok=True)

    targets = [
        ("N7F-901 (Crystal)", os.path.join(STRUCT_DIR, "N7F901_crystal_pH74.pdbqt"), os.path.join(POSES_DPP4, "N7F901_redocked_DPP4.pdbqt")),
        ("Linagliptin", os.path.join(LIG_PH74_DIR, "10096344_Linagliptin_pH74.pdbqt"), os.path.join(POSES_DPP4, "10096344_Linagliptin_DPP4.pdbqt")),
        ("Sitagliptin", os.path.join(LIG_PH74_DIR, "4369359_Sitagliptin_pH74.pdbqt"), os.path.join(POSES_DPP4, "4369359_Sitagliptin_DPP4.pdbqt")),
        ("Teneligliptin", os.path.join(LIG_PH74_DIR, "11949652_Teneligliptin_pH74.pdbqt"), os.path.join(POSES_DPP4, "11949652_Teneligliptin_DPP4.pdbqt")),
        ("Naringenin", os.path.join(LIG_PH74_DIR, "439246_Naringenin_pH74.pdbqt"), os.path.join(POSES_DPP4, "439246_Naringenin_DPP4.pdbqt")),
        ("Berberine", os.path.join(LIG_PH74_DIR, "2353_Berberine_pH74.pdbqt"), os.path.join(POSES_DPP4, "2353_Berberine_DPP4.pdbqt")),
    ]

    audit_rows = []
    for name, lig_qt, ref_pose in targets:
        if not os.path.exists(lig_qt):
            matches = glob.glob(os.path.join(LIG_PH74_DIR, f"*{name.split()[0]}*.pdbqt"))
            lig_qt = matches[0] if matches else lig_qt
        safe = name.split()[0].replace("-", "")

        # Seed 42 Conformer A (re-dock deterministically to ensure identical binary & box)
        out_a42 = os.path.join(tmp_dir, f"{safe}_altA_seed42.pdbqt")
        m_a42 = shift_and_dock_vina(VINA_EXE, rec_a, lig_qt, out_a42, box_4a5s, exhaustiveness=8, seed=42, cpu=4)
        dg_a42 = m_a42[0]["affinity"]

        out_a137 = os.path.join(tmp_dir, f"{safe}_altA_seed137.pdbqt")
        m_a137 = shift_and_dock_vina(VINA_EXE, rec_a, lig_qt, out_a137, box_4a5s, exhaustiveness=8, seed=137, cpu=4)
        dg_a137 = m_a137[0]["affinity"]

        out_a2026 = os.path.join(tmp_dir, f"{safe}_altA_seed2026.pdbqt")
        m_a2026 = shift_and_dock_vina(VINA_EXE, rec_a, lig_qt, out_a2026, box_4a5s, exhaustiveness=8, seed=2026, cpu=4)
        dg_a2026 = m_a2026[0]["affinity"]

        # Seed 42 Conformer B (Val546 altLoc B)
        out_b42 = os.path.join(tmp_dir, f"{safe}_altB_seed42.pdbqt")
        m_b42 = shift_and_dock_vina(VINA_EXE, rec_b, lig_qt, out_b42, box_4a5s, exhaustiveness=8, seed=42, cpu=4)
        dg_b42 = m_b42[0]["affinity"]

        rmsd_s137 = hungarian_heavy_atom_rmsd(out_a42, out_a137, is_pdbqt_a=True, is_pdbqt_b=True)
        rmsd_s2026 = hungarian_heavy_atom_rmsd(out_a42, out_a2026, is_pdbqt_a=True, is_pdbqt_b=True)
        rmsd_altb = hungarian_heavy_atom_rmsd(out_a42, out_b42, is_pdbqt_a=True, is_pdbqt_b=True)

        dgs_3 = [dg_a42, dg_a137, dg_a2026]
        row = {
            "Compound": name,
            "AltA_Seed42_dG": round(dg_a42, 3),
            "AltA_Seed137_dG": round(dg_a137, 3),
            "AltA_Seed2026_dG": round(dg_a2026, 3),
            "Seed_Mean_dG": round(float(np.mean(dgs_3)), 3),
            "Seed_SD_dG": round(float(np.std(dgs_3, ddof=1)), 3),
            "Max_Seed_RMSD_A": round(max(rmsd_s137, rmsd_s2026), 3),
            "AltB_Seed42_dG": round(dg_b42, 3),
            "AltB_minus_AltA_dG": round(dg_b42 - dg_a42, 3),
            "AltB_vs_AltA_RMSD_A": round(rmsd_altb, 3),
        }
        audit_rows.append(row)
        print(f"   {name:18s}: AltA(3-seed)={row['Seed_Mean_dG']:.3f}±{row['Seed_SD_dG']:.3f} (maxRMSD={row['Max_Seed_RMSD_A']:.3f} Å) | AltB={row['AltB_Seed42_dG']:.3f} (Δ={row['AltB_minus_AltA_dG']:+.3f}, RMSD={row['AltB_vs_AltA_RMSD_A']:.3f} Å)", flush=True)

    df_out = pd.DataFrame(audit_rows)
    out_csv = os.path.join(RESULTS_DIR, "protocol_b_altloc_and_multiseed_audit.csv")
    df_out.to_csv(out_csv, index=False)
    # Also save fixed-pose Val546 audit
    pd.DataFrame(fixed_rows).to_csv(os.path.join(RESULTS_DIR, "protocol_b_val546_fixed_pose_audit.csv"), index=False)
    print(f"[STEP 3b] Saved 4A5S AltLoc B + 3-Seed Audit: {out_csv}")
    return df_out


if __name__ == "__main__":
    step1_update_protocol_b_boiled_egg()
    step2_audit_dpp4_roc_benchmark()
    step3_4a5s_altloc_and_multiseed()
