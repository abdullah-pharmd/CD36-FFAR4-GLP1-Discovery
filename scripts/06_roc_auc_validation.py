"""
06_roc_auc_validation.py
Empirical Benchmarking of AutoDock Vina Protocols Against Size-Balanced ChEMBL
Experimental Actives vs. Confirmed Inactives:
  1. Human FFAR4 (PDB: 8G59, pH 6.8) vs ChEMBL5339 (15 Actives pChEMBL >= 6.5 vs 15 Inactives > 10 uM)
  2. Human DPP-4 (PDB: 4A5S, pH 7.4) vs ChEMBL284 (30 Non-Covalent Actives pChEMBL >= 7.0 vs 30 Confirmed Inactives >= 50 uM)
Computes ROC-AUC, Enrichment Factor (EF20%), and Mann-Whitney U p-values for both raw ΔG and LE.
"""
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ttest_ind
from sklearn.metrics import roc_auc_score, roc_curve
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structural_prep_utils import prepare_ligand_pdbqt, shift_and_dock_vina

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VINA_EXE = os.path.join(BASE_DIR, "scripts", "vina.exe")
BENCH_DIR = os.path.join(BASE_DIR, "results", "benchmark_poses")
os.makedirs(BENCH_DIR, exist_ok=True)


def run_benchmark(
    pool_json, rec_pdbqt, box, ph, out_csv, out_fig,
    title_label, short_label, act_thresh_lbl, inact_thresh_lbl,
    exhaustiveness=6, force_redock=False, source_prefix="ChEMBL"
):
    with open(pool_json) as f:
        pool = json.load(f)

    items = [(1, x) for x in pool["actives"]] + [(0, x) for x in pool["inactives"]]
    rows = []
    for idx, (label, item) in enumerate(items, 1):
        cid = item["chembl_id"]
        smi = item["smiles"]
        m_base = Chem.MolFromSmiles(smi)
        if m_base is None:
            continue
        n_heavy = m_base.GetNumHeavyAtoms()
        mol = Chem.AddHs(m_base)
        etkdg = AllChem.ETKDGv3()
        etkdg.randomSeed = 42
        AllChem.EmbedMolecule(mol, etkdg)
        AllChem.MMFFOptimizeMolecule(mol)

        sdf_p = os.path.join(BENCH_DIR, f"{title_label[:4]}_{cid}.sdf")
        lig_qt = os.path.join(BENCH_DIR, f"{title_label[:4]}_{cid}_pH{ph}.pdbqt")
        pose_qt = os.path.join(BENCH_DIR, f"{title_label[:4]}_{cid}_out.pdbqt")
        w = Chem.SDWriter(sdf_p)
        w.write(mol)
        w.close()

        prepare_ligand_pdbqt(sdf_p, lig_qt, ph=ph)
        if (not force_redock) and os.path.exists(pose_qt) and os.path.getsize(pose_qt) > 200:
            aff = np.nan
            with open(pose_qt, "r", encoding="utf-8", errors="ignore") as fp:
                for line in fp:
                    if "REMARK VINA RESULT:" in line:
                        aff = float(line.split()[3])
                        break
        else:
            modes = shift_and_dock_vina(VINA_EXE, rec_pdbqt, lig_qt, pose_qt, box, exhaustiveness=exhaustiveness, seed=42, cpu=4)
            aff = modes[0]["affinity"] if modes else np.nan
        le = round(-aff / n_heavy, 3) if (not np.isnan(aff) and n_heavy > 0) else np.nan
        rows.append({
            "Target": title_label, "ChEMBL_ID": cid, "True_Active": label,
            "Experimental_Val": item["val"], "MW_Da": item["mw"], "Heavy_Atoms": n_heavy,
            "Vina_Affinity_kcal_mol": aff, "Ligand_Efficiency": le
        })
        print(f"  [{idx:02d}/{len(items)}] {cid} ({'ACTIVE' if label else 'INACT '}, MW={item['mw']}): dG={aff:.2f} | LE={le:.3f}", flush=True)

    df = pd.DataFrame(rows).dropna()
    df.to_csv(out_csv, index=False)

    y = df["True_Active"].values
    scores_aff = -df["Vina_Affinity_kcal_mol"].values
    scores_le = df["Ligand_Efficiency"].values

    auc_aff = roc_auc_score(y, scores_aff)
    auc_le = roc_auc_score(y, scores_le)
    act_df = df[df["True_Active"] == 1]
    inact_df = df[df["True_Active"] == 0]
    _, p_mw = ttest_ind(act_df["MW_Da"], inact_df["MW_Da"])
    _, p_a = mannwhitneyu(act_df["Vina_Affinity_kcal_mol"], inact_df["Vina_Affinity_kcal_mol"], alternative="less")
    _, p_l = mannwhitneyu(act_df["Ligand_Efficiency"], inact_df["Ligand_Efficiency"], alternative="greater")

    top_n = max(1, int(0.20 * len(df)))
    ef20_aff = (df.sort_values("Vina_Affinity_kcal_mol", ascending=True).head(top_n)["True_Active"].sum() / top_n) / (df["True_Active"].mean())
    ef20_le = (df.sort_values("Ligand_Efficiency", ascending=False).head(top_n)["True_Active"].sum() / top_n) / (df["True_Active"].mean())

    metrics = {
        "Target": title_label,
        "N_Actives": int(df["True_Active"].sum()),
        "N_Inactives": int((df["True_Active"] == 0).sum()),
        "Mean_Active_MW": round(float(act_df["MW_Da"].mean()), 1),
        "Mean_Inactive_MW": round(float(inact_df["MW_Da"].mean()), 1),
        "MW_Balance_p": round(float(p_mw), 3),
        "Mean_Active_dG": round(float(act_df["Vina_Affinity_kcal_mol"].mean()), 2),
        "Mean_Inactive_dG": round(float(inact_df["Vina_Affinity_kcal_mol"].mean()), 2),
        "Mean_Active_LE": round(float(act_df["Ligand_Efficiency"].mean()), 3),
        "Mean_Inactive_LE": round(float(inact_df["Ligand_Efficiency"].mean()), 3),
        "ROC_AUC_Affinity": round(float(auc_aff), 3),
        "ROC_AUC_LE": round(float(auc_le), 3),
        "EF_20pct": round(float(ef20_aff), 2),
        "EF_20pct_LE": round(float(ef20_le), 2),
        "MannWhitney_p_Affinity": float(f"{p_a:.4e}"),
        "MannWhitney_p_LE": float(f"{p_l:.4e}"),
        "Exhaustiveness": exhaustiveness
    }
    print(f"\n[{title_label} BENCHMARK SUMMARY] {json.dumps(metrics)}", flush=True)

    # 3-Panel Publication ROC + Raw dG Boxplot + LE Boxplot (Zero Title Collisions!)
    fpr, tpr, _ = roc_curve(y, scores_aff)
    fpr_le, tpr_le, _ = roc_curve(y, scores_le)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=300)

    ax1.plot(fpr, tpr, color="#0B3C5D", lw=2.3, label=f"Vina ΔG (AUC = {auc_aff:.3f}, p = {p_a:.3g})")
    ax1.plot(fpr_le, tpr_le, color="#16A34A", lw=2.3, linestyle="-.", label=f"Ligand Efficiency LE (AUC = {auc_le:.3f}, p = {p_l:.3g})")
    ax1.plot([0, 1], [0, 1], "k--", lw=1.2, alpha=0.6, label="Random Baseline (AUC = 0.500)")
    ax1.set_xlabel("False Positive Rate (1 − Specificity)", fontweight="bold")
    ax1.set_ylabel("True Positive Rate (Sensitivity)", fontweight="bold")
    ax1.set_title(f"A. Active vs. Inactive ROC\n({short_label}; N={len(df)}, p_MW={p_mw:.2f})", fontweight="bold", fontsize=10.5)
    ax1.legend(loc="lower right", fontsize=8.2, framealpha=0.94)
    ax1.grid(True, linestyle=":", alpha=0.5)

    x_lbls = [f"{source_prefix} Actives\n({act_thresh_lbl})", f"{source_prefix} Inactives\n({inact_thresh_lbl})"]

    act_dgs = act_df["Vina_Affinity_kcal_mol"]
    inact_dgs = inact_df["Vina_Affinity_kcal_mol"]
    bp = ax2.boxplot([act_dgs, inact_dgs], tick_labels=x_lbls, patch_artist=True, widths=0.5)
    bp["boxes"][0].set_facecolor("#BBF7D0")
    bp["boxes"][1].set_facecolor("#FECACA")
    for i_p, grp in enumerate([act_dgs, inact_dgs], 1):
        jitter = np.random.RandomState(42).normal(i_p, 0.05, size=len(grp))
        ax2.scatter(jitter, grp, color="#1E293B", s=26, alpha=0.82, zorder=3)
    ax2.set_ylabel("AutoDock Vina Binding Energy ΔG (kcal/mol)", fontweight="bold")
    ax2.set_title(f"B. Raw Binding Energy Separation\n(EF20% = {ef20_aff:.2f}x, Mann–Whitney p = {p_a:.3g})", fontweight="bold", fontsize=10.5)
    ax2.grid(axis="y", linestyle=":", alpha=0.5)

    act_les = act_df["Ligand_Efficiency"]
    inact_les = inact_df["Ligand_Efficiency"]
    bp3 = ax3.boxplot([act_les, inact_les], tick_labels=x_lbls, patch_artist=True, widths=0.5)
    bp3["boxes"][0].set_facecolor("#BBF7D0")
    bp3["boxes"][1].set_facecolor("#FECACA")
    for i_p, grp in enumerate([act_les, inact_les], 1):
        jitter = np.random.RandomState(42).normal(i_p, 0.05, size=len(grp))
        ax3.scatter(jitter, grp, color="#1E293B", s=26, alpha=0.82, zorder=3)
    ax3.set_ylabel("Ligand Efficiency (LE = −ΔG / N_heavy)", fontweight="bold")
    ax3.set_title(f"C. Ligand Efficiency (LE) Separation\n(EF20%_LE = {ef20_le:.2f}x, Mann–Whitney p = {p_l:.3g})", fontweight="bold", fontsize=10.5)
    ax3.grid(axis="y", linestyle=":", alpha=0.5)

    plt.tight_layout(w_pad=2.4)
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.close()
    return metrics


def main():
    os.makedirs(os.path.join(BASE_DIR, "results", "figures"), exist_ok=True)
    os.makedirs(os.path.join(BASE_DIR, "protocol_b_hepatic_cd36_dpp4", "figures"), exist_ok=True)
    force_dpp4 = "--refresh-dpp4" in sys.argv
    dpp4_only = "--dpp4-only" in sys.argv

    box_8g59 = {"center_x": 98.09, "center_y": 123.78, "center_z": 117.74, "size_x": 22.0, "size_y": 22.0, "size_z": 22.0}
    box_4a5s = {"center_x": 14.14, "center_y": 35.87, "center_z": 54.84, "size_x": 22.0, "size_y": 22.0, "size_z": 22.0}

    summary_path = os.path.join(BASE_DIR, "results", "roc_summary_all.json")
    summary = {}
    if os.path.exists(summary_path):
        with open(summary_path) as f:
            summary = json.load(f)

    if not dpp4_only:
        m_ffar4 = run_benchmark(
            os.path.join(BASE_DIR, "results", "network", "ffar4_chembl_benchmark_pool.json"),
            os.path.join(BASE_DIR, "structures", "8G59_FFAR4_clean.pdbqt"),
            box_8g59, 6.8,
            os.path.join(BASE_DIR, "results", "roc_validation_ffar4_chembl5339.csv"),
            os.path.join(BASE_DIR, "results", "figures", "fig_roc_validation_ffar4.png"),
            "FFAR4_8G59_CHEMBL5339", "FFAR4 8G59 / ChEMBL5339",
            "pChEMBL ≥ 6.5", "> 10,000 nM",
            exhaustiveness=6, force_redock=False
        )
        summary["FFAR4"] = m_ffar4

    m_dpp4 = run_benchmark(
        os.path.join(BASE_DIR, "protocol_b_hepatic_cd36_dpp4", "results", "network", "dpp4_chembl_benchmark_pool.json"),
        os.path.join(BASE_DIR, "structures", "4A5S_DPP4_pH74.pdbqt"),
        box_4a5s, 7.4,
        os.path.join(BASE_DIR, "protocol_b_hepatic_cd36_dpp4", "results", "roc_validation_dpp4_chembl284.csv"),
        os.path.join(BASE_DIR, "protocol_b_hepatic_cd36_dpp4", "figures", "fig_roc_validation_dpp4.png"),
        "DPP4_4A5S_PUBCHEM1803", "DPP-4 4A5S / PubChem GeneID:1803",
        "≤ 100 nM", "≥ 50 µM / HTS-Inact",
        exhaustiveness=8, force_redock=force_dpp4, source_prefix="PubChem"
    )
    summary["DPP4"] = m_dpp4
    with open(summary_path, "w") as fj:
        json.dump(summary, fj, indent=2)


if __name__ == "__main__":
    main()
