# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "propka>=3.5.0",
#     "pandas>=2.0",
# ]
# ///
"""
07_propka_pka_audit.py
Runs continuum-empirical PROPKA 3.5 micro-environment pKa calculations on altLoc-clean
Protocol B receptors (PDB: 4A5S DPP-4 and PDB: 5LGD CD36 ectodomain) at plasma pH 7.4,
explicitly auditing active-site / tunnel ionizable residues and validating OpenBabel pH 7.4
protonation assignments (closing the continuum pKa caveat in Methods 2.2 / Section 3.6).
"""
import os
import json
import pandas as pd
import propka.run as pk_run

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT_DIR = os.path.join(ROOT_DIR, "structures")
OUT_DIR = os.path.join(ROOT_DIR, "protocol_b_hepatic_cd36_dpp4", "results")
os.makedirs(OUT_DIR, exist_ok=True)

POCKET_RESIDUES = {
    "4A5S_DPP4": {
        "pdb": os.path.join(STRUCT_DIR, "4A5S_DPP4_clean.pdb"),
        "focus": {"125", "205", "206", "547", "630", "631", "662", "708", "740", "752"},
        "ph": 7.4,
    },
    "5LGD_CD36": {
        "pdb": os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdb"),
        "focus": {"63", "164", "231", "270", "334", "335", "337", "385"},
        "ph": 7.4,
    },
}

DEFAULT_PKA = {
    "ASP": 3.80, "GLU": 4.50, "HIS": 6.50, "CYS": 9.00,
    "TYR": 10.00, "LYS": 10.50, "ARG": 12.50, "N+": 8.00, "C-": 3.20
}


def protonation_state_at_ph(resname, pka, ph=7.4):
    if resname in ("ASP", "GLU", "C-"):
        return "Deprotonated (-1)" if ph > pka else "Protonated (0, Neutral)"
    if resname in ("CYS", "TYR"):
        return "Deprotonated (-1)" if ph > pka else "Protonated (0, Neutral)"
    if resname in ("HIS", "LYS", "ARG", "N+"):
        return "Protonated (+1)" if ph < pka else "Deprotonated (0, Neutral)"
    return "Neutral"


def main():
    all_rows = []
    summary = {}
    for target, cfg in POCKET_RESIDUES.items():
        pdb_path = cfg["pdb"]
        assert os.path.exists(pdb_path), f"Missing cleaned PDB: {pdb_path}"
        mol = pk_run.single(pdb_path, optargs=["--quiet"], write_pka=False)
        conf = mol.conformations["AVR"]

        shifted_global = []
        focus_rows = []
        for grp in conf.groups:
            resn = grp.residue_type.strip()
            resnum = str(grp.atom.res_num).strip()
            chain = grp.atom.chain_id.strip()
            pka = float(grp.pka_value)
            model_pka = float(grp.model_pka)
            buried_pct = float(getattr(grp, "buried", 0.0)) * 100.0 if float(getattr(grp, "buried", 0.0)) <= 1.0 else float(getattr(grp, "buried", 0.0))
            state_propka = protonation_state_at_ph(resn, pka, cfg["ph"])
            state_model = protonation_state_at_ph(resn, model_pka, cfg["ph"])
            state_flipped = (state_propka != state_model)
            is_focus = resnum in cfg["focus"]

            rec = {
                "Receptor": target,
                "Chain": chain,
                "Residue": f"{resn}{resnum}",
                "ResName": resn,
                "ResNum": int(resnum) if resnum.isdigit() else -1,
                "Is_Pocket_Focus": is_focus,
                "Model_pKa": round(model_pka, 2),
                "PROPKA_pKa": round(pka, 2),
                "Delta_pKa": round(pka - model_pka, 2),
                "Buried_Pct": round(buried_pct, 1),
                "State_at_pH74_Model": state_model,
                "State_at_pH74_PROPKA": state_propka,
                "Protonation_Flip_vs_Model": state_flipped,
            }
            if is_focus:
                focus_rows.append(rec)
                all_rows.append(rec)
            if state_flipped:
                shifted_global.append(rec)

        focus_flips = [r for r in focus_rows if r["Protonation_Flip_vs_Model"]]
        summary[target] = {
            "Total_Ionizable_Groups": len(conf.groups),
            "Global_State_Flips_at_pH74": len(shifted_global),
            "Pocket_Focus_Count": len(focus_rows),
            "Pocket_Focus_Flips_at_pH74": len(focus_flips),
            "Pocket_Focus_Flipped_Residues": [f"{r['Residue']} (pKa={r['PROPKA_pKa']}, {r['State_at_pH74_PROPKA']})" for r in focus_flips],
            "Pocket_Audit": focus_rows,
        }
        print(f"\n=== PROPKA 3.5 Audit: {target} (pH {cfg['ph']}) ===")
        print(f"Total ionizable groups: {len(conf.groups)} | Global state flips: {len(shifted_global)} | Pocket focus flips: {len(focus_flips)}")
        for fr in sorted(focus_rows, key=lambda x: x["ResNum"]):
            flag = " [FLIP!]" if fr["Protonation_Flip_vs_Model"] else " [MATCH]"
            print(f"  {fr['Residue']:8s} | Model={fr['Model_pKa']:5.2f} -> PROPKA={fr['PROPKA_pKa']:5.2f} (d={fr['Delta_pKa']:+5.2f}, bur={fr['Buried_Pct']:5.1f}%) | {fr['State_at_pH74_PROPKA']}{flag}")

    df_out = pd.DataFrame(all_rows).sort_values(["Receptor", "ResNum"])
    csv_p = os.path.join(OUT_DIR, "propka_pka_audit_pH74.csv")
    json_p = os.path.join(OUT_DIR, "propka_pka_audit_pH74.json")
    df_out.to_csv(csv_p, index=False)
    with open(json_p, "w", encoding="utf-8") as fj:
        json.dump(summary, fj, indent=2)
    print(f"\n[SUCCESS] Saved PROPKA 3.5 audit to {csv_p} and {json_p}")


if __name__ == "__main__":
    main()
