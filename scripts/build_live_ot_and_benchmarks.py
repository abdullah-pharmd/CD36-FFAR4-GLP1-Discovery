import os, sys, json, time, requests
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, inchi

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NET_DIR = os.path.join(BASE_DIR, "results", "network")
PROT_B_NET = os.path.join(BASE_DIR, "protocol_b_hepatic_cd36_dpp4", "results", "network")
os.makedirs(NET_DIR, exist_ok=True)
os.makedirs(PROT_B_NET, exist_ok=True)

# Step 1: Query OpenTargets GraphQL for 4 standardized phenotypes
qt = """
query GetTargets($id: String!, $size: Int!) {
  disease(efoId: $id) {
    id
    name
    associatedTargets(page: {index: 0, size: $size}) {
      rows { target { approvedSymbol } score }
    }
  }
}
"""
ot_rows = []
for did in ["MONDO_0011122", "MONDO_0005148", "EFO_1001248", "EFO_1000653"]:
    res = requests.post("https://api.platform.opentargets.org/api/v4/graphql", json={"query": qt, "variables": {"id": did, "size": 300}}, timeout=25).json()
    d = res["data"]["disease"]
    for row in d["associatedTargets"]["rows"]:
        ot_rows.append({"Disease_ID": did, "Disease_Name": d["name"], "Gene_Symbol": row["target"]["approvedSymbol"], "OpenTargets_Score": round(float(row["score"]), 4)})
df_ot = pd.DataFrame(ot_rows)
df_ot.to_csv(os.path.join(NET_DIR, "opentargets_disease_genes_live.csv"), index=False)
df_ot.to_csv(os.path.join(PROT_B_NET, "opentargets_disease_genes_live.csv"), index=False)
ot_gene_set = set(df_ot["Gene_Symbol"])
print(f"[OpenTargets Live] Pulled {len(df_ot)} associations across 4 phenotypes ({len(ot_gene_set)} unique HGNC symbols).", flush=True)

# Step 2: Check MW-matched ChEMBL Active vs Inactive pools for FFAR4 (CHEMBL5339) and DPP4 (CHEMBL284)
def build_mw_balanced_benchmark(target_id, pmin_active, mw_low=260, mw_high=460, n_each=15):
    def fetch_pool(url):
        pool, seen_ik = [], set()
        r = requests.get(url, timeout=25).json()
        for a in r.get("activities", []):
            smi = a.get("canonical_smiles")
            cid = a.get("molecule_chembl_id")
            if not smi or "." in smi:
                continue
            m = Chem.MolFromSmiles(smi)
            if not m:
                continue
            mw = Descriptors.MolWt(m)
            rotb = Lipinski.NumRotatableBonds(m)
            if not (mw_low <= mw <= mw_high and rotb <= 7):
                continue
            ik14 = Chem.MolToInchiKey(m)[:14]
            if ik14 in seen_ik:
                continue
            seen_ik.add(ik14)
            pool.append({"chembl_id": cid, "smiles": smi, "mw": round(mw, 2), "rotb": rotb, "ik14": ik14, "val": a.get("pchembl_value") or a.get("standard_value")})
        return pool

    act_pool = fetch_pool(f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={target_id}&pchembl_value__gte={pmin_active}&limit=300")
    inact_pool = fetch_pool(f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={target_id}&standard_relation__in=>,>=&standard_value__gte=10000&limit=400")
    act_iks = {x["ik14"] for x in act_pool}
    inact_pool = [x for x in inact_pool if x["ik14"] not in act_iks]

    # Nearest-MW greedy 1-to-1 pairing so Active & Inactive MW distributions are balanced (zero size bias!)
    act_sorted = sorted(act_pool, key=lambda x: x["mw"])
    matched_a, matched_i = [], []
    used_i = set()
    # Try middle/diverse actives to avoid extreme edge MWs
    for a in act_sorted:
        best_j, best_diff = None, 999
        for j, i_item in enumerate(inact_pool):
            if j in used_i:
                continue
            diff = abs(a["mw"] - i_item["mw"])
            if diff < best_diff:
                best_diff, best_j = diff, j
        if best_j is not None and best_diff <= 25.0:
            used_i.add(best_j)
            matched_a.append(a)
            matched_i.append(inact_pool[best_j])
        if len(matched_a) >= n_each:
            break

    t_stat, p_mw = ttest_ind([x["mw"] for x in matched_a], [x["mw"] for x in matched_i])
    print(f"[{target_id} Benchmark] Matched pairs={len(matched_a)} | Mean Act MW={np.mean([x['mw'] for x in matched_a]):.1f} vs Inact MW={np.mean([x['mw'] for x in matched_i]):.1f} (MW balance p={p_mw:.3f})", flush=True)
    return matched_a, matched_i

fa_act, fa_inact = build_mw_balanced_benchmark("CHEMBL5339", 6.5, mw_low=260, mw_high=440, n_each=15)
dp_act, dp_inact = build_mw_balanced_benchmark("CHEMBL284", 7.0, mw_low=250, mw_high=460, n_each=15)

with open(os.path.join(NET_DIR, "ffar4_chembl_benchmark_pool.json"), "w") as f:
    json.dump({"actives": fa_act, "inactives": fa_inact}, f, indent=2)
with open(os.path.join(PROT_B_NET, "dpp4_chembl_benchmark_pool.json"), "w") as f:
    json.dump({"actives": dp_act, "inactives": dp_inact}, f, indent=2)
