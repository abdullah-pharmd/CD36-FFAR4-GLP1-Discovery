import os, sys, json, time, requests
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NET_DIR = os.path.join(BASE_DIR, "results", "network")
PROT_B_NET = os.path.join(BASE_DIR, "protocol_b_hepatic_cd36_dpp4", "results", "network")
os.makedirs(NET_DIR, exist_ok=True)
os.makedirs(PROT_B_NET, exist_ok=True)

HDR = {"Accept": "application/json"}
COVALENT_SMARTS = [
    Chem.MolFromSmarts(s)
    for s in ["[CX2]#[NX1]", "[BX3]", "[N-]=[N+]=[N-]", "C(=O)C=C", "C1OC1", "C=CC(=O)"]
]


def get_json_retry(url, timeout=45):
    for attempt in range(4):
        try:
            r = requests.get(url, headers=HDR, timeout=timeout)
            if r.status_code == 200:
                return r.json()
        except Exception:
            pass
        time.sleep(1.5 * (attempt + 1))
    return {}


def is_clean_noncovalent(smi):
    if not smi or "." in smi:
        return None
    m = Chem.MolFromSmiles(smi)
    if not m:
        return None
    for pat in COVALENT_SMARTS:
        if pat and m.HasSubstructMatch(pat):
            return None
    return m


def fetch_opentargets():
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
        res = requests.post(
            "https://api.platform.opentargets.org/api/v4/graphql",
            json={"query": qt, "variables": {"id": did, "size": 300}},
            timeout=25
        ).json()
        d = res["data"]["disease"]
        for row in d["associatedTargets"]["rows"]:
            ot_rows.append({
                "Disease_ID": did,
                "Disease_Name": d["name"],
                "Gene_Symbol": row["target"]["approvedSymbol"],
                "OpenTargets_Score": round(float(row["score"]), 4)
            })
    df_ot = pd.DataFrame(ot_rows)
    df_ot.to_csv(os.path.join(NET_DIR, "opentargets_disease_genes_live.csv"), index=False)
    df_ot.to_csv(os.path.join(PROT_B_NET, "opentargets_disease_genes_live.csv"), index=False)
    print(f"[OpenTargets Live] Pulled {len(df_ot)} associations ({df_ot['Gene_Symbol'].nunique()} unique HGNC symbols).", flush=True)


def build_mw_balanced_benchmark(
    target_id, pmin_active, mw_low=260, mw_high=460, n_each=15,
    min_inact_nm=10000, max_pages=3, strict_noncov=False, strict_units=False
):
    def fetch_pool(start_url, is_act):
        pool, seen_ik = [], set()
        url = start_url
        for _ in range(max_pages):
            js = get_json_retry(url, timeout=45)
            acts = js.get("activities", [])
            if not acts:
                break
            for a in acts:
                if strict_units:
                    if a.get("target_organism") != "Homo sapiens" or a.get("standard_units") != "nM":
                        continue
                    if a.get("standard_type") not in ("IC50", "Ki", "Kd"):
                        continue
                    if not is_act:
                        if a.get("standard_relation") not in (">", ">=", "="):
                            continue
                        try:
                            if float(a.get("standard_value")) < min_inact_nm:
                                continue
                        except Exception:
                            continue
                smi = a.get("canonical_smiles")
                cid = a.get("molecule_chembl_id")
                m = is_clean_noncovalent(smi) if strict_noncov else (Chem.MolFromSmiles(smi) if (smi and "." not in smi) else None)
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
                pool.append({
                    "chembl_id": cid, "smiles": smi, "mw": round(mw, 2), "rotb": rotb,
                    "ik14": ik14, "val": a.get("pchembl_value") if is_act else a.get("standard_value")
                })
            nxt = js.get("page_meta", {}).get("next")
            if not nxt or len(pool) >= 250:
                break
            url = "https://www.ebi.ac.uk" + nxt
        return pool

    act_pool = fetch_pool(
        f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={target_id}&pchembl_value__gte={pmin_active}&limit=250",
        True
    )
    inact_pool = fetch_pool(
        f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={target_id}&standard_value__gte=30000&limit=500",
        False
    )
    act_iks = {x["ik14"] for x in act_pool}
    inact_pool = [x for x in inact_pool if x["ik14"] not in act_iks]

    act_sorted = sorted(act_pool, key=lambda x: x["mw"])
    matched_a, matched_i = [], []
    used_i = set()
    for a in act_sorted:
        best_j, best_diff = None, 999
        for j, i_item in enumerate(inact_pool):
            if j in used_i:
                continue
            diff = abs(a["mw"] - i_item["mw"])
            if diff < best_diff:
                best_diff, best_j = diff, j
        if best_j is not None and best_diff <= 20.0:
            used_i.add(best_j)
            matched_a.append(a)
            matched_i.append(inact_pool[best_j])
        if len(matched_a) >= n_each:
            break

    _, p_mw = ttest_ind([x["mw"] for x in matched_a], [x["mw"] for x in matched_i])
    print(
        f"[{target_id} Benchmark] Matched={len(matched_a)} vs {len(matched_i)} | "
        f"Act MW={np.mean([x['mw'] for x in matched_a]):.1f} vs Inact MW={np.mean([x['mw'] for x in matched_i]):.1f} (p_MW={p_mw:.3f})",
        flush=True
    )
    return matched_a, matched_i


def build_pubchem_dpp4_fallback(n_each=30):
    import io
    r = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/geneid/1803/concise/CSV", timeout=60)
    df = pd.read_csv(io.BytesIO(r.content), low_memory=False)
    df = df[df["CID"].notna()].copy()
    df["CID"] = df["CID"].astype(int)
    df["Target GeneID"] = pd.to_numeric(df["Target GeneID"], errors="coerce")
    df = df[df["Target GeneID"] == 1803]
    if "RNAi" in df.columns:
        df = df[df["RNAi"].isna() | (df["RNAi"].astype(str) == "0") | (df["RNAi"].astype(str) == "")]
    df["val_uM"] = pd.to_numeric(df["Activity Value [uM]"], errors="coerce")
    active_or_weak_cids = set(df[(df["Activity Outcome"] == "Active") | ((df["val_uM"].notna()) & (df["val_uM"] < 15.0))]["CID"].unique())
    act_df = df[(df["Activity Outcome"] == "Active") & (df["Activity Name"].isin(["IC50", "Ki", "Kd"])) & (df["val_uM"] > 0) & (df["val_uM"] <= 0.1)].copy()
    inact_df = df[(~df["CID"].isin(active_or_weak_cids)) & ((df["Activity Outcome"] == "Inactive") | ((df["Activity Name"].isin(["IC50", "Ki", "Kd"])) & (df["val_uM"] >= 50.0)))].copy()
    act_cids = act_df.groupby("CID")["val_uM"].median().sort_values().index.tolist()[:600]
    inact_quant = inact_df[inact_df["val_uM"] >= 50.0]["CID"].unique().tolist()
    inact_conf = [c for c in inact_df["CID"].unique().tolist() if c not in set(inact_quant)]
    inact_cids = (inact_quant + inact_conf)[:600]

    def fetch_smi(cids):
        out = {}
        for i in range(0, len(cids), 100):
            chunk = ",".join(str(x) for x in cids[i:i + 100])
            js = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{chunk}/property/SMILES,CanonicalSMILES,ConnectivitySMILES/JSON", timeout=30).json()
            for p in js.get("PropertyTable", {}).get("Properties", []):
                s = p.get("SMILES") or p.get("CanonicalSMILES") or p.get("ConnectivitySMILES")
                if s:
                    out[int(p["CID"])] = s
        return out

    smi_a, smi_i = fetch_smi(act_cids), fetch_smi(inact_cids)
    act_vals = {int(c): round(-np.log10(v * 1e-6), 2) for c, v in act_df.groupby("CID")["val_uM"].median().items()}
    inact_vals = {int(c): (round(float(v) * 1000.0, 1) if pd.notna(v) else "Inactive_HTS") for c, v in inact_df.groupby("CID")["val_uM"].max().items()}

    def clean_p(cids, smap, vmap):
        out, seen = [], set()
        for c in cids:
            m = is_clean_noncovalent(smap.get(c))
            if not m:
                continue
            mw, rotb = Descriptors.MolWt(m), Lipinski.NumRotatableBonds(m)
            if not (250 <= mw <= 460 and rotb <= 7):
                continue
            ik14 = Chem.MolToInchiKey(m)[:14]
            if ik14 in seen:
                continue
            seen.add(ik14)
            out.append({"chembl_id": f"CID{c}", "cid": c, "smiles": smap[c], "mw": round(mw, 2), "rotb": rotb, "ik14": ik14, "val": vmap.get(c, "Inactive")})
        return out

    pa, pi = clean_p(act_cids, smi_a, act_vals), clean_p(inact_cids, smi_i, inact_vals)
    pi = [x for x in pi if x["ik14"] not in {a["ik14"] for a in pa}]
    ma, mi, used_i = [], [], set()
    for a in sorted([x for x in pa if x["mw"] >= 280], key=lambda x: x["mw"]):
        best_j, best_diff = None, 999
        for j, it in enumerate(pi):
            if j in used_i:
                continue
            d = abs(a["mw"] - it["mw"])
            if d < best_diff:
                best_diff, best_j = d, j
        if best_j is not None and best_diff <= 15.0:
            used_i.add(best_j)
            ma.append(a)
            mi.append(pi[best_j])
        if len(ma) >= n_each:
            break
    return ma, mi


def main():
    if "--dpp4-only" not in sys.argv:
        fetch_opentargets()
        fa_act, fa_inact = build_mw_balanced_benchmark("CHEMBL5339", 6.5, mw_low=260, mw_high=440, n_each=15, min_inact_nm=10000)
        if fa_act and fa_inact:
            with open(os.path.join(NET_DIR, "ffar4_chembl_benchmark_pool.json"), "w") as f:
                json.dump({"actives": fa_act, "inactives": fa_inact}, f, indent=2)

    # Expanded N=60 (30 actives <= 100 nM vs 30 confirmed non-binders >= 50 uM / HTS-inactive) strict non-covalent DPP-4 benchmark
    dp_act, dp_inact = build_mw_balanced_benchmark(
        "CHEMBL284", 7.0, mw_low=300, mw_high=460, n_each=30,
        min_inact_nm=50000, max_pages=6, strict_noncov=True, strict_units=True
    )
    source_tag = "ChEMBL284"
    if not (dp_act and dp_inact):
        print("[INFO] EBI ChEMBL unavailable/empty; falling back to live PubChem BioAssay DPP4 (GeneID: 1803)...", flush=True)
        dp_act, dp_inact = build_pubchem_dpp4_fallback(n_each=30)
        source_tag = "PubChem_BioAssay_GeneID1803"
    if dp_act and dp_inact:
        with open(os.path.join(PROT_B_NET, "dpp4_chembl_benchmark_pool.json"), "w") as f:
            json.dump({"source": source_tag, "actives": dp_act, "inactives": dp_inact}, f, indent=2)
    else:
        print("[WARN] Skipped writing dpp4_chembl_benchmark_pool.json due to empty API response.")


if __name__ == "__main__":
    main()
