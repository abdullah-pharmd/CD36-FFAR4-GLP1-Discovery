"""
05_live_chembl_string_network.py
100% Live API Network Pharmacology for BOTH Protocol A & Protocol B (ZERO hardcoded target dicts or FDRs!):
  1. Resolves every study ligand's exact InChIKey to ChEMBL Molecule ID (fast O(1) lookup + fallback).
  2. Queries ChEMBL REST API for strictly verified Homo sapiens SINGLE PROTEIN targets
     (standard_units='nM', standard_type in IC50/EC50/Ki/Kd/AC50, relation in =/</<=, value <= 10000 nM)
     plus approved clinical mechanisms (mechanism.json).
  3. Batch-resolves ChEMBL Target IDs -> Official HGNC Gene Symbols (SINGLE PROTEIN only).
  4. Intersects empirical ChEMBL compound targets against OpenTargets GraphQL disease genes (894 unique symbols).
  5. Quantifies empirical ChEMBL target coverage gaps for CD36, FFAR4, and DPP4.
  6. Runs TWO strictly separated STRING v12.0 stages (required_score = 400):
       - Stage 1 (empirical): 100% of pure empirical ChEMBL ∩ OpenTargets overlap genes (ZERO spiked nodes, ZERO head(20) truncation).
       - Stage 2 (bridge): Empirical overlap + pre-specified mechanistic effectors strictly for topological connectivity bridging.
"""
import os
import sys
import json
import time
import requests
import pandas as pd
import numpy as np
import networkx as nx
from rdkit import Chem

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NET_A = os.path.join(BASE_DIR, "results", "network")
NET_B = os.path.join(BASE_DIR, "protocol_b_hepatic_cd36_dpp4", "results", "network")
os.makedirs(NET_A, exist_ok=True)
os.makedirs(NET_B, exist_ok=True)

TARGET_GENE_CACHE = {}
VALID_TYPES = {"IC50", "EC50", "Ki", "Kd", "AC50"}
VALID_RELS = {"=", "<", "<="}


def get_json_retry(url, timeout=35):
    for attempt in range(4):
        try:
            r = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
            if r.status_code == 200:
                return r.json()
        except Exception:
            pass
        time.sleep(1.2 * (attempt + 1))
    return {}


def batch_resolve_single_protein_targets(target_ids):
    missing = [t for t in sorted(set(target_ids)) if t not in TARGET_GENE_CACHE]
    for i in range(0, len(missing), 40):
        chunk = missing[i:i + 40]
        js = get_json_retry(
            f"https://www.ebi.ac.uk/chembl/api/data/target.json?target_chembl_id__in={','.join(chunk)}&limit=100",
            timeout=35
        )
        for t in js.get("targets", []):
            tid = t.get("target_chembl_id")
            if t.get("target_type") == "SINGLE PROTEIN" and t.get("organism") == "Homo sapiens":
                syms = set()
                for comp in t.get("target_components", []):
                    for syn in comp.get("target_component_synonyms", []):
                        if syn.get("syn_type") == "GENE_SYMBOL":
                            syms.add(syn["component_synonym"].upper())
                TARGET_GENE_CACHE[tid] = sorted(syms)
            else:
                TARGET_GENE_CACHE[tid] = []
        for tid in chunk:
            TARGET_GENE_CACHE.setdefault(tid, [])
        time.sleep(0.1)


def filter_strict_quantitative(df):
    """Retain strictly curated quantitative single-protein records (pChEMBL >= 5.0, i.e. <= 10 uM, or Approved Mechanism)."""
    if df.empty:
        return df
    pc_num = pd.to_numeric(df["pChEMBL"], errors="coerce")
    is_app = df["Evidence"].astype(str).str.contains("Approved", case=False, na=False)
    is_strict = ((pc_num >= 5.0) & pc_num.notna()) | is_app
    return df[is_strict].copy()


def query_chembl_targets_for_compounds(comp_df, cache_file, force_refresh=False):
    existing_map = {}
    if os.path.exists(cache_file):
        prev_df = pd.read_csv(cache_file)
        if not force_refresh:
            strict_df = filter_strict_quantitative(prev_df)
            strict_path = cache_file.replace("_live.csv", "_strict_live.csv")
            strict_df.to_csv(strict_path, index=False)
            return strict_df
        for _, pr in prev_df.iterrows():
            existing_map.setdefault(pr["Compound"], set()).add(pr["Molecule_ChEMBL_ID"])

    raw_records = []
    all_tids = set()
    for _, r in comp_df.iterrows():
        cname = r["Compound"]
        smi = r["SMILES"]
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        ik = Chem.MolToInchiKey(mol)
        mcids = sorted(existing_map.get(cname, set()))
        if not mcids:
            # Fast direct O(1) lookup by exact InChIKey first
            direct = get_json_retry(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{ik}.json", timeout=20)
            if direct.get("molecule_chembl_id"):
                mcids = [direct["molecule_chembl_id"]]
            else:
                res_m = get_json_retry(
                    f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key__startswith={ik[:14]}&limit=3",
                    timeout=30
                )
                mcids = [m["molecule_chembl_id"] for m in res_m.get("molecules", [])]
        if not mcids:
            print(f"  [ChEMBL] {cname}: No ChEMBL ID found", flush=True)
            continue

        t_set = {}
        for mcid in mcids:
            mech_r = get_json_retry(
                f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id={mcid}",
                timeout=25
            )
            for m in mech_r.get("mechanisms", []):
                tid = m.get("target_chembl_id")
                if tid:
                    t_set[tid] = {"mcid": mcid, "pchembl": np.nan, "evidence": f"Approved Mechanism ({m.get('action_type', '')})"}

            act_r = get_json_retry(
                f"https://www.ebi.ac.uk/chembl/api/data/activity.json?"
                f"molecule_chembl_id={mcid}&target_organism=Homo+sapiens&standard_units=nM&standard_value__lte=10000&limit=300",
                timeout=35
            )
            for a in act_r.get("activities", []):
                tid = a.get("target_chembl_id")
                val = a.get("standard_value")
                stype = a.get("standard_type")
                srel = a.get("standard_relation")
                pc = a.get("pchembl_value")
                if tid and val is not None and stype in VALID_TYPES and srel in VALID_RELS:
                    try:
                        if float(val) <= 10000.0:
                            if tid not in t_set or (pc and (pd.isna(t_set[tid]["pchembl"]) or float(pc) > float(t_set[tid]["pchembl"]))):
                                t_set[tid] = {"mcid": mcid, "pchembl": pc, "evidence": f"Bioassay <=10uM ({stype}, nM)"}
                    except ValueError:
                        pass

        all_tids.update(t_set.keys())
        raw_records.append((cname, mcids[0], t_set))
        print(f"  [ChEMBL] {cname:32s} ({mcids[0]}): {len(t_set)} raw nM/mechanism targets", flush=True)
        time.sleep(0.05)

    batch_resolve_single_protein_targets(all_tids)
    rows = []
    for cname, _, t_set in raw_records:
        for tid, meta in t_set.items():
            for g in TARGET_GENE_CACHE.get(tid, []):
                rows.append({
                    "Compound": cname,
                    "Molecule_ChEMBL_ID": meta["mcid"],
                    "Target_ChEMBL_ID": tid,
                    "Gene_Symbol": g,
                    "pChEMBL": meta["pchembl"],
                    "Evidence": meta["evidence"]
                })

    df_out = pd.DataFrame(rows, columns=["Compound", "Molecule_ChEMBL_ID", "Target_ChEMBL_ID", "Gene_Symbol", "pChEMBL", "Evidence"])
    if len(df_out) > 0:
        df_out.to_csv(cache_file, index=False)
    else:
        print(f"[WARN] Skipped overwriting {cache_file} because query returned 0 rows.", flush=True)
    return df_out


def run_string_live(seed_genes, out_dir, prefix):
    seed_genes = sorted(set(seed_genes))
    print(f"\n[{prefix}] Querying Live STRING v12.0 API on {len(seed_genes)} genes...", flush=True)
    net_res = requests.post(
        "https://string-db.org/api/json/network",
        data={"identifiers": "%0d".join(seed_genes), "species": 9606, "required_score": 400, "caller_identity": "gcuf_pharmd_research"},
        timeout=30
    ).json()
    edges = []
    seen_e = set()
    seed_set = set(seed_genes)
    for it in net_res:
        a, b, sc = it.get("preferredName_A"), it.get("preferredName_B"), float(it.get("score", 0.0))
        if a in seed_set and b in seed_set and a != b:
            pair = tuple(sorted([a, b]))
            if pair not in seen_e:
                seen_e.add(pair)
                edges.append({"Protein_A": pair[0], "Protein_B": pair[1], "STRING_Score": round(sc, 3)})
    df_e = pd.DataFrame(edges, columns=["Protein_A", "Protein_B", "STRING_Score"])
    df_e.to_csv(os.path.join(out_dir, f"{prefix}_string_live_edges.csv"), index=False)

    G = nx.Graph()
    for n in seed_genes:
        G.add_node(n)
    for _, er in df_e.iterrows():
        G.add_edge(er["Protein_A"], er["Protein_B"], weight=float(er["STRING_Score"]))

    deg = dict(G.degree())
    bet = nx.betweenness_centrality(G, normalized=True)
    clo = nx.closeness_centrality(G)
    eig = nx.eigenvector_centrality(G, max_iter=1000) if G.number_of_edges() > 0 else {n: 0.0 for n in G.nodes()}
    c_rows = [{
        "Gene_Symbol": n, "Degree": deg[n], "Betweenness": round(bet[n], 4),
        "Closeness": round(clo[n], 4), "Eigenvector": round(eig[n], 4)
    } for n in G.nodes()]
    df_c = pd.DataFrame(c_rows).sort_values(["Degree", "Betweenness"], ascending=[False, False])
    df_c.to_csv(os.path.join(out_dir, f"{prefix}_hub_centrality_live.csv"), index=False)

    enr_res = requests.post(
        "https://string-db.org/api/json/enrichment",
        data={"identifiers": "%0d".join(seed_genes), "species": 9606, "caller_identity": "gcuf_pharmd_research"},
        timeout=30
    ).json()
    kegg_rows = [
        {
            "Term_ID": x.get("term"), "Pathway_Description": x.get("description"),
            "Category": x.get("category"), "Gene_Count": x.get("number_of_genes"),
            "Background_Genes": x.get("number_of_genes_in_background"),
            "p_value": x.get("p_value"), "FDR": x.get("fdr"),
            "Genes": ",".join(x.get("inputGenes", []))
        }
        for x in enr_res if x.get("category") in ("KEGG", "Process")
    ]
    df_k = pd.DataFrame(kegg_rows).sort_values("FDR") if kegg_rows else pd.DataFrame(
        columns=["Term_ID", "Pathway_Description", "Category", "Gene_Count", "Background_Genes", "p_value", "FDR", "Genes"]
    )
    df_k.to_csv(os.path.join(out_dir, f"{prefix}_string_live_enrichment.csv"), index=False)

    connected = [n for n, d in deg.items() if d >= 1]
    singletons = sorted([n for n, d in deg.items() if d == 0])
    print(f"  -> Nodes: {len(seed_genes)} total ({len(connected)} connected, {len(singletons)} singletons) | Edges: {len(df_e)}")
    if not df_k.empty:
        k_only = df_k[df_k["Category"] == "KEGG"].head(8)
        print(k_only[["Term_ID", "Pathway_Description", "Gene_Count", "FDR"]].to_string(index=False))
    return {
        "total_nodes": len(seed_genes),
        "connected_nodes": len(connected),
        "singleton_nodes": len(singletons),
        "singleton_list": singletons,
        "edges": len(df_e),
        "kegg_fdr05_count": int(sum((df_k["Category"] == "KEGG") & (df_k["FDR"] < 0.05))) if not df_k.empty else 0
    }


def main():
    force = "--refresh" in sys.argv
    ot_df = pd.read_csv(os.path.join(NET_A, "opentargets_disease_genes_live.csv"))
    ot_set = set(ot_df["Gene_Symbol"])

    # Protocol A Compounds
    df_a = pd.read_csv(os.path.join(BASE_DIR, "results", "admet_and_physicochemical_properties.csv"))
    print("=== QUERYING STRICT SINGLE-PROTEIN nM CHEMBL TARGETS (PROTOCOL A) ===", flush=True)
    chem_a = query_chembl_targets_for_compounds(df_a, os.path.join(NET_A, "chembl_experimental_targets_live.csv"), force_refresh=force)
    a_genes = set(chem_a["Gene_Symbol"])
    overlap_a = sorted(a_genes.intersection(ot_set))

    emp_stats_a = run_string_live(overlap_a, NET_A, "protocol_a_empirical")
    anchors_a = ["CD36", "FFAR4", "FFAR1", "GLP1R", "GCG", "DPP4", "PPARG", "PPARA", "PRKAA1", "AKT1"]
    bridge_stats_a = run_string_live(sorted(set(overlap_a + anchors_a)), NET_A, "protocol_a")

    ffar4_hits_a = sorted(chem_a[chem_a["Gene_Symbol"] == "FFAR4"]["Compound"].unique().tolist())
    cd36_hits_a = sorted(chem_a[chem_a["Gene_Symbol"] == "CD36"]["Compound"].unique().tolist())
    gap_report_a = {
        "Total_Unique_ChEMBL_Targets": len(a_genes),
        "Total_OpenTargets_Phenotype_Genes": len(ot_set),
        "Empirical_ChEMBL_OpenTargets_Overlap": len(overlap_a),
        "Overlap_List": overlap_a,
        "CD36_in_Empirical_ChEMBL_Hits": "CD36" in a_genes,
        "CD36_Compound_Hits": cd36_hits_a,
        "FFAR4_in_Empirical_ChEMBL_Hits": "FFAR4" in a_genes,
        "FFAR4_Compound_Hits": ffar4_hits_a,
        "Empirical_STRING_Stats": emp_stats_a,
        "Bridge_Genes_Added": sorted(set(anchors_a) - set(overlap_a)),
        "Bridge_STRING_Stats": bridge_stats_a
    }
    with open(os.path.join(NET_A, "empirical_target_gap_report.json"), "w") as f:
        json.dump(gap_report_a, f, indent=2)
    print("\n[Protocol A Gap Report]", json.dumps(gap_report_a, indent=2))

    # Protocol B Compounds
    pb_csv = os.path.join(BASE_DIR, "protocol_b_hepatic_cd36_dpp4", "results", "protocol_b_cd36_dpp4_complete_results.csv")
    if os.path.exists(pb_csv):
        df_b = pd.read_csv(pb_csv)
        print("\n=== QUERYING STRICT SINGLE-PROTEIN nM CHEMBL TARGETS (PROTOCOL B) ===", flush=True)
        chem_b = query_chembl_targets_for_compounds(df_b, os.path.join(NET_B, "chembl_experimental_targets_live.csv"), force_refresh=force)
        b_genes = set(chem_b["Gene_Symbol"])
        overlap_b = sorted(b_genes.intersection(ot_set))

        emp_stats_b = run_string_live(overlap_b, NET_B, "protocol_b_empirical")
        anchors_b = ["CD36", "DPP4", "GLP1R", "GCG", "ACSL4", "GPX4"]
        bridge_stats_b = run_string_live(sorted(set(overlap_b + anchors_b)), NET_B, "protocol_b")

        dpp4_hits_b = sorted(chem_b[chem_b["Gene_Symbol"] == "DPP4"]["Compound"].unique().tolist())
        cd36_hits_b = sorted(chem_b[chem_b["Gene_Symbol"] == "CD36"]["Compound"].unique().tolist())
        gap_report_b = {
            "Total_Unique_ChEMBL_Targets": len(b_genes),
            "Total_OpenTargets_Phenotype_Genes": len(ot_set),
            "Empirical_ChEMBL_OpenTargets_Overlap": len(overlap_b),
            "Overlap_List": overlap_b,
            "CD36_in_Empirical_ChEMBL_Hits": "CD36" in b_genes,
            "CD36_Compound_Hits": cd36_hits_b,
            "DPP4_in_Empirical_ChEMBL_Hits": "DPP4" in b_genes,
            "DPP4_Compound_Hits": dpp4_hits_b,
            "Empirical_STRING_Stats": emp_stats_b,
            "Bridge_Genes_Added": sorted(set(anchors_b) - set(overlap_b)),
            "Bridge_STRING_Stats": bridge_stats_b
        }
        with open(os.path.join(NET_B, "empirical_target_gap_report.json"), "w") as f:
            json.dump(gap_report_b, f, indent=2)
        print("\n[Protocol B Gap Report]", json.dumps(gap_report_b, indent=2))


if __name__ == "__main__":
    main()
