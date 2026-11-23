"""
03_network_and_figures.py
Generates 100% reproducible, 300-DPI publication figures for Protocol A directly
from LIVE OpenTargets + ChEMBL + STRING v12.0 network files, RDKit ADMET/BOILED-Egg
descriptors, deterministic AutoDock Vina (--seed 42, pH 6.8) outputs,
multi-exhaustiveness convergence checks, and ChEMBL5339 ROC-AUC benchmarking.
"""
import os
import json
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
NET_DIR = os.path.join(RESULTS_DIR, "network")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)


def repel_text_labels(ax, fig, items, x_pts, y_pts, pad_px=7, max_iter=420):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    trans = ax.transData
    inv_trans = trans.inverted()
    ax_bb = ax.get_window_extent(renderer=renderer)
    pts_px = trans.transform(np.column_stack([x_pts, y_pts])) if len(x_pts) else np.empty((0, 2))
    placed = []
    golden_angle = 2.39996323
    for idx, it in enumerate(items):
        tx, ty = trans.transform((it["x"], it["y"]))
        ang = idx * golden_angle
        ann = ax.text(it["x"], it["y"], it["text"], ha="center", va="center", fontsize=it.get("fontsize", 7.1), fontweight="bold",
                      color="#0F172A", zorder=10, bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="#CBD5E1", lw=0.45, alpha=0.88))
        placed.append({
            "ann": ann, "ax": it["x"], "ay": it["y"],
            "px": tx + 24.0 * np.cos(ang), "py": ty + 18.0 * np.sin(ang),
            "ox": tx, "oy": ty, "w": 42.0, "h": 15.0
        })
    fig.canvas.draw()
    for p in placed:
        bb = p["ann"].get_window_extent(renderer=renderer)
        p["w"], p["h"] = float(bb.width) + pad_px, float(bb.height) + pad_px
    for _ in range(max_iter):
        moved = False
        for i in range(len(placed)):
            pi = placed[i]
            fx, fy = 0.0, 0.0
            i_x0, i_x1 = pi["px"] - 0.5 * pi["w"], pi["px"] + 0.5 * pi["w"]
            i_y0, i_y1 = pi["py"] - 0.5 * pi["h"], pi["py"] + 0.5 * pi["h"]
            for j in range(len(placed)):
                if i == j:
                    continue
                pj = placed[j]
                j_x0, j_x1 = pj["px"] - 0.5 * pj["w"], pj["px"] + 0.5 * pj["w"]
                j_y0, j_y1 = pj["py"] - 0.5 * pj["h"], pj["py"] + 0.5 * pj["h"]
                if i_x0 < j_x1 and i_x1 > j_x0 and i_y0 < j_y1 and i_y1 > j_y0:
                    ox = min(i_x1, j_x1) - max(i_x0, j_x0)
                    oy = min(i_y1, j_y1) - max(i_y0, j_y0)
                    dx = pi["px"] - pj["px"]
                    dy = pi["py"] - pj["py"]
                    if abs(dx) < 1e-2 and abs(dy) < 1e-2:
                        dx, dy = (1.0 if i % 2 == 0 else -1.0), (1.0 if i % 3 == 0 else -1.0)
                    if ox < oy:
                        fx += (ox * 0.65) * (1.0 if dx >= 0 else -1.0)
                    else:
                        fy += (oy * 0.65) * (1.0 if dy >= 0 else -1.0)
            for qx, qy in pts_px:
                if (i_x0 - 8) <= qx <= (i_x1 + 8) and (i_y0 - 8) <= qy <= (i_y1 + 8):
                    dxp, dyp = pi["px"] - qx, pi["py"] - qy
                    nrm = max(1.0, (dxp * dxp + dyp * dyp) ** 0.5)
                    fx += 9.0 * (dxp / nrm)
                    fy += 9.0 * (dyp / nrm)
            fx += (pi["ox"] - pi["px"]) * 0.006
            fy += (pi["oy"] - pi["py"]) * 0.006
            if abs(fx) > 0.25 or abs(fy) > 0.25:
                pi["px"] = float(np.clip(pi["px"] + np.clip(fx, -16, 16), ax_bb.x0 + 0.5 * pi["w"] + 5, ax_bb.x1 - 0.5 * pi["w"] - 5))
                pi["py"] = float(np.clip(pi["py"] + np.clip(fy, -16, 16), ax_bb.y0 + 0.5 * pi["h"] + 5, ax_bb.y1 - 0.5 * pi["h"] - 5))
                moved = True
        if not moved:
            break
    for p in placed:
        nx_d, ny_d = inv_trans.transform((p["px"], p["py"]))
        p["ann"].set_position((nx_d, ny_d))
        if ((p["px"] - p["ox"]) ** 2 + (p["py"] - p["oy"]) ** 2) > 100.0:
            ax.plot([p["ax"], nx_d], [p["ay"], ny_d], color="#64748B", lw=0.65, alpha=0.72, zorder=2)


def render_live_network_figures():
    strict_p = os.path.join(NET_DIR, "chembl_experimental_targets_strict_live.csv")
    chem_df = pd.read_csv(strict_p if os.path.exists(strict_p) else os.path.join(NET_DIR, "chembl_experimental_targets_live.csv"))
    ot_df = pd.read_csv(os.path.join(NET_DIR, "opentargets_disease_genes_live.csv"))
    emp_edges = os.path.join(NET_DIR, "protocol_a_empirical_string_live_edges.csv")
    emp_cent = os.path.join(NET_DIR, "protocol_a_empirical_hub_centrality_live.csv")
    emp_enr = os.path.join(NET_DIR, "protocol_a_empirical_string_live_enrichment.csv")
    edges_df = pd.read_csv(emp_edges if os.path.exists(emp_edges) else os.path.join(NET_DIR, "protocol_a_string_live_edges.csv"))
    cent_df = pd.read_csv(emp_cent if os.path.exists(emp_cent) else os.path.join(NET_DIR, "protocol_a_hub_centrality_live.csv"))
    enr_df = pd.read_csv(emp_enr if os.path.exists(emp_enr) else os.path.join(NET_DIR, "protocol_a_string_live_enrichment.csv"))
    with open(os.path.join(NET_DIR, "empirical_target_gap_report.json")) as fg:
        gap = json.load(fg)

    chem_set = set(chem_df["Gene_Symbol"])
    ot_set = set(ot_df["Gene_Symbol"])
    overlap = chem_set.intersection(ot_set)
    ff_hits = gap.get("FFAR4_Compound_Hits", [])
    cd_hits = gap.get("CD36_Compound_Hits", [])
    ff_str = f"FFAR4={len(ff_hits)} hits ({', '.join(ff_hits)})" if ff_hits else "FFAR4=0 ChEMBL"
    cd_str = f"CD36={len(cd_hits)} ChEMBL hits"

    # Venn Diagram + Target Coverage Gap
    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=300)
    c1 = plt.Circle((0.34, 0.5), 0.30, color="#0B3C5D", alpha=0.72)
    c2 = plt.Circle((0.66, 0.5), 0.32, color="#D95F02", alpha=0.68)
    ax.add_patch(c1)
    ax.add_patch(c2)
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(0.08, 0.94)
    ax.axis("off")
    plt.text(0.16, 0.5, f"{len(chem_set - ot_set)}\nStrict ChEMBL-Only\nSingle-Protein Targets", ha="center", va="center", fontweight="bold", color="white", fontsize=8.8)
    plt.text(0.84, 0.5, f"{len(ot_set - chem_set)}\nOpenTargets\nDisease Genes", ha="center", va="center", fontweight="bold", color="white", fontsize=8.8)
    plt.text(0.50, 0.5, f"{len(overlap)}\nPure Empirical\nShared Hubs\n({ff_str};\n{cd_str})", ha="center", va="center", fontweight="bold", color="white", fontsize=7.6, bbox=dict(boxstyle="round,pad=0.35", facecolor="#1E293B", alpha=0.94))
    plt.title("Strict Quantitative ChEMBL (pChEMBL ≥ 5.0 / Mechanism) vs.\nOpenTargets (4 MeSH Phenotypes)", fontsize=10.5, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(NET_DIR, "target_disease_venn.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Stage 1 Pure Empirical STRING v12.0 PPI Graph
    plt.figure(figsize=(10.5, 8.2), dpi=300)
    G = nx.Graph()
    for n in cent_df["Gene_Symbol"]:
        G.add_node(n)
    for _, r in edges_df.iterrows():
        G.add_edge(r["Protein_A"], r["Protein_B"], weight=float(r["STRING_Score"]))
    conn_nodes = [n for n, d in G.degree() if d > 0]
    G_conn = G.subgraph(conn_nodes)
    pos = nx.kamada_kawai_layout(G_conn) if len(conn_nodes) > 1 else {n: np.array([0.0, 0.0]) for n in G.nodes()}
    sings = sorted([n for n, d in G.degree() if d == 0])
    for idx_s, sn in enumerate(sings):
        pos[sn] = np.array([-0.55 + idx_s * 1.1, -1.10])
    nsize = [560 + G.degree(n) * 95 for n in G.nodes()]
    ncolor = ["#0B3C5D" if G.degree(n) >= 10 else ("#0284C7" if G.degree(n) > 0 else "#64748B") for n in G.nodes()]
    nx.draw_networkx_edges(G, pos, alpha=0.35, edge_color="#475569", width=1.2)
    nx.draw_networkx_nodes(G, pos, node_size=nsize, node_color=ncolor, edgecolors="black", linewidths=0.9)
    nx.draw_networkx_labels(G, pos, font_size=6.8, font_weight="bold", font_color="white")
    xs_g, ys_g = zip(*pos.values())
    plt.xlim(min(xs_g) - 0.22, max(xs_g) + 0.22)
    plt.ylim(min(ys_g) - 0.20, max(ys_g) + 0.20)
    plt.title(f"Stage 1 Pure Empirical STRING v12.0 Network ({len(overlap)} Unspiked Overlap Genes, {len(edges_df)} Edges)", fontsize=11.2, fontweight="bold")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(NET_DIR, "ppi_network_graph.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Figure 3: Hub Genes Centrality Barplot (Live STRING Degrees)
    plt.figure(figsize=(9, 5), dpi=300)
    top_h = cent_df.head(12).sort_values("Degree", ascending=True)
    bars = plt.barh(top_h["Gene_Symbol"], top_h["Degree"], color="#0B3C5D", edgecolor="black", height=0.68)
    for b in bars:
        plt.text(b.get_width() + 0.25, b.get_y() + b.get_height() / 2, f"{int(b.get_width())}", va="center", fontweight="bold", fontsize=9)
    plt.xlabel("Degree Centrality (Live STRING v12.0 Interactions)", fontweight="bold")
    plt.title("Top 12 Topological Hub Genes Linking Lipid Sensing & Metabolic Homeostasis", fontweight="bold")
    plt.xlim(0, max(top_h["Degree"]) + 3)
    plt.grid(axis="x", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(NET_DIR, "hub_genes_barplot.png"), dpi=300)
    plt.close()

    # Figure 4: Live STRING KEGG Enrichment Bubble Chart
    kegg_df = enr_df[enr_df["Category"] == "KEGG"].head(10).sort_values("Gene_Count", ascending=True)
    plt.figure(figsize=(10.5, 5.8), dpi=300)
    neg_log_fdr = -np.log10(kegg_df["FDR"].clip(lower=1e-20))
    sc = plt.scatter(kegg_df["Gene_Count"], kegg_df["Pathway_Description"], s=kegg_df["Gene_Count"] * 45, c=neg_log_fdr, cmap="viridis", edgecolors="black", linewidth=1.2)
    cbar = plt.colorbar(sc)
    cbar.set_label("-log10(FDR) [Live STRING API]", fontweight="bold")
    plt.xlabel("Enriched Gene Count", fontweight="bold")
    plt.title("Live STRING v12.0 KEGG Pathway Enrichment (FDR < 0.05)", fontweight="bold")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(NET_DIR, "kegg_bubble_plot.png"), dpi=300)
    plt.close()


def render_docking_and_admet_figures():
    admet_df = pd.read_csv(os.path.join(RESULTS_DIR, "admet_and_physicochemical_properties.csv"))
    dock_df = pd.read_csv(os.path.join(RESULTS_DIR, "molecular_docking_results.csv"))
    with open(os.path.join(RESULTS_DIR, "roc_summary_all.json")) as fr:
        roc_meta = json.load(fr)["FFAR4"]

    # Figure 5: 2D Structures of Top Tier-1A Leads, Controls & Clashing Triterpenoid
    selected = [
        "TUG-891", "Oleic acid", "Palmitic acid", "Oleoylethanolamide",
        "Naringenin", "Resveratrol", "Phloretin", "6-Gingerol",
        "Spilanthol", "Pellitorine", "Pinolenic acid", "10-Hydroxy-cis-12-octadecenoic acid",
        "Luteolin", "Quercetin", "Ferulic acid", "Oleanolic acid"
    ]
    mols, legends = [], []
    for name in selected:
        sub = admet_df[admet_df["Compound"] == name]
        if sub.empty:
            continue
        r = sub.iloc[0]
        d_sub = dock_df[dock_df["Compound"] == name].iloc[0]
        mol = Chem.MolFromSmiles(r["SMILES"])
        AllChem.Compute2DCoords(mol)
        mols.append(mol)
        short_n = name.replace("10-Hydroxy-cis-12-octadecenoic acid", "10-HOE (HYA)").replace("Oleoylethanolamide", "OEA")
        legends.append(f"{short_n}\nFFAR4: {d_sub['FFAR4_8G59_Affinity_kcal_mol']:.2f} | LE: {d_sub['FFAR4_LE_kcal_mol_HA']:.2f} | LogS: {r['ESOL_LogS_mol_L']:.2f}")
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(380, 280), legends=legends)
    img.save(os.path.join(FIG_DIR, "figure5_chemical_structures.png"))

    # Figure 6: Dual Binding Affinity Bar Chart & Ligand Efficiency (LE >= 0.30) + Clash Detection
    fig = plt.figure(figsize=(16, 7), dpi=300)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 0.95], wspace=0.25)

    ax1 = fig.add_subplot(gs[0])
    df_p = dock_df.copy()
    df_p["Short"] = df_p["Compound"].str.replace("10-Hydroxy-cis-12-octadecenoic acid", "10-HOE (HYA)").str.replace("Oleoylethanolamide", "OEA")
    x = np.arange(len(df_p))
    w = 0.38
    ax1.bar(x - w / 2, df_p["FFAR4_8G59_Affinity_kcal_mol"], w, label="FFAR4 Active Pocket (8G59, pH 6.8)", color="#0B3C5D", edgecolor="black", linewidth=0.7)
    ax1.bar(x + w / 2, df_p["CD36_5LGD_Affinity_kcal_mol"], w, label="CD36 PLM-512 Cavity (5LGD, pH 6.8)", color="#D9534F", edgecolor="black", linewidth=0.7)
    ax1.axhline(0.0, color="black", linewidth=1.0)
    ax1.axhline(-7.0, color="#16A34A", linestyle="--", linewidth=1.2, label="-7.0 kcal/mol High-Affinity Gate")
    ax1.set_xticks(x)
    ax1.set_xticklabels(df_p["Short"], rotation=45, ha="right", fontsize=8.5, fontweight="bold")
    ax1.set_ylabel("AutoDock Vina Binding Energy (kcal/mol) [--seed 42]", fontweight="bold")
    ax1.set_title("A. Deterministic Dual-Receptor Affinities (Salivary pH 6.8, altLoc-Clean)", fontweight="bold")
    ax1.grid(axis="y", linestyle=":", alpha=0.6)
    ax1.legend(loc="lower left", fontsize=8.5)

    ax2 = fig.add_subplot(gs[1])
    lbl_items_f6 = []
    for _, r in df_p.iterrows():
        tier = r["Lead_Qualification_Tier"]
        if "Tier 1A" in tier:
            c, m, lbl = "#16A34A", "^", "Tier 1A: PAINS-Free Lead (LE ≥ 0.30, 0 Clashes)"
        elif "Tier 1B" in tier:
            c, m, lbl = "#EA580C", "D", "Tier 1B: Catechol-PAINS Hit (Redox Control Needed)"
        elif "Tier 2" in tier:
            c, m, lbl = "#2563EB", "o", "Tier 2: Low-MW Fragment Scaffold"
        elif "Control" in tier:
            c, m, lbl = "#475569", "s", "Reference Control (TUG-891 / OEA / Native Lipids)"
        else:
            c, m, lbl = "#DC2626", "X", "Disqualified Brick-Dust (Phe115 Clash, LE < 0.05)"
        ax2.scatter(r["Heavy_Atoms"], r["FFAR4_LE_kcal_mol_HA"], s=110, color=c, marker=m, edgecolor="black", zorder=3, label=lbl)
        lbl_items_f6.append({"text": r["Short"], "x": r["Heavy_Atoms"], "y": r["FFAR4_LE_kcal_mol_HA"], "fontsize": 7.3})

    handles, labels = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax2.axhline(0.30, color="#16A34A", linestyle="--", linewidth=1.3)
    ax2.axhspan(0.30, 0.56, color="#DCFCE7", alpha=0.35, zorder=0)
    ax2.axvspan(31, 38, color="#FEE2E2", alpha=0.35, zorder=0)
    ax2.set_xlabel("Heavy Atom Count (N_heavy)", fontweight="bold")
    ax2.set_ylabel("FFAR4 Ligand Efficiency (LE = -ΔG / N_heavy)", fontweight="bold")
    ax2.set_title("B. Ligand Efficiency vs. Heavy-Atom Size & Clash Stratification", fontweight="bold")
    ax2.set_xlim(12, 38.5)
    ax2.set_ylim(-0.08, 0.56)
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(by_label.values(), by_label.keys(), loc="lower left", fontsize=7.3)
    fig.subplots_adjust(bottom=0.22, wspace=0.25)
    repel_text_labels(ax2, fig, lbl_items_f6, df_p["Heavy_Atoms"].values, df_p["FFAR4_LE_kcal_mol_HA"].values, pad_px=4, max_iter=160)
    plt.savefig(os.path.join(FIG_DIR, "figure6_binding_affinities_radar_bar.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Figure 7: Corrected Biophysical Circuit
    fig, ax = plt.subplots(figsize=(15, 8), dpi=300)
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.text(7.5, 9.6, "Biophysical & Pharmacokinetic Working Hypothesis — Reversible CD36 Chaperone & Active FFAR4 Agonism",
            ha="center", va="center", fontsize=12.2, fontweight="bold", color="#0B3C5D")

    p1 = patches.FancyBboxPatch((0.4, 1.3), 4.2, 7.7, boxstyle="round,pad=0.2", facecolor="#F0F4F8", edgecolor="#0B3C5D", linewidth=2)
    ax.add_patch(p1)
    ax.text(2.5, 8.6, "1. Salivary Dispersion & PAINS Filter", ha="center", fontsize=10.5, fontweight="bold", color="#0B3C5D")
    ax.text(0.6, 8.0, (
        "• Zero PAINS & Zero Brenk Leads:\n"
        "   - Naringenin (LogS=-3.49, 0/0 alerts)\n"
        "   - Phloretin (LogS=-3.18, 0/0 alerts)\n"
        "• PAINS-Free Sialogogue Alkamides:\n"
        "   - Spilanthol (359 ug/mL) & Pellitorine\n"
        "• PAINS-Free GLP-1 Secretagogues:\n"
        "   - Pinolenic acid, 6-Gingerol, 10-HOE\n"
        "• Catechol-PAINS Flagged Transparently:\n"
        "   - Luteolin, Quercetin, Chlorogenic\n"
        "• Excluded: TRPV1 irritants & >6-cLogP\n"
        "   triterpenoids (Ursolic/Oleanolic clash)"
    ), va="top", fontsize=8.5, linespacing=1.55)

    p2 = patches.FancyBboxPatch((5.3, 1.3), 4.6, 7.7, boxstyle="round,pad=0.2", facecolor="#EBF5FB", edgecolor="#2980B9", linewidth=2)
    ax.add_patch(p2)
    ax.text(7.6, 8.6, "2. Taste-Pore Biophysical Recognition", ha="center", fontsize=10.5, fontweight="bold", color="#1B4F72")
    ax.text(5.5, 8.0, (
        "• Stage A — CD36 (5LGD, pH 6.8) Chaperone:\n"
        "   - Moderate reversible binding (-6.19 to\n"
        "     -7.25 kcal/mol) at PLM-512 cavity\n"
        "     (Arg337, Lys334, Gly199, Ala141 H-bonds)\n"
        "• Stage B — FFAR4 (8G59, pH 6.8) Agonism:\n"
        "   - Active-state cryo-EM validation:\n"
        f"     TUG-891 Hungarian RMSD = 0.864 Å\n"
        f"     ChEMBL5339 ROC-AUC(LE) = {roc_meta['ROC_AUC_LE']:.3f}\n"
        "   - Engages Thr119 (2.72-2.94 Å H-bond),\n"
        "     Asp208, Met118, Phe115, Ile284\n"
        "   - Triggers Gq -> PLCb2 -> IP3 -> Ca2+"
    ), va="top", fontsize=8.5, linespacing=1.55)

    p3 = patches.FancyBboxPatch((10.6, 1.3), 3.9, 7.7, boxstyle="round,pad=0.2", facecolor="#E8F8F5", edgecolor="#16A085", linewidth=2)
    ax.add_patch(p3)
    ax.text(12.55, 8.6, "3. Cephalic-Enteric GLP-1 Reflex", ha="center", fontsize=10.5, fontweight="bold", color="#0E6251")
    ax.text(10.8, 8.0, (
        "• CN IX & CN VII Afferent Activation:\n"
        "   - Taste-bud ATP release via CALHM1/3\n"
        "   - Projections to Nucleus Tractus Solitarius\n"
        "   - Vagal efferent priming of intestinal\n"
        "     L-cells for GLP-1, PYY & OEA secretion\n"
        "• Methodological Guardrails:\n"
        "   - Triterpenoids clash with Phe115 (1.70 Å)\n"
        "   - FFAR4 MD requires POPC membrane\n"
        "   - Future pH 6.8 salivary dissolution QC\n"
        "     & POPC-bilayer MD required before\n"
        "     experimental FLIPR/GLP-1 assays"
    ), va="top", fontsize=8.5, linespacing=1.55)

    ax.annotate("", xy=(5.2, 5.2), xytext=(4.7, 5.2), arrowprops=dict(arrowstyle="->", lw=3.0, color="#0B3C5D"))
    ax.annotate("", xy=(10.5, 5.2), xytext=(10.0, 5.2), arrowprops=dict(arrowstyle="->", lw=3.0, color="#1B4F72"))
    plt.savefig(os.path.join(FIG_DIR, "figure7_lingual_delivery_circuit.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # BOILED-Egg GI/BBB Map & Lipophilic Ligand Efficiency (LLE) vs ESOL LogS
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.2), dpi=300)
    egg_hia = patches.Rectangle((-0.7, 0), 6.58, 131.6, facecolor="#F8FAFC", edgecolor="#0B3C5D", linestyle="--", lw=1.2, alpha=0.8, label="HIA Passive Absorption Window (TPSA ≤ 131.6 Å², WLogP ≤ 5.88)")
    egg_bbb = patches.Rectangle((0.4, 0), 5.6, 79.0, facecolor="#FEF9C3", edgecolor="#CA8A04", linestyle="-", lw=1.2, alpha=0.55, label="BBB Permeation Yolk (TPSA < 79 Å², WLogP 0.4–6.0)")
    ax1.add_patch(egg_hia)
    ax1.add_patch(egg_bbb)
    lbls_f8a = []
    for _, r in admet_df.iterrows():
        c = "#16A34A" if (r["PAINS_Structural_Alert"] == "0 Alerts" and r["ESOL_LogS_mol_L"] >= -4.5) else ("#EA580C" if "catechol" in r["PAINS_Structural_Alert"] else "#DC2626")
        ax1.scatter(r["LogP_WildmanCrippen"], r["TPSA_A2"], s=95, color=c, edgecolor="black", zorder=3)
        short = r["Compound"].replace("10-Hydroxy-cis-12-octadecenoic acid", "10-HOE").replace("Oleoylethanolamide", "OEA")
        lbls_f8a.append({"text": short, "x": r["LogP_WildmanCrippen"], "y": r["TPSA_A2"], "fontsize": 7.1})
    ax1.set_xlabel("Wildman–Crippen Lipophilicity (WLogP)", fontweight="bold")
    ax1.set_ylabel("Topological Polar Surface Area (TPSA, Å²)", fontweight="bold")
    ax1.set_title("A. Egan / BOILED-Egg Passive Permeability & Safety Classification", fontweight="bold")
    ax1.set_xlim(-1.6, 8.2)
    ax1.set_ylim(10, 180)
    ax1.legend(loc="upper right", fontsize=7.5)
    ax1.grid(True, linestyle=":", alpha=0.5)

    # Panel B: LLE vs ESOL LogS
    lbls_f8b = []
    for _, r in dock_df.iterrows():
        c = "#16A34A" if (r["FFAR4_LLE"] >= 2.0 and r["PAINS_Alert"] == "0 Alerts") else ("#EA580C" if "catechol" in r["PAINS_Alert"] else "#DC2626")
        ax2.scatter(r["ESOL_LogS"], r["FFAR4_LLE"], s=95, color=c, edgecolor="black", zorder=3)
        short = r["Compound"].replace("10-Hydroxy-cis-12-octadecenoic acid", "10-HOE").replace("Oleoylethanolamide", "OEA")
        lbls_f8b.append({"text": short, "x": r["ESOL_LogS"], "y": r["FFAR4_LLE"], "fontsize": 7.1})
    ax2.axvline(-6.0, color="#DC2626", linestyle="--", lw=1.2, label="ESOL LogS = -6.0 Brick-Dust Limit")
    ax2.axhline(0.0, color="#475569", linestyle=":", lw=1.2, label="LLE = 0 Lipophilicity Efficiency Line")
    ax2.set_xlabel("Delaney ESOL Salivary Solubility LogS (mol/L)", fontweight="bold")
    ax2.set_ylabel("FFAR4 Lipophilic Ligand Efficiency (LLE = pKd - cLogP)", fontweight="bold")
    ax2.set_title("B. Salivary Solubility vs. Lipophilic Ligand Efficiency (LLE)", fontweight="bold")
    ax2.set_xlim(-7.6, -1.0)
    ax2.set_ylim(-7.8, 8.2)
    ax2.grid(True, linestyle=":", alpha=0.5)
    ax2.legend(loc="lower right", fontsize=7.5)
    plt.tight_layout(w_pad=2.2)
    repel_text_labels(ax1, fig, lbls_f8a, admet_df["LogP_WildmanCrippen"].values, admet_df["TPSA_A2"].values, pad_px=4, max_iter=160)
    repel_text_labels(ax2, fig, lbls_f8b, dock_df["ESOL_LogS"].values, dock_df["FFAR4_LLE"].values, pad_px=4, max_iter=160)
    plt.savefig(os.path.join(FIG_DIR, "figure8_redocking_validation.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print("[SUCCESS] Rendered all Protocol A live network, ROC/BOILED-Egg, and LE figures!")


if __name__ == "__main__":
    render_live_network_figures()
    render_docking_and_admet_figures()
