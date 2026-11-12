import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

BASE_DIR = r"C:\Users\Abdullah\.gemini\antigravity\scratch\cd36_ffar4_glp1_study"
RESULTS_DIR = os.path.join(BASE_DIR, "results", "network")

edges_path = os.path.join(RESULTS_DIR, "ppi_network_edges.csv")
centrality_path = os.path.join(RESULTS_DIR, "hub_genes_centrality.csv")

df_edges = pd.read_csv(edges_path)
df_centrality = pd.read_csv(centrality_path)

G = nx.Graph()
for _, r in df_edges.iterrows():
    G.add_edge(r["Protein_A"], r["Protein_B"], weight=float(r["STRING_Score"]))

deg_dict = dict(df_centrality.set_index("Gene_Symbol")["Degree"])
hub_dict = dict(df_centrality.set_index("Gene_Symbol")["Is_Hub"])

plt.figure(figsize=(14, 12), dpi=300)
pos = nx.spring_layout(G, k=0.45, seed=42, iterations=100)

# Colors and sizes
node_colors = []
node_sizes = []
for node in G.nodes():
    d = deg_dict.get(node, G.degree(node))
    is_hub = hub_dict.get(node, "No")
    if node in ["CD36", "FFAR4"]:
        node_colors.append("#e41a1c")  # Primary Receptors - Red
        node_sizes.append(1800)
    elif node in ["GLP1R", "DPP4", "GCG", "PYY", "CCK"]:
        node_colors.append("#377eb8")  # Incretin Axis - Blue
        node_sizes.append(1400)
    elif is_hub == "Yes":
        node_colors.append("#ff7f00")  # Major Hubs - Orange
        node_sizes.append(1200)
    else:
        node_colors.append("#4daf4a")  # Effector nodes - Green
        node_sizes.append(600)

# Draw edges with transparency
nx.draw_networkx_edges(G, pos, alpha=0.25, edge_color="gray", width=0.8)

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, edgecolors="black", linewidths=1.2)

# Draw labels
labels = {node: node for node in G.nodes()}
nx.draw_networkx_labels(G, pos, labels=labels, font_size=9, font_weight="bold", font_family="sans-serif")

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Primary Targets (CD36, FFAR4)', markerfacecolor='#e41a1c', markersize=14, markeredgecolor='black'),
    Line2D([0], [0], marker='o', color='w', label='Incretin & Satiety Axis (GLP1R, DPP4, GCG)', markerfacecolor='#377eb8', markersize=12, markeredgecolor='black'),
    Line2D([0], [0], marker='o', color='w', label='Metabolic Hub Regulators (PPARA, PPARG, AKT1, MTOR)', markerfacecolor='#ff7f00', markersize=11, markeredgecolor='black'),
    Line2D([0], [0], marker='o', color='w', label='Secondary Signaling Effectors', markerfacecolor='#4daf4a', markersize=8, markeredgecolor='black')
]
plt.legend(handles=legend_elements, loc="upper right", fontsize=10, frameon=True, framealpha=0.9)

plt.title("Protein-Protein Interaction (PPI) Network of CD36 + FFAR4 Post-GLP-1 Satiety Hubs\n(STRING v12.0 Confidence >= 0.400)", fontsize=13, fontweight="bold", pad=20)
plt.axis("off")
plt.tight_layout()

out_path = os.path.join(RESULTS_DIR, "ppi_network_plot.png")
plt.savefig(out_path, dpi=300)
plt.close()
print(f"[SUCCESS] Saved publication-grade PPI network plot to: {out_path}")
