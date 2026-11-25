# Hepatic CD36 and Vascular DPP-4 Structure-Based Drug Discovery

**Author:** Muhammad Abdullah Khan ([ORCID: 0009-0009-8867-6526](https://orcid.org/0009-0009-8867-6526))  
**Affiliation:** Department of Pharmaceutics, Faculty of Pharmaceutical Sciences, Government College University Faisalabad (GCUF), Faisalabad 38000, Pakistan  
**Contact:** abdullah.khan.pharmd@gmail.com  

## Overview

This repository archives the code, prepared structures, docked poses, continuum electrostatic logs, benchmarking tables, and explicit-water molecular dynamics outputs supporting the computational pharmacology manuscript:

> *Structure-Based Prioritization of Hepatic CD36 Lipid-Tunnel Modulators and Non-Covalent DPP-4 Inhibitors at Plasma pH 7.4: Unspiked Empirical Network Auditing, PROPKA Electrostatic Verification, Size-Matched Docking Validation, and 50 ns Explicit-Water Molecular Dynamics* (submitted to *Molecular Diversity*, Springer Nature).

The study evaluates 24 clinical gliptins, hepatoprotective phytochemicals, and native crystallographic controls against human soluble CD36 (PDB: 5LGD, 2.07 Å) and human DPP-4 (PDB: 4A5S, 1.62 Å) at physiological plasma pH 7.4. The workflow separates empirical single-protein ChEMBL/Open Targets/STRING v12.0 network enrichment from hypothesis-driven bridge graphs, verifies pocket micro-pKa states with continuum PROPKA 3.5, validates AutoDock Vina v1.2.7 via symmetry-corrected redocking, alternate-rotamer/multi-seed stability checks, and a size-matched N = 60 PubChem BioAssay benchmark, applies exact Daina–Zoete BOILED-Egg and transporter-aware pharmacokinetic gating, and tests dynamic pose retention across 400 ns (8 × 50 ns) of explicit-water OpenMM 8 molecular dynamics with end-point OBC2 MM-GBSA rescoring.

All results in this repository represent falsifiable dry-lab hypotheses awaiting experimental biochemical and cellular validation.

## Repository Structure

| Path | Description |
| :--- | :--- |
| `protocol_b_hepatic_cd36_dpp4/` | Primary hepatic CD36 (5LGD) and vascular DPP-4 (4A5S) plasma pH 7.4 docking pipeline, PROPKA 3.5 logs, N = 60 PubChem BioAssay validation, Val546 altLoc B and 3-seed audits, docked PDBQT poses, PyMOL sessions, and figure generation scripts. |
| `protocol_b_hepatic_cd36_dpp4/colab_md/` | Input bundles, OpenMM 8 / OpenFF Sage 2.2.0 / AmberTools simulation notebooks, RMSD/RMSF/H-bond/water-bridge logs, and OBC2 MM-GBSA summaries for all 8 × 50 ns trajectories. |
| `protocol_b_hepatic_cd36_dpp4/submission_molecular_diversity/` | Complete 8-sheet Supplementary Data Workbook (`04_Supplementary_Tables.xlsx`), Graphical Abstract, and 300-DPI LZW TIFF figures (`figures_300dpi_tiff/`). |
| `scripts/` | Shared structural preparation utilities (`structural_prep_utils.py`, strict alternate-conformation cleaning, Hungarian heavy-atom RMSD, exact BOILED-Egg geometry) and companion salivary pH 6.8 FFAR4 (8G59) / CD36 (5LGD) conformer-sensitivity scripts (`07_audit_8g59_phe115_and_fix_admet.py`). |
| `structures/`, `ligands/`, `results/` | Cleaned receptor structures, protonated ligand libraries, and backing CSV/JSON outputs for the shared structural utilities and companion FFAR4/CD36 sensitivity dataset. |

## Summary of Primary Validation & 400 ns MD Outputs

| Validation / Analysis | Key Quantitative Outcome |
| :--- | :--- |
| Empirical Target Audit (ChEMBL v33 ∩ Open Targets v4) | 125 single-protein ChEMBL targets (pChEMBL ≥ 5.0 or approved mechanism) ∩ 894 disease genes = 21 unspiked STRING v12.0 overlap proteins (44 edges; 4 DPP-4 gliptins, 0 strict CD36 hits). |
| Continuum PROPKA 3.5 (pH 7.4) | 0 protonation flips across 113 ionizable groups in 5LGD; 8/9 active-site anchors unchanged in 4A5S (buried catalytic His740 shifted to pKa 8.00). |
| Crystallographic Redocking & Conformer Stability | N7F-901 in 4A5S: −10.71 kcal/mol, 0.780 Å Hungarian heavy-atom RMSD; 0 clashes across all 24 poses against 4A5S Val546 altLoc B (min distance 2.720 Å). Palmitate (CID 985) in 5LGD: 0.156 Å perpendicular tunnel-axis error (3.098 Å heavy-atom RMSD across 14 torsions). |
| N = 60 Non-Covalent DPP-4 Benchmark (PubChem GeneID: 1803) | Raw Vina ROC-AUC = 0.686 (95% CI 0.546–0.813, p = 0.0068, PR-AUC = 0.711, EF20% = 1.83×); Bemis–Murcko scaffold-collapsed AUC = 0.665 (p = 0.0219); rotatable-bond confounder (mean RotB 3.20 vs 4.77) disclosed explicitly. |
| Prioritized Tier 1A Oral Hypotheses (0 PAINS, Veber Pass) | Linagliptin (−9.81 / −8.79 kcal/mol), Sitagliptin (−9.21 / −8.21), Teneligliptin (−8.99 / −7.37), and Naringenin (−8.19 / −7.25, ESOL 88.1 µg/mL); portal-hepatic Berberine separated via P-gp/OCT1 rules; Silybin A/B placed in phytosome-dependent Tier 2A. |
| 8 × 50 ns Explicit-Water OpenMM MD + OBC2 MM-GBSA | Compact protein backbones across all 8 complexes (Cα RMSD 1.25–1.94 Å). Tightest ligand retention in DPP-4: Naringenin (2.78 ± 0.30 Å, ΔG = −21.61 ± 0.17 kcal/mol). Tightest locking in CD36: Linagliptin (1.69 ± 0.25 Å, ΔG = −34.79 ± 0.27 kcal/mol, ARG63 78.8% HB / LYS334 89.6% WB), with GLU335 convergence across palmitate, teneligliptin, and naringenin. |

### Companion Exploratory Dataset (Lingual FFAR4 / CD36, Salivary pH 6.8)
The root `scripts/` and `results/` folders also archive an exploratory structural sensitivity dataset evaluating human FFAR4 (PDB: 8G59) and CD36 at salivary pH 6.8 (TUG-891 redocking RMSD = 0.864 Å), including empirical redocking across the 50/50 Phe115 alternate conformer (`results/audit_8g59_phe115_conformer_ab.csv`) demonstrating rotamer-dependent steric clashes for bulky triterpenoids versus conformer-stable compact scaffolds. This companion study remains in methodological revision and is separate from the submitted hepatic CD36/DPP-4 manuscript.

## Reproducibility

Install pinned Python dependencies from `requirements.txt`, along with OpenBabel (`obabel` on `PATH` or via `OBABEL_EXE`) and AutoDock Vina v1.2.7 (`vina` on `PATH` or `scripts/vina.exe`):

```bash
pip install -r requirements.txt

# 1. Re-run N=60 PubChem BioAssay bootstrap/scaffold & Val546 altLoc B / 3-seed audits
python protocol_b_hepatic_cd36_dpp4/03_hardening_audits_protocol_b.py

# 2. Regenerate primary hepatic CD36 / DPP-4 figures and summary tables from backing CSV/JSON outputs
python protocol_b_hepatic_cd36_dpp4/02_figures_and_manuscript_protocol_b.py

# 3. Regenerate Graphical Abstract
python protocol_b_hepatic_cd36_dpp4/submission_molecular_diversity/02_generate_graphical_abstract.py

# 4. Run companion 8G59 Phe115 Conformer A/B sensitivity & ADMET verification script
python scripts/07_audit_8g59_phe115_and_fix_admet.py
```

The 8 × 50 ns explicit-water OpenMM 8 MD simulations (AMBER ff14SB / TIP3P / OpenFF Sage 2.2.0 / AM1-BCC, 4 fs HMR, 310 K / 1 bar) can be reproduced using the input bundle and Colab GPU notebook in `protocol_b_hepatic_cd36_dpp4/colab_md/`.

## Software & Authorship

All analyses, scripts, structural preparations, and simulations were executed and verified by Muhammad Abdullah Khan using open-source scientific software (RDKit, OpenBabel, PROPKA 3.5, AutoDock Vina v1.2.7, PyMOL Open-Source, OpenMM 8, OpenFF Sage 2.2.0, AmberTools, MDTraj, SciPy, NetworkX, and Matplotlib).
