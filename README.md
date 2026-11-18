# Dual-Protocol Computational Drug Discovery Repository: Lingual FFAR4/CD36 (Protocol A) & Hepatic CD36/DPP-4 in MASLD (Protocol B)

**Author:** Muhammad Abdullah Khan  
**Affiliation:** Department of Pharmaceutics, Faculty of Pharmaceutical Sciences, Government College University Faisalabad (GCUF), Pakistan  
**Contact:** `abdullah.khan.pharmd@gmail.com`  
**Repository:** [https://github.com/abdullah-pharmd/CD36-FFAR4-GLP1-Discovery](https://github.com/abdullah-pharmd/CD36-FFAR4-GLP1-Discovery)

---

## 1. Scientific Scope & Architectural Rigor

This repository contains two complementary, reproducible in-silico pipelines designed to address post-GLP-1 discontinuation weight rebound and metabolic dysfunction-associated steatotic liver disease (MASLD):

1. **Protocol A (Lingual Fat-Taste Axis — Salivary pH 6.8):** Active-state human **FFAR4/GPR120** agonism (`PDB: 8G59`, Chain R, `2.64 Å`; `11` duplicate `altLoc B` atoms including orthosteric `Phe115` stripped; `424` polar `HD` donors) paired with reversible, palmitate-like **CD36** lipid-chaperone binding (`PDB: 5LGD`, Chain A, `2.07 Å`; `619` polar `HD` donors) across `22` ligands (`4` non-covalent controls, lingual sialogogue alkamides, GLP-1-secretagogue PUFAs, dietary polyphenols, and a `4`-triterpenoid size/clash negative-control panel).
2. **Protocol B (Hepatic MASLD / Ferroptosis & Plasma Incretin Axis — Plasma pH 7.4):** Dual steric/catalytic inhibition of the hepatic **CD36** `PLM-512` lipid-translocation tunnel (`PDB: 5LGD`, Chain A, `3,218` heavy atoms, `619` `HD` donors) and human **DPP-4** (`PDB: 4A5S`, Chain A, `1.62 Å`; `48` duplicate `altLoc B` atoms including `Val546` stripped, `5,960` heavy atoms, `1,087` `HD` donors) across `24` compounds (clinical non-covalent/covalent gliptins, biochemical CD36 antagonist/lipid controls, milk-thistle flavonolignans, caffeic acid esters, coumarins, flavonoids, and berberine).

All results satisfy strict methodological audit standards:
- **Zero Circular Target Dictionaries:** Live querying of **ChEMBL REST API** (`<= 10 uM` human single-protein assays), **Open Targets GraphQL v4** (`894` disease genes across `MONDO_0011122`, `MONDO_0005148`, `EFO_1001248`, `EFO_1000653`), and **STRING v12.0** (`required_score >= 400`).
- **True `altLoc` Chain Surgery & Physiological Ionization:** Alternate side-chain duplicates (`altLoc B`) removed prior to physiological protonation (`obabel -p 6.8` for Protocol A; `-p 7.4` for Protocol B).
- **Deterministic Docking & Hungarian Symmetry RMSD:** Reproducible AutoDock Vina v1.2.7 (`--seed 42`) with `exhaustiveness = 12` convergence checks and bipartite element-matched Hungarian RMSD (`scipy.optimize.linear_sum_assignment`).
- **Empirical ChEMBL ROC-AUC Benchmarking:** MW-matched experimental Active vs. Inactive discrimination testing (`CHEMBL5339` for FFAR4; non-covalent `CHEMBL284` `340–480 Da` congeneric panel for DPP-4).
- **Explicit Med-Chem & Translational Filtering:** Lipinski, Veber, Egan/BOILED-Egg GI/BBB windows, Delaney ESOL `LogS`, PAINS, Brenk structural alerts, Ligand Efficiency (`LE`), Lipophilic Ligand Efficiency (`LLE`), and 3D polar H-bond / steric clash (`< 2.2 Å`) classification.

---

## 2. Repository Structure

```text
├── data/                                   # OpenTargets & curated compound benchmark inputs
├── ligands/                                # Protocol A SDF, neutral PDBQT & salivary pH 6.8 PDBQT ligands
├── scripts/
│   ├── structural_prep_utils.py            # Shared altLoc stripper, pH protonation, Hungarian RMSD, 3D contacts & ADMET
│   ├── prepare_all_receptors_strict.py     # Strict PDB receptor auditing & preparation (8G59, 5LGD, 4A5S)
│   ├── build_live_ot_and_benchmarks.py     # Live OpenTargets GraphQL v4 & size-matched ChEMBL ROC pool builder
│   ├── 01_curate_and_admet.py              # Protocol A RDKit/ESOL/BOILED-Egg/PAINS/Brenk profiling
│   ├── 02_prepare_and_dock.py              # Protocol A deterministic pH 6.8 docking (exh=8 + exh=12 convergence)
│   ├── 03_network_and_figures.py           # Protocol A 300-DPI publication figures
│   ├── 04_extract_atomic_contacts.py       # Protocol A 3D H-bond/contact/clash extractor & final CSV builder
│   ├── 05_live_chembl_string_network.py    # Live ChEMBL + OpenTargets + STRING v12.0 PPI & KEGG/GO enrichment
│   └── 06_roc_auc_validation.py            # ChEMBL5339 (FFAR4) & ChEMBL284 (DPP-4) ROC-AUC benchmark runner
├── results/                                # Protocol A CSVs, redocking/ROC JSONs, figures/, network/, poses/
├── manuscript/                             # Protocol A Quarto source, verified references.bib, Solo & Co-Authored HTML/PDF/DOCX
└── protocol_b_hepatic_cd36_dpp4/
    ├── run_protocol_b_pipeline.py          # Protocol B altLoc-clean, plasma pH 7.4 docking (exh=4 + exh=12 refinement)
    ├── 02_figures_and_manuscript_protocol_b.py # Dynamic zero-hardcoding figure + Solo-Author DOCX/HTML/PDF builder
    ├── ligands/                            # Protocol B SDF & plasma pH 7.4 PDBQT ligands
    ├── results/                            # Complete 24-compound CSV, convergence CSV, redocking JSON, live network/
    ├── figures/                            # 6 publication figures (300 DPI) including ChEMBL284 ROC plot
    └── manuscript/                         # Protocol B Solo-Author Word (.docx), HTML, and PDF
```

---

## 3. Verifiable Computational Findings

### Protocol A — Lingual FFAR4 (`8G59`) + CD36 (`5LGD`) at Salivary pH 6.8
- **Redocking & Convergence:** Native `TUG-891` (`8G59`) achieved **`-11.58 kcal/mol`**, **`Hungarian RMSD = 0.864 Å`**, with `0.000–0.069 Å` pose drift between `exh=8` and `exh=12`. Native `PLM-512` (`5LGD`) docked at `-6.13 kcal/mol` (`Hungarian RMSD = 3.098 Å`).
- **Live ChEMBL + OpenTargets + STRING v12.0:** `228` ChEMBL (`<= 10 uM`) targets overlapped `46 / 894` OpenTargets genes (`CD36 = false`, `FFAR4 = true [TUG-891 only]`); `139` STRING edges enriched `PPAR (hsa03320 FDR = 1.81e-05)`, `AMPK (hsa04152 FDR = 1.00e-03)`, and `MASLD (hsa04932 FDR = 2.40e-03)`.
- **ChEMBL5339 ROC-AUC Benchmark (`n = 15` Actives vs `15` Inactives, `p_MW = 0.941`):** `ROC-AUC = 0.733 (p = 0.0155)` by affinity and **`ROC-AUC = 0.756 (p = 0.0090)`** by Ligand Efficiency (`EF20% = 1.67x`).
- **Lead Stratification & Steric Clash Disproof:** Rigid triterpenoids (`Ursolic`, `Oleanolic`, `Corosolic acid`, `Gymnemagenin`) clashed against `Phe115` (`1.70–2.05 Å`) and `Ile280` (`2.17 Å`) inside `8G59` (`-1.12 to +1.49 kcal/mol`, `LE <= 0.034`), whereas PAINS-free Tier-1A leads led by **Naringenin** (`FFAR4 -9.27, LE 0.464, LLE 4.29; CD36 -7.25; Dual -8.46`), **Resveratrol** (`-8.53, LE 0.502`), **Phloretin** (`-8.24, LE 0.412, Thr119 2.72 Å`), **Pinolenic acid** (`-7.92`), **6-Gingerol** (`-7.85`), **Spilanthol** (`-7.65, Thr119 2.94 Å`), and **Pellitorine** (`-7.38`) exhibited zero clashes.

### Protocol B — Hepatic CD36 (`5LGD`) + DPP-4 (`4A5S`) at Plasma pH 7.4
- **Redocking & Convergence:** Native `N7F-901` (`4A5S`) achieved **`-10.71 kcal/mol`**, **`Hungarian RMSD = 0.780 Å`**, restoring H-bonds to `GLU206 (2.76 Å)`, `GLU205 (2.80 Å)`, `TYR662 (2.88 Å)`, `TYR631 (3.17 Å)`, and `TYR547 (3.36 Å)`. Exhaustiveness refinement (`exh=12`) confirmed stability (`max |ΔΔG| = 0.020 kcal/mol`, `max pose RMSD = 0.134 Å`).
- **Live ChEMBL + OpenTargets + STRING v12.0:** `242` ChEMBL targets overlapped `55 / 894` OpenTargets genes (`DPP4 = true`, `CD36 = false`); `145` STRING edges enriched `Adipocytokine (hsa04920 FDR = 2.42e-07)`, `Insulin resistance (hsa04931 FDR = 4.99e-05)`, `PPAR (hsa03320 FDR = 3.30e-04)`, `AMPK (hsa04152 FDR = 7.20e-04)`, **`Ferroptosis (hsa04216 FDR = 7.20e-04)`**, and **`MASLD (hsa04932 FDR = 1.50e-03)`**.
- **ChEMBL284 Non-Covalent ROC-AUC Benchmark (`n = 15` Actives vs `14` Congeneric Inactives, `340–480 Da`):** Raw Vina yielded `ROC-AUC = 0.538 (p = 0.3717)` vs **`ROC-AUC = 0.638 (p = 0.1067)`** for Ligand Efficiency (`EF20% = 1.16x`).
- **Honest Med-Chem Stratification:** Tier 1A PAINS-free Veber-compliant leads identified **Teneligliptin** (`Dual -8.18; DPP-4 -8.99, LE 0.300; CD36 -7.38`), **Naringenin** (`Dual -7.72; DPP-4 -8.19, LE 0.409; CD36 -7.25, LE 0.363; 0 PAINS / 0 Brenk`), and **Berberine** (`Dual -7.57`), distinct from catechol-PAINS Tier 1B hits (`Rosmarinic`, `Chlorogenic`, `Luteolin`, `Quercetin`), high-TPSA (`> 155 Å²`) phytosome-dependent Tier 2A flavonolignans (`Silybin A -8.76`, `Silybin B -8.66`, `Isosilybin A -8.55`, `Silychristin -8.46`), and disqualified oligomer `Lithospermic acid` (`LE <= 0.229`).

---

## 4. Reproducing Both Pipelines

```bash
pip install rdkit pandas numpy scipy scikit-learn networkx matplotlib seaborn python-docx requests

# Shared strict receptor prep, live API & ROC benchmarks
python scripts/prepare_all_receptors_strict.py
python scripts/build_live_ot_and_benchmarks.py
python scripts/05_live_chembl_string_network.py
python scripts/06_roc_auc_validation.py

# Protocol A (Lingual FFAR4 + CD36, salivary pH 6.8)
python scripts/01_curate_and_admet.py
python scripts/02_prepare_and_dock.py
python scripts/04_extract_atomic_contacts.py
python scripts/03_network_and_figures.py
python manuscript/render_academic.py

# Protocol B (Hepatic CD36 + DPP-4, plasma pH 7.4)
python protocol_b_hepatic_cd36_dpp4/run_protocol_b_pipeline.py
python protocol_b_hepatic_cd36_dpp4/02_figures_and_manuscript_protocol_b.py
```
