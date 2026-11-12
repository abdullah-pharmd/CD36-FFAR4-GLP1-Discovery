# Dual CD36 & FFAR4 (GPR120) In-Silico Discovery Pipeline
**Reproducible Network Pharmacology, Validated Molecular Docking, and Oromucosal Delivery Rationale for Post-GLP-1 Weight Rebound**

**Author:** Muhammad Abdullah Khan  
**Affiliation:** Faculty of Pharmaceutical Sciences, Government College University Faisalabad (GCUF), Pakistan  
**Contact:** `abdullah.khan.pharmd@gmail.com`

---

## 1. Scientific Rationale & Scope
While injectable GLP-1 receptor agonists (semaglutide, tirzepatide) achieve substantial weight loss, clinical trials document rapid weight regain upon treatment cessation. Recent clinical proof-of-concept has shown that targeting lingual and enteric lipid-sensing receptors—specifically **CD36 (Fatty Acid Translocase)** and **FFAR4 (GPR120)**—stimulates endogenous GLP-1, CCK, and PYY release.

This repository contains an open-source, fully reproducible **computational drug discovery and network pharmacology pipeline** screening 25 natural dietary phytochemicals (alkamides, flavonoids, polyphenols, and pentacyclic triterpenoids) alongside 4 reference standards against human **FFAR4** (`PDB: 8G59`, 2.80 Å Cryo-EM) and human **CD36** (`PDB: 5LGD`, 2.07 Å X-ray).

---

## 2. Repository Structure

```text
├── data/                         # Input chemical libraries and disease-gene sets
├── structures/                   # Cleaned receptor PDB/PDBQT files (8G59, 5LGD, 8ID6) & native ligands
│   └── ligands_pdbqt/            # 3D protonated ligand structures (29 compounds)
├── scripts/
│   ├── 01_curate_and_admet.py        # RDKit physicochemical profiling (MW, LogP, TPSA, Lipinski/Veber)
│   ├── 02_network_pharmacology.py    # STRING v12 PPI construction & NetworkX hub centrality
│   ├── 03_prepare_and_dock.py        # AutoDock Vina v1.2.7 execution & TUG-891 redocking control
│   ├── 04_extract_atomic_contacts.py # 3D Euclidean heavy-atom contact mapping (<= 4.0 Å)
│   └── generate_new_figures.py       # Publication figure generation (300 DPI)
├── results/
│   ├── curated_29_compounds.csv          # RDKit descriptors & drug-likeness table
│   ├── all_docking_affinities.csv        # Complete dual-receptor AutoDock Vina scores
│   ├── docked_pose_residue_contacts.csv  # Exact interatomic residue distances (<= 4.0 Å)
│   ├── network/                          # PPI topology metrics, KEGG enrichment & plots
│   ├── poses/                            # Docked 3D ligand poses (.pdbqt) for CD36 and FFAR4
│   └── figures/                          # 2D chemical structures, binding plots & redocking validation
└── manuscript/                   # Quarto source (.qmd), BibTeX references, and compiled PDF/DOCX
```

---

## 3. Key Computational Findings (100% Verifiable from Repository Files)

1. **Crystallographic Redocking Validation (`PDB: 8G59`):**
   - Redocking native co-crystallized agonist **TUG-891** into the FFAR4 orthosteric pocket achieved an affinity of **`-11.73 kcal/mol`** and a heavy-atom **`RMSD = 1.530 Å`** (`< 2.0 Å` validation threshold).
2. **Dual CD36 + FFAR4 Agonist Candidates (Alkamides & Polyphenols):**
   - **Piperine:** FFAR4 `-10.34 kcal/mol` (`Asp208 2.90 Å`, `Thr119 3.13 Å`, `Trp198 3.23 Å`) | CD36 `-5.46 kcal/mol`
   - **Chlorogenic acid:** FFAR4 `-9.62 kcal/mol` (`Thr119 2.81 Å`, `Asp208 2.97 Å`) | CD36 `-5.97 kcal/mol`
   - **Naringenin:** FFAR4 `-9.53 kcal/mol` (`Ile284 3.09 Å`, `Thr119 3.38 Å`) | CD36 `-6.15 kcal/mol`
   - **Quercetin:** FFAR4 `-9.30 kcal/mol` (`Ile284 3.24 Å`, `Phe211 3.35 Å`) | CD36 `-5.94 kcal/mol`
   - **Pellitorine:** FFAR4 `-7.62 kcal/mol` (`Thr119 3.35 Å`, `Asp208 3.44 Å`) | CD36 `-5.52 kcal/mol`
3. **Selective CD36 Entrance Tunnel Blockers (Pentacyclic Triterpenoids):**
   - Rigid pentacyclic triterpenoids (**Glycyrrhetinic acid** `-6.51 kcal/mol`, **Corosolic acid** `-6.34 kcal/mol`, **Maslinic acid** `-6.17 kcal/mol`, **Ursolic acid** `-6.13 kcal/mol`) selectively cap the CD36 entrance tunnel (`Lys166 3.33 Å`, `Asn139 2.94 Å`, `Arg173 2.96 Å`) while being sterically excluded from FFAR4's narrow transmembrane pocket (`+0.24` to `-1.56 kcal/mol`).

---

## 4. Reproducing the Analysis

```bash
pip install rdkit pandas numpy networkx matplotlib seaborn requests
python scripts/01_curate_and_admet.py
python scripts/02_network_pharmacology.py
python scripts/04_extract_atomic_contacts.py
```
