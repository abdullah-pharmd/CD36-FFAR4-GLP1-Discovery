# Reproducible Computational Pharmacology Suite: Hepatic CD36 / Vascular DPP-4 Prioritization (Protocol B) & Exploratory Lingual FFAR4 / CD36 Structural Sensitivity Study (Protocol A)

**Author:** Muhammad Abdullah Khan ([ORCID: 0009-0009-8867-6526](https://orcid.org/0009-0009-8867-6526))  
**Affiliation:** Department of Pharmaceutics, Faculty of Pharmaceutical Sciences, Government College University Faisalabad (GCUF), Pakistan  
**Contact:** `abdullah.khan.pharmd@gmail.com`  
**Repository:** [https://github.com/abdullah-pharmd/CD36-FFAR4-GLP1-Discovery](https://github.com/abdullah-pharmd/CD36-FFAR4-GLP1-Discovery)

> **Scope & Evidence Boundary:** All findings in this repository are **dry-lab computational hypotheses and methodological benchmarks**. No wet-lab IC50, SPR/ITC binding, cellular lipid-uptake, GLP-1 secretion, or clinical outcome assays were performed here. Both studies are independent undergraduate research manuscripts currently in preparation / revision (not yet peer-reviewed).

---

## 1. Repository Architecture & Current Project Status

| Track | Target Axis & Physiological pH | Crystal / Cryo-EM Structures | Status | Primary Scientific Contribution |
| :--- | :--- | :--- | :--- | :--- |
| **Protocol B (Flagship)**<br>`protocol_b_hepatic_cd36_dpp4/` | Hepatic **CD36** lipid-translocation tunnel + vascular **DPP-4** catalytic pocket (**plasma pH 7.4**) | `PDB: 5LGD` (Chain A, `2.07 Å`, 400 res)<br>`PDB: 4A5S` (Chain A, `1.62 Å`, 728 res) | **Solo-Author Manuscript Ready for Submission** *(Molecular Diversity; all 8/8 complexes × 50 ns = 400 ns OpenMM MD + OBC2 MM-GBSA complete)* | Unspiked two-stage ChEMBL–Open Targets–STRING target-gap audit (`0` empirical CD36 hits vs `4` DPP-4 gliptins), continuum PROPKA 3.5 verification, `4A5S` Val546 Conformer A/B & 3-seed stability, N=60 PubChem DPP-4 ROC + scaffold/baseline audits, exact Daina–Zoete BOILED-Egg + OCT1/P-gp transporter gating, and 400 ns explicit-water OpenMM MD + OBC2 MM-GBSA. |
| **Protocol A (Companion)**<br>`scripts/`, `results/`, `manuscript/` | Lingual **FFAR4/GPR120** orthosteric pocket + exploratory **CD36** chaperone docking (**salivary pH 6.8**) | `PDB: 8G59` (Chain R, `2.64 Å`)<br>`PDB: 5LGD` (Chain A, `2.07 Å`) | **Exploratory Structural Sensitivity Study (In Revision — Not Submitted)** | Self-redocking validation (`TUG-891`, `0.864 Å` Hungarian RMSD), `CHEMBL5339` ranking benchmark (`n = 30`), and explicit diagnostic audits exposing `LE ≥ 0.30` zero-specificity threshold failure, `8G59 Phe115` 50/50 alternate-conformer sensitivity, and solubility unit verification. |

---

## 2. Protocol B (Flagship) — Verified Computational Findings (`protocol_b_hepatic_cd36_dpp4/`)

1. **Unspiked Two-Stage Empirical Network Audit (Zero Circularity):**
   - Strict human single-protein ChEMBL curation (`pChEMBL ≥ 5.0` or approved mechanism) across all `24` curated compounds recovered `125` unique targets; intersection with `894` Open Targets GraphQL v4 cardiometabolic genes (`MONDO_0011122`, `MONDO_0005148`, `EFO_1001248`, `EFO_1000653`) yielded **`21` pure empirical overlap proteins** (`44` STRING v12.0 edges, `16` connected nodes, `5` singletons, `9` KEGG pathways at `FDR < 0.05`).
   - Empirical bioassays confirmed DPP-4 potency for `4` clinical gliptins (`Linagliptin`, `Sitagliptin`, `Teneligliptin`, `Vildagliptin`) and **`0` strict quantitative hits for CD36** (`CD36_in_Empirical_ChEMBL_Hits = false`), demonstrating why purely ligand-based network pipelines miss CD36. Stage 2 bridge graphs (`+5` unhit effectors: `ACSL4, CD36, GCG, GLP1R, GPX4`) are strictly separated for topological visualization and excluded from enrichment testing.
2. **Continuum PROPKA 3.5 Electrostatics (Plasma pH 7.4):**
   - `5LGD` (`113` ionizable groups, `8` pocket anchors): `0` global or pocket protonation flips.
   - `4A5S` (`238` ionizable groups): `8/9` key pocket residues preserved aqueous ionization states (`Glu205/Glu206` anionic, active-site tyrosines neutral); buried catalytic triad `His740` shifted upward (`pKa = 8.00`) via coupling to `Asp708`, transparently documented as non-perturbing for non-covalent leads.
3. **Crystal Redocking, Conformer A/B Sensitivity, Multi-Seed Stability & N=60 DPP-4 Benchmark:**
   - **Native Redocking (`AutoDock Vina v1.2.7`, `exh = 8`, `seed = 42`):** `N7F-901` (`4A5S`) achieved `−10.71 kcal/mol`, **`0.780 Å` symmetry-corrected Hungarian heavy-atom RMSD** (`0.329 Å` CoM shift). Valence-explicit PubChem palmitate (`CID 985`, `TORSDOF 14`) in `5LGD` achieved `−6.13 kcal/mol` with `0.156 Å` perpendicular tunnel-axis error and `3.098 Å` heavy-atom RMSD along 14 flexible alkyl torsions (explicitly motivating dynamic MD validation).
   - **`4A5S Val546` 50/50 Alternate Conformer & 3-Seed Audit:** Fixed-pose auditing across all `24` docked DPP-4 poses showed `0` steric clashes against Conformer B (`min dist = 2.720 Å`; `N7F-901` sits `3.092 Å` from both A and B), verified alongside Conformer B redocking and 3-seed stability (`seeds 42, 137, 2026`).
   - **N=60 PubChem BioAssay (`GeneID: 1803`) Non-Covalent DPP-4 Benchmark (`30` actives ≤ 100 nM vs `30` inactives ≥ 50 µM / HTS-inactive; `p_MW = 0.779`):**
     - Raw Vina `ΔG`: **`ROC-AUC = 0.686`** (`95% CI 0.546–0.813`, `p = 0.0068`, `PR-AUC = 0.711`, `Top-12 = 11/12` actives, `EF20% = 1.83×`).
     - Bemis–Murcko scaffold-collapsed (`51` unique scaffolds: `25/30` active, `26/30` inactive): **`ROC-AUC = 0.665`** (`95% CI 0.508–0.811`, `p = 0.0219`).
     - Beats simple `MW` (`0.485`) and `Heavy_Atoms` (`0.471`) baselines (`p < 0.025` paired bootstrap); `LE` (`0.671`) performs indistinguishably on the MW-matched panel (`ΔAUC = +0.015`, `95% CI −0.118 to +0.149`) while filtering `> 500 Da` oligomers in the 24-compound screen.
     - Tier 1A DPP-4 gate (`ΔG ≤ −8.0 & LE ≥ 0.280`): `Sensitivity = 73.3% (22/30)`, `Specificity = 43.3% (13/30)`.
     - **Disclosed Confounder:** PubChem actives are more conformationally rigid than flexible HTS inactives (`mean RotB 3.20 vs. 4.77`, `p = 0.0001`), coupling into Vina's torsion term (`Eq. 2`) and attenuating static separation within the narrower `RotB 3–5` subset (`n = 38`, `AUC = 0.557`), which directly motivates explicit-water 50 ns OpenMM MD + OBC2 MM-GBSA rescoring.
4. **Pharmacokinetic & Med-Chem Triage (Exact Daina–Zoete BOILED-Egg Ellipses + Transporter Gating):**
   - **Tier 1A Dual-Target Hypotheses (0 PAINS, Veber Pass, High HIA):** `Linagliptin` (`DPP-4 −9.81 / CD36 −8.79`), `Sitagliptin` (`−9.21 / −8.21`, `BBB−`), `Teneligliptin` (`−8.99 / −7.37`), and oral PAINS-free phytochemical **`Naringenin`** (`−8.19 / −7.25`, `LE 0.409 / 0.363`, `ESOL 88.1 µg/mL`, `0 PAINS / 0 Brenk`).
   - **Tier 1B-Cation:** `Berberine` (`−8.39 / −6.76`) explicitly excluded from passive BOILED-Egg absorption due to permanent quaternary nitrogen (`ABCB1/P-gp` intestinal efflux + sinusoidal `OCT1` hepatic import).
   - **Tier 2A (Phytosome-Dependent):** `Silybin A`, `Silybin B`, `Silychristin`, `Isosilybin A` (`0 PAINS`, `TPSA > 155 Å²`).
   - **Disqualified Oligomers:** `Lithospermic acid` & `Salvianolic acid B` (`MW > 500 Da`, `LE < 0.250`).
5. **Completed Explicit-Water 50 ns OpenMM MD + OBC2 MM-GBSA (`protocol_b_hepatic_cd36_dpp4/colab_md/`, 8 complexes = `400 ns` aggregate):**
   - AMBER `ff14SB` + TIP3P (`0.15 M NaCl`) + OpenFF Sage `2.2.0` + AmberTools `sqm` AM1-BCC charges, `4 fs` HMR (`3.0 amu`), `310 K / 1 bar`, `250` equilibrated OBC2 MM-GBSA frames (`25–50 ns`), and Baker–Hubbard (`≤ 3.5 Å, ≥ 120°`) H-bond & water-bridge auditing in exact PDB numbering:
   - **DPP-4 (`4A5S`):** `N7F-901` (`Cα 1.94±0.20 Å, Lig 3.28±0.42 Å, ΔG −61.66±0.38 kcal/mol`), `Teneligliptin` (`1.66±0.10, 4.42±0.69 Å, −26.13±0.45`), `Berberine` (`1.42±0.10, 5.78±2.03 Å, −23.20±0.26`), `Naringenin` (`1.62±0.14`, **tightest ligand RMSD `2.78±0.30 Å`**, `−21.61±0.17`).
   - **Hepatic CD36 (`5LGD`):** `PLM-512` (`1.53±0.14, 2.63±0.52 Å, −42.13±0.44 kcal/mol; ARG337 72.0% HB / GLU335 70.4% WB`), **`Linagliptin`** (`1.66±0.15`, **locked at `1.69±0.25 Å`**, `−34.79±0.27; ARG63 78.8% HB / LYS334 89.6% WB`), `Teneligliptin` (`1.57±0.12, 6.55±0.68 Å, −32.80±0.30; GLU335 32.0% HB / 84.0% WB`), `Naringenin` (`1.25±0.14, 6.19±0.54 Å, −24.87±0.18; GLU335 69.6% HB / TYR149 75.2% WB`).

---

## 3. Protocol A (Companion Exploratory Study — In Revision) — Transparent Audit Status

- **What Holds Up:** Native `TUG-891` self-redocks cleanly in `8G59` (`−11.58 kcal/mol`, `Hungarian RMSD = 0.864 Å`), and labeled `CHEMBL5339` actives (`n = 15`) outrank labeled inactives (`n = 15`) on average (`ROC-AUC = 0.733` raw, `0.756` LE).
- **Documented Methodological Issues Corrected or Under Revision (Not Submitted as Efficacy Claims):**
  1. **Repaired 1,000× Solubility Unit Bug:** `results/admet_and_physicochemical_properties.csv` previously stored `µg/mL` values under a `mg/mL` column header; regenerated with explicit `Estimated_Solubility_mg_mL` and `Estimated_Solubility_ug_mL` columns, unit-conversion invariant checks, isomeric SMILES, full InChIKeys, and exact Daina–Zoete BOILED-Egg ellipses.
  2. **Zero-Specificity Cutoff Retired:** On the 30-compound `CHEMBL5339` panel, `LE ≥ 0.30` admitted `15/15` actives and `15/15` inactives (`0%` specificity; highest LE `0.488` belonged to a labeled inactive), and paired `LE − Raw` AUC gain (`+0.022`, `95% CI −0.156 to +0.198`) is unproven. `"Efficiency-Validated"` wording has been retired from data tables pending held-out scaffold calibration and assay-level provenance rebuilding.
  3. **`8G59 Phe115` 50/50 Alternate-Conformer Sensitivity (Empirically Verified via Conformer B Redocking):** Because `Phe115` in `8G59` has `0.50 / 0.50` occupancy across `altLoc A` and `altLoc B`, fresh redocking into Conformer B ([audit_8g59_phe115_conformer_ab.csv](results/audit_8g59_phe115_conformer_ab.csv)) completely eliminated every `< 2.2 Å` clash for all four triterpenoids and swung their energies by `−5.22 to −7.04 kcal/mol` (`Ursolic acid −1.12 → −7.45`, `Oleanolic acid +0.04 → −7.00`, `Corosolic acid +0.52 → −6.04`, `Gymnemagenin +1.49 → −3.73 kcal/mol`), whereas compact top leads stayed conformer-stable (`TUG-891 −11.58/−11.31`, `Naringenin −9.27/−8.90`, `Phloretin −8.24/−8.48`, `Resveratrol −8.53/−8.26 kcal/mol`). Single-conformer steric clashes are therefore reported as conformer-sensitive rather than biological disproof.
  4. **Biological Scope & Bitter-Taste (`TAS2R`) Boundary:** Rigid docking in active-state `8G59` cannot establish GPCR agonism, reversible CD36 lipid transport, GLP-1 release, or clinical weight maintenance, and oral polyphenols such as naringenin and phloretin carry reported micromolar bitter-taste receptor (`TAS2R`) activity requiring direct sensory/exposure evaluation.

---

## 4. Reproducing the Workspace

Install pinned Python dependencies (`requirements.txt`) plus OpenBabel (`obabel` on `PATH` or via `OBABEL_EXE`) and AutoDock Vina v1.2.7 (`scripts/vina.exe` or system `vina`):

```bash
pip install -r requirements.txt

# Protocol B — Flagship Hepatic CD36 (5LGD) + DPP-4 (4A5S) Pipeline & Audits
python protocol_b_hepatic_cd36_dpp4/03_hardening_audits_protocol_b.py
python protocol_b_hepatic_cd36_dpp4/02_figures_and_manuscript_protocol_b.py
python protocol_b_hepatic_cd36_dpp4/submission_molecular_diversity/02_generate_graphical_abstract.py
python protocol_b_hepatic_cd36_dpp4/submission_molecular_diversity/03_build_submission_package.py
python protocol_b_hepatic_cd36_dpp4/submission_molecular_diversity/01_stylometric_and_ai_tell_audit.py

# Protocol A — ADMET Repair & 8G59 Phe115 Conformer A/B Sensitivity Audit
python scripts/07_audit_8g59_phe115_and_fix_admet.py
```

---

## 5. Authorship & Software Disclosure
All scripts, analyses, and manuscripts in this repository are sole-authored by **Muhammad Abdullah Khan** (independent undergraduate computational research). Open-source scientific tools (RDKit, OpenBabel, PROPKA 3.5, AutoDock Vina v1.2.7, PyMOL Open-Source, OpenMM 8, OpenFF Sage 2.2.0, AmberTools, MDTraj) executed all calculations; generative AI coding and editing assistants were used for script scaffolding, numerical cross-audits, and prose refinement under the author's full verification and responsibility.
