# CD36 + FFAR4 Dual-Target Molecular Docking Screening Results

**Software:** AutoDock Vina v1.2.7 | **Search Space:** Grid 22×22×22 Å (FFAR4) & 24×24×24 Å (CD36) | **Exhaustiveness:** 4
**Hardware:** Core i5 Laptop (2 workers × 2 threads) | **Execution Time:** 67.1 min

## 1. Redocking Validation Control (TUG-891 into 8G59)

- **Crystallized Reference Pose:** `structures/8G59_TUG891.pdb` (YN9, 27 heavy atoms)
- **Redocked Mode 1 Binding Affinity:** **-11.93 kcal/mol**
- **Heavy-Atom RMSD vs. Cryo-EM Pose:** **1.530 Å** (Criterion: RMSD < 2.0 Å — **PASSED**)
- **Orthosteric Pocket Residues Engaged:** Arg99, Thr119, Trp198, Phe211, Ile280, Asn291

## 2. Complete Dual-Receptor Docking Leaderboard

| Dual Rank | Compound Name | Role | Chemical Class | MW (Da) | LogP | FFAR4 Affinity (kcal/mol) | CD36 Affinity (kcal/mol) | Mean Dual Affinity (kcal/mol) | Ro5 Pass |
|:---:|:---|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | **Gymnemic acid I** | Test Candidate (CD36/FFAR4) | Triterpenoid Glucuronide Saponin | 806.99 | 3.03 | **-12.5** | **-8.0** | **-10.22** | No |
| 2 | **TUG-891** | Positive Control (FFAR4) | Synthetic Agonist (Reference Control) | 364.42 | 5.4 | **-11.7** | **-6.1** | **-8.87** | Yes |
| 3 | **Piperine** | Test Candidate (Lingual CD36/FFAR4) | Piperidine Alkaloid | 285.34 | 3.0 | **-10.3** | **-5.5** | **-7.90** | Yes |
| 4 | **Naringenin** | Test Candidate (CD36/FFAR4) | Flavanone | 272.26 | 2.51 | **-9.5** | **-6.1** | **-7.84** | Yes |
| 5 | **Chlorogenic acid** | Test Candidate (CD36/FFAR4) | Caffeoylquinic Acid Ester | 354.31 | -0.65 | **-9.6** | **-6.0** | **-7.79** | Yes |
| 6 | **Quercetin** | Test Candidate (CD36/FFAR4) | Flavonol | 302.24 | 1.99 | **-9.3** | **-5.9** | **-7.62** | Yes |
| 7 | **Sulfo-N-succinimidyl oleate (SSO)** | Positive Control (CD36) | Covalent Translocation Blocker (Reference Control) | 459.61 | 4.5 | **-9.3** | **-5.8** | **-7.58** | Yes |
| 8 | **Nitro-oleic acid (NO2-OA)** | Positive Control (CD36) | Electrophilic Nitro-Fatty Acid | 327.46 | 5.36 | **-8.1** | **-6.5** | **-7.30** | Yes |
| 9 | **Berberine** | Test Candidate (CD36/FFAR4) | Quaternary Isoquinoline Alkaloid | 336.37 | 3.1 | **-8.7** | **-5.6** | **-7.16** | Yes |
| 10 | **Curcumin** | Test Candidate (CD36/FFAR4) | Curcuminoid Polyphenol | 368.39 | 3.37 | **-8.9** | **-5.2** | **-7.03** | Yes |
| 11 | **Fraxin** | Test Candidate (CD36/FFAR4) | Coumarin Glucoside | 370.31 | -1.31 | **-7.5** | **-6.4** | **-6.91** | Yes |
| 12 | **Phloretin** | Test Candidate (CD36/FFAR4) | Dihydrochalcone | 274.27 | 2.32 | **-8.5** | **-5.3** | **-6.90** | Yes |
| 13 | **Resveratrol** | Test Candidate (CD36/FFAR4) | Stilbenoid | 228.25 | 2.97 | **-8.8** | **-5.0** | **-6.89** | Yes |
| 14 | **6-Gingerol** | Test Candidate (CD36/FFAR4) | Gingerol Phenol | 294.39 | 3.23 | **-8.2** | **-5.4** | **-6.80** | Yes |
| 15 | **6-Shogaol** | Test Candidate (CD36/FFAR4) | Dehydrated Gingerol Phenol | 276.38 | 4.04 | **-8.1** | **-5.3** | **-6.74** | Yes |
| 16 | **Capsaicin** | Test Candidate (Lingual CD36/FFAR4) | Capsaicinoid | 305.42 | 3.79 | **-8.8** | **-4.4** | **-6.63** | Yes |
| 17 | **Pellitorine** | Test Candidate (Lingual CD36/FFAR4) | Isobutylamide (Alkamide) | 223.36 | 3.45 | **-7.6** | **-5.5** | **-6.57** | Yes |
| 18 | **Hydroxy-alpha-sanshool** | Test Candidate (Lingual CD36/FFAR4) | Polyunsaturated Alkamide | 263.38 | 2.9 | **-8.3** | **-4.7** | **-6.47** | Yes |
| 19 | **Spilanthol** | Test Candidate (Lingual CD36/FFAR4) | Isobutylamide (Alkamide) | 221.34 | 3.23 | **-7.9** | **-3.9** | **-5.88** | Yes |
| 20 | **Oleic acid** | Positive Control (FFAR4 / Natural Ligand) | Endogenous Monounsaturated Fatty Acid | 282.47 | 6.11 | **-7.3** | **-3.4** | **-5.34** | Yes |
| 21 | **Ursolic acid** | Test Candidate (CD36/FFAR4) | Pentacyclic Triterpenoid (Ursane) | 456.71 | 7.09 | **-1.6** | **-6.1** | **-3.85** | Yes |
| 22 | **Maslinic acid** | Test Candidate (CD36/FFAR4) | Pentacyclic Triterpenoid (Oleanane) | 472.71 | 6.2 | **-1.2** | **-6.2** | **-3.67** | Yes |
| 23 | **Glycyrrhetinic acid** | Test Candidate (CD36/FFAR4) | Pentacyclic Triterpenoid (Oleanane) | 470.69 | 6.41 | **-0.7** | **-6.5** | **-3.61** | Yes |
| 24 | **Oleanolic acid** | Test Candidate (CD36/FFAR4) | Pentacyclic Triterpenoid (Oleanane) | 456.71 | 7.23 | **-0.4** | **-6.1** | **-3.28** | Yes |
| 25 | **Asiatic acid** | Test Candidate (CD36/FFAR4) | Pentacyclic Triterpenoid (Ursane) | 488.71 | 5.03 | **-0.2** | **-5.9** | **-3.09** | Yes |
| 26 | **Corosolic acid** | Test Candidate (CD36/FFAR4) | Pentacyclic Triterpenoid (Ursane) | 472.71 | 6.06 | **0.2** | **-6.3** | **-3.05** | Yes |
| 27 | **Betulinic acid** | Test Candidate (CD36/FFAR4) | Pentacyclic Triterpenoid (Lupane) | 456.71 | 7.09 | **-0.1** | **-5.5** | **-2.81** | Yes |
| 28 | **Punicalagin** | Test Candidate (CD36/FFAR4) | Ellagitannin Macromolecule | 1084.72 | 1.8 | **3.5** | **-8.7** | **-2.59** | No |
| 29 | **Gymnemagenin** | Test Candidate (CD36/FFAR4) | Triterpenoid Aglycone (Sapogenin) | 506.72 | 3.03 | **1.0** | **-5.6** | **-2.28** | No |

## 3. High-Priority Dual Lead Scaffolds

### Lead #1: Gymnemic acid I (Triterpenoid Glucuronide Saponin)
- **Dual Affinity:** Mean = **-10.22 kcal/mol** (FFAR4: **-12.5**, CD36: **-8.0**)
- **Physicochemical:** MW = 806.99 Da, LogP = 3.03, Lipinski Ro5 Pass: **No**
- **Mechanistic Role:** Gymnema sylvestre; blocks sweet/lipid lingual taste receptors, stimulates GLP-1

### Lead #2: Piperine (Piperidine Alkaloid)
- **Dual Affinity:** Mean = **-7.90 kcal/mol** (FFAR4: **-10.3**, CD36: **-5.5**)
- **Physicochemical:** MW = 285.34 Da, LogP = 3.0, Lipinski Ro5 Pass: **Yes**
- **Mechanistic Role:** Piper nigrum (Black pepper); bioenhancer, TRPA1/TRPV1 agonist, stimulates satiety

### Lead #3: Naringenin (Flavanone)
- **Dual Affinity:** Mean = **-7.84 kcal/mol** (FFAR4: **-9.5**, CD36: **-6.1**)
- **Physicochemical:** MW = 272.26 Da, LogP = 2.51, Lipinski Ro5 Pass: **Yes**
- **Mechanistic Role:** Citrus fruits; activates PPAR-alpha/gamma, suppresses hepatic lipid accumulation

### Lead #4: Chlorogenic acid (Caffeoylquinic Acid Ester)
- **Dual Affinity:** Mean = **-7.79 kcal/mol** (FFAR4: **-9.6**, CD36: **-6.0**)
- **Physicochemical:** MW = 354.31 Da, LogP = -0.65, Lipinski Ro5 Pass: **Yes**
- **Mechanistic Role:** Green coffee / fruits; delays glucose absorption, stimulates GLP-1 release

### Lead #5: Quercetin (Flavonol)
- **Dual Affinity:** Mean = **-7.62 kcal/mol** (FFAR4: **-9.3**, CD36: **-5.9**)
- **Physicochemical:** MW = 302.24 Da, LogP = 1.99, Lipinski Ro5 Pass: **Yes**
- **Mechanistic Role:** Onion/apples; antioxidant, enhances intestinal GLP-1 secretion via GPR40/FFAR1

