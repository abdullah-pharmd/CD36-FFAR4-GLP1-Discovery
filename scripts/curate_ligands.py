import os
import sys
import time
import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Crippen

DATA_DIR = r"C:\Users\Abdullah\.gemini\antigravity\scratch\cd36_ffar4_glp1_study\data"
STRUCTURES_DIR = r"C:\Users\Abdullah\.gemini\antigravity\scratch\cd36_ffar4_glp1_study\structures"
SDF_DIR = os.path.join(STRUCTURES_DIR, "ligands_sdf")

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(SDF_DIR, exist_ok=True)

# Curated panel of 29 molecules: 4 controls (2 FFAR4, 2 CD36) + 25 test phytochemicals / bioactive scaffolds
COMPOUND_SPECS = [
    # Positive Controls - FFAR4
    {
        "name": "TUG-891",
        "cid": 57522038,
        "class": "Synthetic Agonist (Reference Control)",
        "role": "Positive Control (FFAR4)",
        "source_mechanism": "Selective FFAR4 (GPR120) agonist, co-crystallized in PDB 8G59"
    },
    {
        "name": "Oleic acid",
        "cid": 445639,
        "class": "Endogenous Monounsaturated Fatty Acid",
        "role": "Positive Control (FFAR4 / Natural Ligand)",
        "source_mechanism": "Natural FFAR4 agonist, co-crystallized in PDB 8ID6"
    },
    # Positive Controls - CD36
    {
        "name": "Sulfo-N-succinimidyl oleate (SSO)",
        "cid": 46782948,
        "class": "Covalent Translocation Blocker (Reference Control)",
        "role": "Positive Control (CD36)",
        "source_mechanism": "Irreversible CD36 tunnel blocker modifying Lys164, inhibits LCFA uptake"
    },
    {
        "name": "Nitro-oleic acid (NO2-OA)",
        "cid": 86656470,
        "class": "Electrophilic Nitro-Fatty Acid",
        "role": "Positive Control (CD36)",
        "source_mechanism": "Endogenous electrophilic lipid, binds Lys164/Lys166 CD36 entrance tunnel"
    },
    # Test Ligands - Pentacyclic Triterpenoids
    {
        "name": "Oleanolic acid",
        "cid": 10494,
        "class": "Pentacyclic Triterpenoid (Oleanane)",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Olive/clove bioactive; activates GPCR lipid sensing, enhances GLP-1 secretion"
    },
    {
        "name": "Ursolic acid",
        "cid": 64945,
        "class": "Pentacyclic Triterpenoid (Ursane)",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Apple peel/rosemary; anti-sarcopenic, preserves skeletal muscle mass via Akt/mTOR"
    },
    {
        "name": "Corosolic acid",
        "cid": 6918774,
        "class": "Pentacyclic Triterpenoid (Ursane)",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Banaba leaf (Lagerstroemia speciosa); GLUT4 activator, insulin-sensitizing"
    },
    {
        "name": "Maslinic acid",
        "cid": 73659,
        "class": "Pentacyclic Triterpenoid (Oleanane)",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Crataegus/olive pomace; anti-inflammatory, lipid metabolism regulator"
    },
    {
        "name": "Betulinic acid",
        "cid": 64971,
        "class": "Pentacyclic Triterpenoid (Lupane)",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Birch bark; AMPK activator, mitigates adipogenesis"
    },
    {
        "name": "Asiatic acid",
        "cid": 119034,
        "class": "Pentacyclic Triterpenoid (Ursane)",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Centella asiatica; hepatoprotective, suppresses lipotoxicity and fibrosis"
    },
    {
        "name": "Glycyrrhetinic acid",
        "cid": 10114,
        "class": "Pentacyclic Triterpenoid (Oleanane)",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Glycyrrhiza glabra (licorice); 11beta-HSD1 modulator, lipid-clearing"
    },
    {
        "name": "Gymnemic acid I",
        "cid": 11953919,
        "class": "Triterpenoid Glucuronide Saponin",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Gymnema sylvestre; blocks sweet/lipid lingual taste receptors, stimulates GLP-1"
    },
    {
        "name": "Gymnemagenin",
        "cid": 10051937,
        "class": "Triterpenoid Aglycone (Sapogenin)",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Gymnema sylvestre core aglycone; high lipophilic membrane permeability"
    },
    # Test Ligands - Alkamides & Lingual Sensory Modulators
    {
        "name": "Spilanthol",
        "cid": 5353001,
        "class": "Isobutylamide (Alkamide)",
        "role": "Test Candidate (Lingual CD36/FFAR4)",
        "source_mechanism": "Acmella oleracea (Jambu); fat-taste and trigeminal sensory agonist, increases salivation"
    },
    {
        "name": "Pellitorine",
        "cid": 5318516,
        "class": "Isobutylamide (Alkamide)",
        "role": "Test Candidate (Lingual CD36/FFAR4)",
        "source_mechanism": "Anacyclus pyrethrum (Akarkara); pungent alkamide, lingual TRP/lipid-taste stimulant"
    },
    {
        "name": "Hydroxy-alpha-sanshool",
        "cid": 10084135,
        "class": "Polyunsaturated Alkamide",
        "role": "Test Candidate (Lingual CD36/FFAR4)",
        "source_mechanism": "Zanthoxylum piperitum (Sichuan pepper); tingles lingual somatosensory/taste papillae"
    },
    {
        "name": "Piperine",
        "cid": 638024,
        "class": "Piperidine Alkaloid",
        "role": "Test Candidate (Lingual CD36/FFAR4)",
        "source_mechanism": "Piper nigrum (Black pepper); bioenhancer, TRPA1/TRPV1 agonist, stimulates satiety"
    },
    {
        "name": "Capsaicin",
        "cid": 1548943,
        "class": "Capsaicinoid",
        "role": "Test Candidate (Lingual CD36/FFAR4)",
        "source_mechanism": "Capsicum annuum; TRPV1 agonist, stimulates enteric GLP-1 secretion and thermogenesis"
    },
    # Test Ligands - Phenolic & Phenylpropanoids
    {
        "name": "6-Gingerol",
        "cid": 442793,
        "class": "Gingerol Phenol",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Zingiber officinale; stimulates gut incretin secretion, improves insulin sensitivity"
    },
    {
        "name": "6-Shogaol",
        "cid": 5281794,
        "class": "Dehydrated Gingerol Phenol",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Dried ginger; PPAR-gamma and AMPK agonist, anti-inflammatory"
    },
    {
        "name": "Curcumin",
        "cid": 969516,
        "class": "Curcuminoid Polyphenol",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Curcuma longa; potent Nrf2/AMPK activator, downregulates hepatic CD36 in MASLD"
    },
    {
        "name": "Berberine",
        "cid": 2353,
        "class": "Quaternary Isoquinoline Alkaloid",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Berberis aristata; AMPK activator, induces enteric L-cell GLP-1 secretion"
    },
    {
        "name": "Chlorogenic acid",
        "cid": 1794427,
        "class": "Caffeoylquinic Acid Ester",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Green coffee / fruits; delays glucose absorption, stimulates GLP-1 release"
    },
    {
        "name": "Quercetin",
        "cid": 5280343,
        "class": "Flavonol",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Onion/apples; antioxidant, enhances intestinal GLP-1 secretion via GPR40/FFAR1"
    },
    {
        "name": "Resveratrol",
        "cid": 445154,
        "class": "Stilbenoid",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Grapes/berries; SIRT1/AMPK activator, preserves muscle mitochondrial integrity"
    },
    {
        "name": "Naringenin",
        "cid": 439246,
        "class": "Flavanone",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Citrus fruits; activates PPAR-alpha/gamma, suppresses hepatic lipid accumulation"
    },
    {
        "name": "Punicalagin",
        "cid": 16129719,
        "class": "Ellagitannin Macromolecule",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Punica granatum (Pomegranate); prebiotic gut metabolite precursor (urolithins)"
    },
    {
        "name": "Fraxin",
        "cid": 5273568,
        "class": "Coumarin Glucoside",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Fraxinus rhynchophylla; antioxidant, protects against hepatic lipotoxicity"
    },
    {
        "name": "Phloretin",
        "cid": 4788,
        "class": "Dihydrochalcone",
        "role": "Test Candidate (CD36/FFAR4)",
        "source_mechanism": "Apple leaves; GLUT2 inhibitor, modulates mucosal lipid and sugar absorption"
    }
]

def fetch_pubchem_record(cid):
    """Fetch JSON details from PubChem PUG REST API."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/ConnectivitySMILES,SMILES,MolecularFormula,MolecularWeight,InChIKey,IUPACName,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount/JSON"
    headers = {"User-Agent": "Antigravity-BioChemPipeline/1.0 (Academic Research; solo-PharmD)"}
    try:
        r = requests.get(url, headers=headers, timeout=15)
        if r.status_code == 200:
            data = r.json()
            props = data["PropertyTable"]["Properties"][0]
            return props
    except Exception as e:
        print(f"Error fetching CID {cid} properties: {e}")
    return {}

def download_sdf(cid, name):
    """Download 3D SDF (or fallback 2D) from PubChem."""
    safe_name = name.split("(")[0].strip().replace(" ", "_").replace("-", "_")
    sdf_path = os.path.join(SDF_DIR, f"{cid}_{safe_name}.sdf")
    if os.path.exists(sdf_path) and os.path.getsize(sdf_path) > 100:
        return sdf_path

    # Try 3D first
    url_3d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    headers = {"User-Agent": "Antigravity-BioChemPipeline/1.0"}
    try:
        r = requests.get(url_3d, headers=headers, timeout=15)
        if r.status_code == 200 and len(r.content) > 200:
            with open(sdf_path, 'wb') as f:
                f.write(r.content)
            return sdf_path
    except Exception:
        pass

    # Fallback to 2D
    url_2d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=2d"
    try:
        r = requests.get(url_2d, headers=headers, timeout=15)
        if r.status_code == 200 and len(r.content) > 100:
            with open(sdf_path, 'wb') as f:
                f.write(r.content)
            return sdf_path
    except Exception as e:
        print(f"Error downloading SDF for CID {cid}: {e}")
    return None

def compute_rdkit_properties(smiles):
    """Compute exact physicochemical properties and Lipinski violations using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    arom_rings = Lipinski.NumAromaticRings(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    fraction_csp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    
    # Lipinski's Rule of 5:
    # 1. MW <= 500
    # 2. LogP <= 5.0
    # 3. HBD <= 5
    # 4. HBA <= 10
    violations = 0
    violation_reasons = []
    if mw > 500:
        violations += 1
        violation_reasons.append("MW>500")
    if logp > 5.0:
        violations += 1
        violation_reasons.append("LogP>5")
    if hbd > 5:
        violations += 1
        violation_reasons.append("HBD>5")
    if hba > 10:
        violations += 1
        violation_reasons.append("HBA>10")
        
    lipinski_violations_str = "; ".join(violation_reasons) if violation_reasons else "None"
    
    # Standard medicinal chemistry convention: <= 1 violation is Lipinski compliant (Ro5 pass)
    lipinski_compliant = "Yes" if violations <= 1 else "No"
    strict_ro5_compliant = "Yes" if violations == 0 else "No"
    
    # Veber filter: Rotatable bonds <= 10 and TPSA <= 140 A^2
    veber_violations = 0
    veber_reasons = []
    if rotb > 10:
        veber_violations += 1
        veber_reasons.append("RotB>10")
    if tpsa > 140:
        veber_violations += 1
        veber_reasons.append("TPSA>140")
    veber_compliant = "Yes" if veber_violations == 0 else "No"
    
    return {
        "MW_RDKit": round(mw, 2),
        "LogP_RDKit": round(logp, 2),
        "TPSA_RDKit": round(tpsa, 2),
        "HBD_RDKit": int(hbd),
        "HBA_RDKit": int(hba),
        "Rotatable_Bonds": int(rotb),
        "Aromatic_Rings": int(arom_rings),
        "Heavy_Atoms": int(heavy_atoms),
        "Fraction_CSP3": round(fraction_csp3, 3),
        "Lipinski_Violations_Count": violations,
        "Lipinski_Violations_Detail": lipinski_violations_str,
        "Lipinski_Compliant": lipinski_compliant,
        "Strict_Ro5_Compliant": strict_ro5_compliant,
        "Veber_Compliant": veber_compliant
    }

def main():
    print(f"=== CURATING 29-LIGAND PANEL (25 TEST + 4 CONTROLS) ===")
    print(f"Target directory: {DATA_DIR}\n")
    
    records = []
    for idx, spec in enumerate(COMPOUND_SPECS, 1):
        name = spec["name"]
        cid = spec["cid"]
        print(f"[{idx:02d}/29] Fetching CID {cid}: {name}...", end=" ", flush=True)
        
        props = fetch_pubchem_record(cid)
        canonical_smiles = props.get("ConnectivitySMILES") or props.get("SMILES")
        isomeric_smiles = props.get("SMILES") or canonical_smiles
        
        if not canonical_smiles:
            print("FAILED TO FETCH SMILES!")
            continue
            
        rd_props = compute_rdkit_properties(isomeric_smiles or canonical_smiles)
        sdf_path = download_sdf(cid, name)
        
        row = {
            "ID": f"LIG_{idx:02d}",
            "Compound_Name": name,
            "Role": spec["role"],
            "Chemical_Class": spec["class"],
            "PubChem_CID": cid,
            "Molecular_Formula": props.get("MolecularFormula", ""),
            "Canonical_SMILES": canonical_smiles,
            "Isomeric_SMILES": isomeric_smiles,
            "InChIKey": props.get("InChIKey", ""),
            "IUPAC_Name": props.get("IUPACName", ""),
            "PubChem_XLogP": props.get("XLogP", ""),
            "PubChem_TPSA": props.get("TPSA", ""),
            "MW": rd_props["MW_RDKit"] if rd_props else props.get("MolecularWeight", ""),
            "LogP": rd_props["LogP_RDKit"] if rd_props else "",
            "TPSA": rd_props["TPSA_RDKit"] if rd_props else "",
            "HBD": rd_props["HBD_RDKit"] if rd_props else "",
            "HBA": rd_props["HBA_RDKit"] if rd_props else "",
            "Rotatable_Bonds": rd_props["Rotatable_Bonds"] if rd_props else "",
            "Aromatic_Rings": rd_props["Aromatic_Rings"] if rd_props else "",
            "Heavy_Atoms": rd_props["Heavy_Atoms"] if rd_props else "",
            "Fraction_CSP3": rd_props["Fraction_CSP3"] if rd_props else "",
            "Lipinski_Violations": rd_props["Lipinski_Violations_Count"] if rd_props else "",
            "Lipinski_Violations_Detail": rd_props["Lipinski_Violations_Detail"] if rd_props else "",
            "Lipinski_Compliant": rd_props["Lipinski_Compliant"] if rd_props else "",
            "Strict_Ro5": rd_props["Strict_Ro5_Compliant"] if rd_props else "",
            "Veber_Compliant": rd_props["Veber_Compliant"] if rd_props else "",
            "SDF_File": os.path.basename(sdf_path) if sdf_path else "",
            "Pharmacological_Mechanism": spec["source_mechanism"]
        }
        records.append(row)
        print(f"OK (MW={row['MW']}, LogP={row['LogP']}, Ro5={row['Lipinski_Compliant']})")
        time.sleep(0.2)
        
    df = pd.DataFrame(records)
    out_csv = os.path.join(DATA_DIR, "01_ligands.csv")
    df.to_csv(out_csv, index=False)
    print(f"\n[SUCCESS] Successfully saved {len(df)} compounds to: {out_csv}")
    
    # Save a clean TSV and an aligned summary table
    summary_md = os.path.join(DATA_DIR, "01_ligands_summary.md")
    with open(summary_md, "w", encoding="utf-8") as f:
        f.write("# CD36 + FFAR4 Curated Ligand Panel & Physicochemical Characterization\n\n")
        f.write(f"**Total Compounds:** {len(df)} (4 Reference Controls + 25 Bioactive Test Compounds)\n\n")
        f.write(f"- **Lipinski Compliant (<=1 violation):** {sum(df['Lipinski_Compliant'] == 'Yes')}/{len(df)} ({(sum(df['Lipinski_Compliant'] == 'Yes')/len(df)*100):.1f}%)\n")
        f.write(f"- **Strict Lipinski Pass (0 violations):** {sum(df['Strict_Ro5'] == 'Yes')}/{len(df)} ({(sum(df['Strict_Ro5'] == 'Yes')/len(df)*100):.1f}%)\n")
        f.write(f"- **Veber Compliant (RotB <= 10, TPSA <= 140 Å²):** {sum(df['Veber_Compliant'] == 'Yes')}/{len(df)} ({(sum(df['Veber_Compliant'] == 'Yes')/len(df)*100):.1f}%)\n\n")
        
        f.write("### Complete Compound Characterization Table\n\n")
        f.write("| ID | Compound Name | Role | Class | CID | MW (Da) | LogP | TPSA (Å²) | HBD | HBA | RotB | Ro5 Detail | Ro5 Pass |\n")
        f.write("|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|\n")
        for _, r in df.iterrows():
            f.write(f"| {r['ID']} | **{r['Compound_Name']}** | {r['Role']} | {r['Chemical_Class']} | {r['PubChem_CID']} | {r['MW']} | {r['LogP']} | {r['TPSA']} | {r['HBD']} | {r['HBA']} | {r['Rotatable_Bonds']} | {r['Lipinski_Violations_Detail']} | **{r['Lipinski_Compliant']}** |\n")
    print(f"[SUCCESS] Saved summary markdown table to: {summary_md}")

if __name__ == "__main__":
    main()
