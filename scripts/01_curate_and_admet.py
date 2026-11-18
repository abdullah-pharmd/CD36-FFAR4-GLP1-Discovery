import os
import json
import time
import urllib.request
import urllib.parse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from structural_prep_utils import compute_admet_descriptors

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIG_DIR = os.path.join(BASE_DIR, "ligands")
SDF_DIR = os.path.join(LIG_DIR, "sdf")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

os.makedirs(SDF_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

COMPOUNDS = [
    {"name": "TUG-891", "category": "Positive Control (Co-Crystal FFAR4 Agonist)", "source": "Synthetic Ortho-Biphenyl FFAR4 Agonist (8G59)"},
    {"name": "Oleic acid", "category": "Positive Control (Endogenous FFAR4/CD36 Lipid)", "source": "Dietary Monounsaturated Fatty Acid (8ID6)"},
    {"name": "Palmitic acid", "category": "Positive Control (Co-Crystal CD36 Tunnel Ligand)", "source": "Dietary Saturated Fatty Acid (5LGD Native PLM-512)"},
    {"name": "Oleoylethanolamide", "category": "Positive Control (Endogenous Satiety Mediator)", "source": "Enteroendocrine Lipid Amide (OEA)"},
    {"name": "Spilanthol", "category": "Lingual Sensory Alkamide", "source": "Acmella oleracea (Jambu / Sialogogue)"},
    {"name": "Pellitorine", "category": "Lingual Sensory Alkamide", "source": "Anacyclus pyrethrum (Akarkara)"},
    {"name": "Pinolenic acid", "category": "Dietary PUFA (GLP-1 Secretagogue)", "source": "Pinus koraiensis (Pine Nut Oil)"},
    {"name": "Alpha-Linolenic acid", "category": "Dietary Omega-3 PUFA", "source": "Linum usitatissimum (Flaxseed / Chia)"},
    {"name": "Bovinic acid", "category": "Dietary Conjugated Linoleic Acid (9Z,11E-CLA)", "source": "Ruminant / Plant Seed Lipid"},
    {"name": "10-Hydroxy-cis-12-octadecenoic acid", "category": "Gut/Lingual Microbiota Lipid Metabolite (HYA)", "source": "Lactobacillus Linoleic Acid Metabolite"},
    {"name": "Phloretin", "category": "Dietary Dihydrochalcone", "source": "Malus domestica (Apple Peel)"},
    {"name": "Naringenin", "category": "Citrus Flavanone", "source": "Citrus paradisi (Grapefruit)"},
    {"name": "Resveratrol", "category": "Dietary Stilbenoid", "source": "Vitis vinifera (Red Grape Skin)"},
    {"name": "6-Gingerol", "category": "Dietary Alkanone", "source": "Zingiber officinale (Ginger)"},
    {"name": "Chlorogenic acid", "category": "Dietary Phenolic Ester", "source": "Coffea arabica (Green Coffee)"},
    {"name": "Ferulic acid", "category": "Dietary Hydroxycinnamic Acid", "source": "Cereal Bran / Citrus"},
    {"name": "Luteolin", "category": "Dietary Flavone", "source": "Apium graveolens (Celery / Parsley)"},
    {"name": "Quercetin", "category": "Dietary Flavonol", "source": "Allium cepa (Red Onion / Capers)"},
    {"name": "Corosolic acid", "category": "Pentacyclic Triterpenoid (Size/LE Benchmark)", "source": "Lagerstroemia speciosa (Banaba Leaf)"},
    {"name": "Ursolic acid", "category": "Pentacyclic Triterpenoid (Size/LE Benchmark)", "source": "Salvia rosmarinus (Rosemary / Apple Peel)"},
    {"name": "Oleanolic acid", "category": "Pentacyclic Triterpenoid (Size/LE Benchmark)", "source": "Olea europaea (Olive Leaf)"},
    {"name": "Gymnemagenin", "category": "Triterpenoid Sapogenin (Sweet-Suppressor Aglycone)", "source": "Gymnema sylvestre (Gurmar)"},
]


def fetch_pubchem_data(name):
    encoded = urllib.parse.quote(name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/CanonicalSMILES,IsomericSMILES,MolecularFormula/JSON"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        data = json.loads(resp.read().decode())
    props = data["PropertyTable"]["Properties"][0]
    cid = props["CID"]
    smiles = props.get("IsomericSMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
    formula = props.get("MolecularFormula", "")
    return cid, smiles, formula


def ensure_3d_sdf(cid, name, smiles):
    safe_name = name.replace(" ", "_").replace("-", "_").replace(",", "")
    sdf_path = os.path.join(SDF_DIR, f"{cid}_{safe_name}.sdf")
    if os.path.exists(sdf_path) and os.path.getsize(sdf_path) > 100:
        return sdf_path
    url_3d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    try:
        req = urllib.request.Request(url_3d, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=20) as resp, open(sdf_path, "wb") as f:
            f.write(resp.read())
        return sdf_path
    except Exception:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        w = Chem.SDWriter(sdf_path)
        w.write(mol)
        w.close()
        return sdf_path


def main():
    rows = []
    print("=" * 75)
    print("STEP 1 (PROTOCOL A): CURATION, ESOL, BOILED-EGG, PAINS & BRENK PROFILING")
    print("=" * 75)
    for idx, item in enumerate(COMPOUNDS, 1):
        name = item["name"]
        cid, smiles, formula = fetch_pubchem_data(name)
        sdf_path = ensure_3d_sdf(cid, name, smiles)
        desc = compute_admet_descriptors(smiles)

        if desc["ESOL_LogS_mol_L"] >= -4.0:
            sal_cls = "Soluble (Rapid Salivary Dissolution)"
        elif desc["ESOL_LogS_mol_L"] >= -6.0:
            sal_cls = "Moderately Soluble (Requires Co-Solvent/Surfactant)"
        else:
            sal_cls = "Poorly Soluble Brick-Dust (Requires Beta-CD Complexation)"

        row = {
            "Compound": name,
            "Category": item["category"],
            "Botanical_or_Chemical_Source": item["source"],
            "PubChem_CID": cid,
            "Molecular_Formula": formula,
            "SMILES": smiles,
            "SDF_File": os.path.basename(sdf_path),
            "Salivary_Biopharmaceutics_Class": sal_cls,
        }
        row.update(desc)
        rows.append(row)
        print(f"[{idx:02d}/22] {name:34s} | MW:{desc['MW_Da']:6.2f} | LogP:{desc['LogP_WildmanCrippen']:5.2f} | LogS:{desc['ESOL_LogS_mol_L']:5.2f} | Egg:{desc['BOILED_Egg_GI_Absorption'][:4]} | PAINS:{desc['PAINS_Structural_Alert'][:10]} | Brenk:{desc['Brenk_Structural_Alert'][:12]}")
        time.sleep(0.15)

    out_csv = os.path.join(RESULTS_DIR, "admet_and_physicochemical_properties.csv")
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"\n[SUCCESS] Saved updated 22-ligand ADMET + BOILED-Egg + Brenk table to: {out_csv}")


if __name__ == "__main__":
    main()
