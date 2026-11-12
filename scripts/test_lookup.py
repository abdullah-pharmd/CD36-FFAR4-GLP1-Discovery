import os
import time
import requests
import pubchempy as pcp

compounds = [
    # Positive Controls - FFAR4
    {"name": "TUG-891", "category": "Positive Control (FFAR4)", "target": "FFAR4"},
    {"name": "Oleic acid", "category": "Positive Control (FFAR4 / Natural agonist)", "target": "FFAR4"},
    # Positive Controls - CD36
    {"name": "Sulfo-N-succinimidyl oleate", "alt_name": "SSO", "category": "Positive Control (CD36)", "target": "CD36"},
    {"name": "Nitro-oleic acid", "alt_name": "9-Nitrooleic acid", "category": "Positive Control (CD36)", "target": "CD36"},
    # Test Ligands - Triterpenoids
    {"name": "Oleanolic acid", "category": "Pentacyclic Triterpenoid", "target": "CD36/FFAR4"},
    {"name": "Ursolic acid", "category": "Pentacyclic Triterpenoid", "target": "CD36/FFAR4"},
    {"name": "Corosolic acid", "category": "Pentacyclic Triterpenoid", "target": "CD36/FFAR4"},
    {"name": "Maslinic acid", "category": "Pentacyclic Triterpenoid", "target": "CD36/FFAR4"},
    {"name": "Betulinic acid", "category": "Pentacyclic Triterpenoid", "target": "CD36/FFAR4"},
    {"name": "Asiatic acid", "category": "Pentacyclic Triterpenoid", "target": "CD36/FFAR4"},
    {"name": "Glycyrrhetinic acid", "category": "Pentacyclic Triterpenoid", "target": "CD36/FFAR4"},
    {"name": "Gymnemic acid I", "alt_name": "Gymnemoside I", "category": "Triterpenoid Saponin", "target": "CD36/FFAR4"},
    {"name": "Gymnemagenin", "category": "Triterpenoid Sapogenin", "target": "CD36/FFAR4"},
    # Test Ligands - Alkamides / Sensory Lingual Modulators
    {"name": "Spilanthol", "alt_name": "Affinin", "category": "Alkamides / Lingual Sensory", "target": "CD36/FFAR4"},
    {"name": "Pellitorine", "category": "Alkamides / Lingual Sensory", "target": "CD36/FFAR4"},
    {"name": "Hydroxy-alpha-sanshool", "alt_name": "alpha-sanshool", "category": "Alkamides / Lingual Sensory", "target": "CD36/FFAR4"},
    {"name": "Piperine", "category": "Alkaloid / TRP modulator", "target": "CD36/FFAR4"},
    {"name": "Capsaicin", "category": "Capsaicinoid / TRP modulator", "target": "CD36/FFAR4"},
    # Test Ligands - Phenolic / Phenylpropanoids
    {"name": "6-Gingerol", "category": "Phenolic / Gingerol", "target": "CD36/FFAR4"},
    {"name": "6-Shogaol", "category": "Phenolic / Shogaol", "target": "CD36/FFAR4"},
    {"name": "Curcumin", "category": "Curcuminoid / Polyphenol", "target": "CD36/FFAR4"},
    {"name": "Berberine", "category": "Isoquinoline Alkaloid", "target": "CD36/FFAR4"},
    {"name": "Chlorogenic acid", "category": "Polyphenol / Hydroxycinnamic", "target": "CD36/FFAR4"},
    {"name": "Quercetin", "category": "Flavonoid / Flavonol", "target": "CD36/FFAR4"},
    {"name": "Resveratrol", "category": "Stilbenoid", "target": "CD36/FFAR4"},
    {"name": "Naringenin", "category": "Flavanone", "target": "CD36/FFAR4"},
    {"name": "Punicalagin", "category": "Ellagitannin", "target": "CD36/FFAR4"},
    {"name": "Fraxin", "category": "Coumarin glucoside", "target": "CD36/FFAR4"},
    {"name": "Phloretin", "category": "Dihydrochalcone", "target": "CD36/FFAR4"}
]

print(f"Total curated compounds: {len(compounds)}")
for c in compounds:
    query = c["name"]
    try:
        res = pcp.get_compounds(query, 'name')
        if not res and "alt_name" in c:
            res = pcp.get_compounds(c["alt_name"], 'name')
        if res:
            comp = res[0]
            print(f"[OK] {c['name']} -> CID: {comp.cid}, Formula: {comp.molecular_formula}, MW: {comp.molecular_weight}")
        else:
            print(f"[NOT FOUND] {c['name']}")
    except Exception as e:
        print(f"[ERR] {c['name']}: {e}")
    time.sleep(0.2)
