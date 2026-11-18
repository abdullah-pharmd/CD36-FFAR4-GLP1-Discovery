import os, sys
sys.path.insert(0, r"C:\Users\Abdullah\.gemini\antigravity\scratch\cd36_ffar4_glp1_study\scripts")
from structural_prep_utils import clean_pdb_strict, prepare_receptor_pdbqt

STRUCT_DIR = r"C:\Users\Abdullah\.gemini\antigravity\scratch\cd36_ffar4_glp1_study\structures"

# 1. 8G59 (Chain R, altLoc clean, salivary pH 6.8)
n8, c8 = clean_pdb_strict(os.path.join(STRUCT_DIR, "8G59.pdb"), os.path.join(STRUCT_DIR, "8G59_FFAR4_clean.pdb"), "R", "YN9", "401", os.path.join(STRUCT_DIR, "TUG891_crystal.pdb"))
prepare_receptor_pdbqt(os.path.join(STRUCT_DIR, "8G59_FFAR4_clean.pdb"), os.path.join(STRUCT_DIR, "8G59_FFAR4_clean.pdbqt"), ph=6.8)

# 2. 5LGD (Chain A, altLoc clean, pH 6.8 & pH 7.4)
n5, c5 = clean_pdb_strict(os.path.join(STRUCT_DIR, "5LGD.pdb"), os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdb"), "A", "PLM", "512", os.path.join(STRUCT_DIR, "PLM512_crystal.pdb"))
prepare_receptor_pdbqt(os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdb"), os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdbqt"), ph=6.8)
prepare_receptor_pdbqt(os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdb"), os.path.join(STRUCT_DIR, "5LGD_CD36_pH74.pdbqt"), ph=7.4)

# 3. 4A5S (Chain A, altLoc clean, plasma pH 7.4)
n4, c4 = clean_pdb_strict(os.path.join(STRUCT_DIR, "4A5S.pdb"), os.path.join(STRUCT_DIR, "4A5S_DPP4_clean.pdb"), "A", "N7F", "901", os.path.join(STRUCT_DIR, "N7F901_crystal.pdb"))
prepare_receptor_pdbqt(os.path.join(STRUCT_DIR, "4A5S_DPP4_clean.pdb"), os.path.join(STRUCT_DIR, "4A5S_DPP4_pH74.pdbqt"), ph=7.4)

for fn in ["8G59_FFAR4_clean.pdbqt", "5LGD_CD36_clean.pdbqt", "5LGD_CD36_pH74.pdbqt", "4A5S_DPP4_pH74.pdbqt"]:
    with open(os.path.join(STRUCT_DIR, fn)) as f:
        lines = [l for l in f if l.startswith("ATOM")]
    hd = sum(1 for l in lines if l[77:79].strip() == "HD")
    print(f"{fn:28s}: total ATOM={len(lines)}, polar HD={hd}")
