"""
02_generate_colab_notebook.py
Generates a clean, self-contained, checkpointed Google Colab Pro Jupyter Notebook
(ProtocolB_OpenMM_100ns_MD_MMGBSA.ipynb) for running explicit-water OpenMM MD
(AMBER ff14SB + TIP3P + OpenFF Sage 2.2.0 + AmberTools AM1-BCC) and native OBC2 MM-GBSA
across all 8 Protocol B hepatic CD36 (5LGD) and vascular DPP-4 (4A5S) complexes.

Defaults to 50 ns/complex (400 ns aggregate; 5 ns Google Drive checkpoints) — optimal for
Colab Pro NVIDIA L4 compute-unit budget (set PRODUCTION_NS=100.0 on A100 if desired).
Pins conda-forge cuda-version=12.6 (cuda-nvrtc 12.6.85, PTX 8.5) to prevent
CUDA_ERROR_UNSUPPORTED_PTX_VERSION on Colab host drivers.
"""
import os
import json

COLAB_DIR = os.path.dirname(os.path.abspath(__file__))
NB_PATH = os.path.join(COLAB_DIR, "ProtocolB_OpenMM_100ns_MD_MMGBSA.ipynb")


def md_cell(text):
    return {"cell_type": "markdown", "metadata": {}, "source": [line + "\n" for line in text.strip().splitlines()]}


def code_cell(code):
    lines = code.strip().splitlines()
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [line + "\n" for line in lines[:-1]] + ([lines[-1]] if lines else []),
    }


RUN_SCRIPT_BODY = r'''import os
import sys
# Ensure micromamba env bin (antechamber/sqm) is ALWAYS on PATH
_env_bin = os.path.dirname(os.path.abspath(sys.executable))
if _env_bin not in os.environ.get("PATH", ""):
    os.environ["PATH"] = _env_bin + ":" + os.environ.get("PATH", "")

import gc
import re
import json
import time
import glob
import zipfile
import warnings
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from rdkit import Chem
import mdtraj as md
import openmm as mm
import openmm.app as app
import openmm.unit as u
from pdbfixer import PDBFixer
from openff.toolkit.topology import Molecule as OffMolecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

warnings.filterwarnings("ignore")

# =============================================================================
# CONFIGURATION (Colab Pro L4 / A100 GPU)
# =============================================================================
INPUT_DIR = os.environ.get("INPUT_DIR", "/content/md_work/inputs")
DRIVE_OUT = os.environ.get("DRIVE_OUT", "/content/drive/MyDrive/CD36_DPP4_ProtocolB_MD/results")
LIVE_LOG_PATH = os.environ.get("LIVE_LOG_PATH", "/content/drive/MyDrive/CD36_DPP4_ProtocolB_MD/md_live.log")
os.makedirs(DRIVE_OUT, exist_ok=True)
if LIVE_LOG_PATH:
    os.makedirs(os.path.dirname(LIVE_LOG_PATH), exist_ok=True)


def log(msg):
    """Prints immediately and fsync-appends to Drive md_live.log so crashes never lose progress."""
    print(msg, flush=True)
    if LIVE_LOG_PATH:
        try:
            ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
            with open(LIVE_LOG_PATH, "a", encoding="utf-8") as lf:
                for line in str(msg).splitlines() or [""]:
                    lf.write(f"[{ts}] {line}\n")
                lf.flush()
                os.fsync(lf.fileno())
        except Exception:
            pass


def load_pdb_res_labels(rec_pdb):
    """Returns ordered list of crystallographic (resname, resseq) pairs from the input receptor PDB."""
    labels, seen = [], set()
    with open(rec_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM  "):
                key = (line[21], int(line[22:26]), line[17:20].strip())
                if key not in seen:
                    seen.add(key)
                    labels.append((key[2], key[1]))
    return labels


def fmt_res_label(res, pdb_res_labels):
    """Maps an MDTraj residue object cleanly to its crystallographic PDB resname+resseq."""
    if res.name == "LIG":
        return "LIG1"
    if 0 <= res.index < len(pdb_res_labels):
        rname, rseq = pdb_res_labels[res.index]
        return f"{rname}{rseq}"
    return f"{res.name}{res.resSeq}"


def _shift_res_token(tok, offset):
    m = re.match(r"^([A-Z]{3})(\d+)(:.*)?$", str(tok).strip())
    if not m or m.group(1) == "LIG":
        return str(tok)
    return f"{m.group(1)}{int(m.group(2)) + offset}{m.group(3) or ''}"


def repair_completed_residue_labels(name, cdir, rec_pdb):
    """Auto-repairs legacy 1-based residue labels (from keepIds=False) to exact PDB crystallographic IDs."""
    pdb_labels = load_pdb_res_labels(rec_pdb)
    if not pdb_labels:
        return
    offset = pdb_labels[0][1] - 1
    if offset == 0:
        return
    rmsf_p = os.path.join(cdir, "rmsf_per_residue.csv")
    if not os.path.exists(rmsf_p):
        return
    df_rmsf = pd.read_csv(rmsf_p)
    if len(df_rmsf) != len(pdb_labels):
        return
    first_lbl = str(df_rmsf["Residue"].iloc[0]).strip()
    expected_legacy = f"{pdb_labels[0][0]}1"
    if first_lbl != expected_legacy:
        return

    log(f"  [REPAIR] Remapping legacy 1-based residue labels (+{offset}) -> crystallographic PDB IDs for {name}...")
    df_rmsf["Residue"] = [f"{rn}{rs}" for rn, rs in pdb_labels]
    df_rmsf.to_csv(rmsf_p, index=False)

    hb_p = os.path.join(cdir, "hbond_occupancy.csv")
    top_hb_str = "Hydrophobic/vdW dominant"
    if os.path.exists(hb_p):
        hb_df = pd.read_csv(hb_p)
        if len(hb_df):
            for col in ("Residue", "Donor", "Acceptor"):
                if col in hb_df.columns:
                    hb_df[col] = hb_df[col].map(lambda x: _shift_res_token(x, offset))
            hb_df.to_csv(hb_p, index=False)
            top_hb_str = "; ".join(f"{r['Residue']} ({r['Occupancy_Pct']:.1f}%)" for _, r in hb_df.head(5).iterrows())

    wb_p = os.path.join(cdir, "water_bridges.csv")
    top_wb_str = "None >= 10%"
    if os.path.exists(wb_p):
        wb_df = pd.read_csv(wb_p)
        if len(wb_df):
            wb_df["Residue"] = wb_df["Residue"].map(lambda x: _shift_res_token(x, offset))
            wb_df.to_csv(wb_p, index=False)
            top_wb_str = "; ".join(f"{r['Residue']} ({r['Water_Bridge_Pct']:.1f}%)" for _, r in wb_df.head(4).iterrows())

    sum_p = os.path.join(cdir, "summary.json")
    if os.path.exists(sum_p):
        s_obj = json.load(open(sum_p))
        s_obj["Top_Direct_HBonds_Eq"] = top_hb_str
        s_obj["Top_Water_Bridges_Eq"] = top_wb_str
        with open(sum_p, "w") as fs:
            json.dump(s_obj, fs, indent=2)
    log(f"  [REPAIR OK] {name}: HB={top_hb_str} | WB={top_wb_str}")

PRODUCTION_NS = float(os.environ.get("PRODUCTION_NS", "50.0"))
CHECKPOINT_INTERVAL_NS = float(os.environ.get("CHECKPOINT_INTERVAL_NS", "5.0"))
SAVE_INTERVAL_PS = float(os.environ.get("SAVE_INTERVAL_PS", "100.0"))
MMGBSA_N_FRAMES = int(os.environ.get("MMGBSA_N_FRAMES", "250"))
EQ_SCALE = float(os.environ.get("EQ_SCALE", "1.0"))
ONLY_COMPLEX = os.environ.get("ONLY_COMPLEX", "").strip()

TIMESTEP_FS = 4.0                                    # 4 fs HMR (hydrogenMass = 3.0 amu; Hopkins et al., JCTC 2015)
HYDROGEN_MASS_AMU = 3.0
TEMPERATURE_K = 310.0
PRESSURE_BAR = 1.0
IONIC_STRENGTH_M = 0.15
PADDING_NM = 1.0

# Exact AMBER mbondi2 atomic radii (nm) and OBC2 (igb=5) scaling factors (Onufriev et al., Proteins 2004)
MBONDI2_PARAMS = {
    "H":  (0.120, 0.85),
    "C":  (0.170, 0.72),
    "N":  (0.155, 0.79),
    "O":  (0.150, 0.85),
    "F":  (0.150, 0.88),
    "P":  (0.185, 0.86),
    "S":  (0.180, 0.96),
    "Cl": (0.170, 0.80),
    "Br": (0.185, 0.80),
    "I":  (0.198, 0.80),
}


def get_platform():
    for pname in ("CUDA", "OpenCL", "CPU"):
        try:
            plat = mm.Platform.getPlatformByName(pname)
            props = {"Precision": "mixed"} if pname == "CUDA" else {}
            s_test = mm.System()
            s_test.addParticle(12.0)
            s_test.addParticle(12.0)
            nb_test = mm.NonbondedForce()
            nb_test.addParticle(0.0, 0.3, 0.1)
            nb_test.addParticle(0.0, 0.3, 0.1)
            s_test.addForce(nb_test)
            ctx_test = mm.Context(s_test, mm.VerletIntegrator(0.001), plat, props)
            ctx_test.setPositions([mm.Vec3(0, 0, 0), mm.Vec3(0.5, 0, 0)])
            _ = ctx_test.getState(getEnergy=True).getPotentialEnergy()
            del ctx_test, s_test
            print(f"[Platform] Verified OpenMM {pname} {props}", flush=True)
            return plat, props
        except Exception as exc:
            print(f"[Platform] Skipping {pname}: {exc}", flush=True)
            continue
    return mm.Platform.getPlatformByName("CPU"), {}


def load_parameterized_ligand(sdf_path, cache_dir):
    sup = Chem.SDMolSupplier(sdf_path, removeHs=False)
    rdmol = sup[0]
    assert rdmol is not None, f"Failed to read {sdf_path}"
    off_mol = OffMolecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    charge_cache = os.path.join(cache_dir, "am1bcc_charges.json")
    if os.path.exists(charge_cache):
        q_list = json.load(open(charge_cache))
        if len(q_list) == off_mol.n_atoms:
            off_mol.partial_charges = np.array(q_list, dtype=float) * u.elementary_charge
            print(f"    Loaded cached AM1-BCC charges from {charge_cache} ({len(q_list)} atoms).", flush=True)
            return off_mol
    t_q = time.time()
    print(f"    Assigning AM1-BCC partial charges (AmberTools sqm) for {os.path.basename(sdf_path)} (q={Chem.GetFormalCharge(rdmol):+d}, {off_mol.n_atoms} atoms)...", flush=True)
    off_mol.assign_partial_charges(partial_charge_method="am1bcc")
    q_vals = [float(q.m) for q in off_mol.partial_charges]
    with open(charge_cache, "w") as fq:
        json.dump(q_vals, fq)
    print(f"    AM1-BCC completed in {time.time() - t_q:.1f} s (cached to Drive).", flush=True)
    return off_mol


def add_obc2_gbsa_force(system, topology):
    """Attaches exact OpenMM native OBC2 (igb=5) polar Generalized Born solvation using mbondi2 radii."""
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    h_on_n = set()
    for b in topology.bonds():
        e1 = b.atom1.element.symbol if b.atom1.element else ""
        e2 = b.atom2.element.symbol if b.atom2.element else ""
        if e1 == "H" and e2 == "N":
            h_on_n.add(b.atom1.index)
        elif e2 == "H" and e1 == "N":
            h_on_n.add(b.atom2.index)

    gb = mm.GBSAOBCForce()
    gb.setNonbondedMethod(mm.GBSAOBCForce.NoCutoff)
    gb.setSoluteDielectric(1.0)
    gb.setSolventDielectric(78.5)
    gb.setSurfaceAreaEnergy(0.0)
    for atom in topology.atoms():
        q, _, _ = nb.getParticleParameters(atom.index)
        sym = atom.element.symbol if atom.element else "C"
        rad, scale = MBONDI2_PARAMS.get(sym, (0.170, 0.72))
        if sym == "H" and atom.index in h_on_n:
            rad = 0.130
        gb.addParticle(q, rad * u.nanometers, scale)
    system.addForce(gb)
    return system


def build_solvated_complex(rec_pdb, lig_sdf, out_dir):
    fixer = PDBFixer(filename=rec_pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    off_mol = load_parameterized_ligand(lig_sdf, out_dir)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=off_mol, forcefield="openff-2.2.0")
    ff = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
    ff.registerTemplateGenerator(smirnoff.generator)

    lig_top = off_mol.to_topology().to_openmm()
    for res in lig_top.residues():
        res.name = "LIG"
    lig_pos = off_mol.conformers[0].to_openmm()

    dry_mod = app.Modeller(fixer.topology, fixer.positions)
    dry_mod.add(lig_top, lig_pos)
    dry_top = dry_mod.topology
    rec_dry_top = fixer.topology
    lig_dry_top = lig_top

    with open(os.path.join(out_dir, "receptor_dry.pdb"), "w") as f:
        app.PDBFile.writeFile(rec_dry_top, fixer.positions, f, keepIds=False)
    with open(os.path.join(out_dir, "ligand_dry.pdb"), "w") as f:
        app.PDBFile.writeFile(lig_dry_top, lig_pos, f, keepIds=False)
    with open(os.path.join(out_dir, "complex_dry.pdb"), "w") as f:
        app.PDBFile.writeFile(dry_top, dry_mod.positions, f, keepIds=False)

    solv_mod = app.Modeller(fixer.topology, fixer.positions)
    solv_mod.add(lig_top, lig_pos)
    solv_mod.addSolvent(
        ff, model="tip3p", padding=PADDING_NM * u.nanometers,
        ionicStrength=IONIC_STRENGTH_M * u.molar, neutralize=True
    )
    solv_pdb = os.path.join(out_dir, "complex_solvated.pdb")
    with open(solv_pdb, "w") as f:
        app.PDBFile.writeFile(solv_mod.topology, solv_mod.positions, f, keepIds=False)

    return solv_mod.topology, solv_mod.positions, ff, off_mol, dry_top, rec_dry_top, lig_dry_top, solv_pdb


def run_complex_md(job, platform, plat_props):
    name = job["name"]
    cdir = os.path.join(DRIVE_OUT, name)
    os.makedirs(cdir, exist_ok=True)
    rec_pdb = os.path.join(INPUT_DIR, "receptors", job["receptor_pdb"])
    lig_sdf = os.path.join(INPUT_DIR, "ligands", job["sdf"])
    summary_json = os.path.join(cdir, "summary.json")
    if os.path.exists(summary_json):
        repair_completed_residue_labels(name, cdir, rec_pdb)
        log(f"[SKIP] {name} already completed on Google Drive ({summary_json})")
        return json.load(open(summary_json))

    t0 = time.time()
    log(f"\n{'='*76}\n[START] {name} | Target Production = {PRODUCTION_NS:.2f} ns (4 fs HMR)\n{'='*76}")

    chk_path = os.path.join(cdir, "prod_checkpoint.chk")
    state_xml_path = os.path.join(cdir, "prod_state.xml")
    sys_xml_path = os.path.join(cdir, "prod_system.xml")
    state_meta_path = os.path.join(cdir, "prod_progress.json")
    chunks_dir = os.path.join(cdir, "chunks")
    os.makedirs(chunks_dir, exist_ok=True)
    dcd_path = os.path.join(cdir, "prod_trajectory.dcd")
    log_path = os.path.join(cdir, "prod_state.csv")
    solv_pdb_existing = os.path.join(cdir, "complex_solvated.pdb")

    steps_per_ns = int(round((1e6) / TIMESTEP_FS))
    save_steps = max(1, int(round((SAVE_INTERVAL_PS * 1000.0) / TIMESTEP_FS)))
    total_prod_steps = max(save_steps * 4, int(round(PRODUCTION_NS * steps_per_ns)))
    chunk_steps = max(save_steps, int(round(CHECKPOINT_INTERVAL_NS * steps_per_ns)))

    existing_chunks = sorted(glob.glob(os.path.join(chunks_dir, "chunk_*.dcd")))
    prog_obj = json.load(open(state_meta_path)) if os.path.exists(state_meta_path) else None
    prog_steps = int(prog_obj["completed_steps"]) if prog_obj and "completed_steps" in prog_obj else 0
    valid_chunks = []
    for cp in existing_chunks:
        m_ck = re.search(r"chunk_(\d+)\.dcd$", cp)
        if m_ck and int(m_ck.group(1)) <= prog_steps and os.path.getsize(cp) > 1024:
            valid_chunks.append(cp)
    expected_n_chunks = int(round(prog_steps / chunk_steps)) if prog_steps > 0 else 0
    chunks_complete_to_prog = (expected_n_chunks > 0) and (len(valid_chunks) == expected_n_chunks)
    has_traj_on_drive = chunks_complete_to_prog or (os.path.exists(dcd_path) and os.path.getsize(dcd_path) > 1024)
    can_resume_xml = (
        os.path.exists(sys_xml_path)
        and prog_steps > 0
        and has_traj_on_drive
        and os.path.exists(solv_pdb_existing)
        and (os.path.exists(state_xml_path) or os.path.exists(chk_path))
    )

    if can_resume_xml:
        completed_steps = prog_steps
        for cp in existing_chunks:
            if cp not in valid_chunks:
                try:
                    os.remove(cp)
                except Exception:
                    pass
        log(f"  [RESUME-XML] Loading serialized system + {len(valid_chunks)} verified DCD chunk(s) at {completed_steps / steps_per_ns:.2f} / {PRODUCTION_NS:.2f} ns")
        fixer = PDBFixer(filename=rec_pdb)
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.4)
        off_mol = load_parameterized_ligand(lig_sdf, cdir)
        lig_top = off_mol.to_topology().to_openmm()
        for res in lig_top.residues():
            res.name = "LIG"
        dry_mod = app.Modeller(fixer.topology, fixer.positions)
        dry_mod.add(lig_top, off_mol.conformers[0].to_openmm())
        dry_top, rec_dry_top, lig_dry_top = dry_mod.topology, fixer.topology, lig_top
        solv_pdb = solv_pdb_existing
        top = app.PDBFile(solv_pdb).topology
        with open(sys_xml_path, "r") as fx:
            system = mm.XmlSerializer.deserialize(fx.read())
        integrator = mm.LangevinMiddleIntegrator(TEMPERATURE_K * u.kelvin, 1.0 / u.picosecond, TIMESTEP_FS * u.femtoseconds)
        sim = app.Simulation(top, system, integrator, platform, plat_props)
        if os.path.exists(state_xml_path):
            sim.loadState(state_xml_path)
        else:
            sim.loadCheckpoint(chk_path)
    else:
        # Clean any orphaned partial chunks when starting fresh from 0 ns
        for stale_p in existing_chunks:
            try:
                os.remove(stale_p)
            except Exception:
                pass
        top, pos, ff, off_mol, dry_top, rec_dry_top, lig_dry_top, solv_pdb = build_solvated_complex(rec_pdb, lig_sdf, cdir)

        smirnoff_chk = SMIRNOFFTemplateGenerator(molecules=off_mol, forcefield="openff-2.2.0")
        ff_vac_chk = app.ForceField("amber14/protein.ff14SB.xml")
        ff_vac_chk.registerTemplateGenerator(smirnoff_chk.generator)
        _sys_chk = ff_vac_chk.createSystem(dry_top, nonbondedMethod=app.NoCutoff, constraints=None)
        add_obc2_gbsa_force(_sys_chk, dry_top)
        del _sys_chk, ff_vac_chk, smirnoff_chk
        log(f"  [PRE-FLIGHT OK] Solvated ({top.getNumAtoms():,} atoms) & native OBC2 MM-GBSA verified.")

        heavy_idx = [
            a.index for a in top.atoms()
            if a.residue.name not in ("HOH", "NA", "CL", "Na+", "Cl-") and a.element is not None and a.element.symbol != "H"
        ]

        system = ff.createSystem(
            top, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers,
            constraints=app.HBonds, rigidWater=True, ewaldErrorTolerance=5e-4,
            hydrogenMass=HYDROGEN_MASS_AMU * u.amu
        )

        if os.path.exists(chk_path) and os.path.exists(state_meta_path) and os.path.exists(dcd_path) and os.path.getsize(dcd_path) > 1024:
            prog = json.load(open(state_meta_path))
            completed_steps = int(prog["completed_steps"])
            log(f"  [RESUME] Found legacy checkpoint + DCD at {completed_steps / steps_per_ns:.2f} / {PRODUCTION_NS:.2f} ns")
            system.addForce(mm.MonteCarloBarostat(PRESSURE_BAR * u.bar, TEMPERATURE_K * u.kelvin, 25))
            integrator = mm.LangevinMiddleIntegrator(TEMPERATURE_K * u.kelvin, 1.0 / u.picosecond, TIMESTEP_FS * u.femtoseconds)
            sim = app.Simulation(top, system, integrator, platform, plat_props)
            sim.loadCheckpoint(chk_path)
        else:
            completed_steps = 0
            restraint = mm.CustomExternalForce("0.5 * k * periodicdistance(x, y, z, x0, y0, z0)^2")
            restraint.addGlobalParameter("k", 1000.0 * u.kilojoules_per_mole / u.nanometer**2)
            restraint.addPerParticleParameter("x0")
            restraint.addPerParticleParameter("y0")
            restraint.addPerParticleParameter("z0")
            for idx in heavy_idx:
                restraint.addParticle(idx, pos[idx].value_in_unit(u.nanometers))
            rest_idx = system.addForce(restraint)

            integrator_nvt = mm.LangevinMiddleIntegrator(50.0 * u.kelvin, 1.0 / u.picosecond, 2.0 * u.femtoseconds)
            sim_nvt = app.Simulation(top, system, integrator_nvt, platform, plat_props)
            sim_nvt.context.setPositions(pos)
            min_iters = max(250, int(round(2500 * EQ_SCALE)))
            nvt_sub = max(200, int(round(10000 * EQ_SCALE)))
            npt_r = max(500, int(round(250000 * EQ_SCALE)))
            npt_u = max(500, int(round(150000 * EQ_SCALE)))
            log(f"  [EQ 1/4] Minimizing energy ({min_iters} steps)...")
            sim_nvt.minimizeEnergy(maxIterations=min_iters)
            sim_nvt.context.setVelocitiesToTemperature(50.0 * u.kelvin, 42)
            log(f"  [EQ 2/4] NVT heating 50 K -> 310 K ({10 * nvt_sub * 2 / 1000:.1f} ps, k=1000)...")
            for temp_t in np.linspace(50.0, TEMPERATURE_K, 10):
                integrator_nvt.setTemperature(float(temp_t) * u.kelvin)
                sim_nvt.step(nvt_sub)
            nvt_state = sim_nvt.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
            del sim_nvt
            gc.collect()

            system.addForce(mm.MonteCarloBarostat(PRESSURE_BAR * u.bar, TEMPERATURE_K * u.kelvin, 25))
            integrator_npt = mm.LangevinMiddleIntegrator(TEMPERATURE_K * u.kelvin, 1.0 / u.picosecond, 2.0 * u.femtoseconds)
            sim_npt = app.Simulation(top, system, integrator_npt, platform, plat_props)
            sim_npt.context.setPeriodicBoxVectors(*nvt_state.getPeriodicBoxVectors())
            sim_npt.context.setPositions(nvt_state.getPositions())
            sim_npt.context.setVelocities(nvt_state.getVelocities())
            log(f"  [EQ 3/4] Restrained NPT equilibration ({npt_r * 2 / 1000:.1f} ps, k=100)...")
            sim_npt.context.setParameter("k", 100.0 * u.kilojoules_per_mole / u.nanometer**2)
            sim_npt.step(npt_r)
            log(f"  [EQ 4/4] Unrestrained NPT equilibration ({npt_u * 2 / 1000:.1f} ps, k=0)...")
            sim_npt.context.setParameter("k", 0.0)
            sim_npt.step(npt_u)

            eq_state = sim_npt.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
            eq_pdb = os.path.join(cdir, "equilibrated_solvated.pdb")
            with open(eq_pdb, "w") as f:
                app.PDBFile.writeFile(top, eq_state.getPositions(), f, keepIds=False)
            del sim_npt
            gc.collect()

            system.removeForce(rest_idx)
            with open(sys_xml_path, "w") as fx:
                fx.write(mm.XmlSerializer.serialize(system))
            integrator = mm.LangevinMiddleIntegrator(TEMPERATURE_K * u.kelvin, 1.0 / u.picosecond, TIMESTEP_FS * u.femtoseconds)
            sim = app.Simulation(top, system, integrator, platform, plat_props)
            sim.context.setPeriodicBoxVectors(*eq_state.getPeriodicBoxVectors())
            sim.context.setPositions(eq_state.getPositions())
            sim.context.setVelocities(eq_state.getVelocities())

    # Write closed, FUSE-synced 5 ns DCD + state CSV chunks so disconnects NEVER lose completed frames!
    while completed_steps < total_prod_steps:
        step_now = min(chunk_steps, total_prod_steps - completed_steps)
        target_steps = completed_steps + step_now
        chunk_tag = f"{target_steps:010d}"
        chunk_dcd = os.path.join(chunks_dir, f"chunk_{chunk_tag}.dcd")
        chunk_csv = os.path.join(chunks_dir, f"state_{chunk_tag}.csv")

        sim.reporters.clear()
        dcd_rep = app.DCDReporter(chunk_dcd, save_steps)
        csv_rep = app.StateDataReporter(
            chunk_csv, save_steps, step=True, time=True, potentialEnergy=True,
            temperature=True, volume=True, speed=True
        )
        sim.reporters.append(dcd_rep)
        sim.reporters.append(csv_rep)

        t_chunk = time.time()
        sim.step(step_now)
        completed_steps = target_steps

        # Explicitly close both chunk reporters so Google Drive FUSE immediately uploads chunk_*.dcd!
        sim.reporters.clear()
        for _rep in (dcd_rep, csv_rep):
            try:
                if hasattr(_rep, "_out") and _rep._out is not None and not _rep._out.closed:
                    _rep._out.flush()
                    os.fsync(_rep._out.fileno())
                    _rep._out.close()
            except Exception:
                pass
        del dcd_rep, csv_rep
        gc.collect()

        sim.saveCheckpoint(chk_path)
        sim.saveState(state_xml_path)
        with open(state_meta_path, "w") as fp:
            json.dump({"completed_steps": completed_steps, "completed_ns": completed_steps / steps_per_ns}, fp)
            fp.flush()
            os.fsync(fp.fileno())
        ns_done = completed_steps / steps_per_ns
        dt_min = (time.time() - t_chunk) / 60.0
        ns_per_day = (step_now / steps_per_ns) / max(dt_min / 1440.0, 1e-9)
        mb_sz = os.path.getsize(chunk_dcd) / (1024 * 1024) if os.path.exists(chunk_dcd) else 0.0
        log(f"  [PROD] {name}: {ns_done:6.2f} / {PRODUCTION_NS:.2f} ns ({dt_min:.2f} min/chunk | {ns_per_day:.1f} ns/day | synced {mb_sz:.1f} MB DCD)")

    fin_state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
    with open(os.path.join(cdir, "final_solvated.pdb"), "w") as f:
        app.PDBFile.writeFile(top, fin_state.getPositions(), f, keepIds=False)

    del sim, system
    gc.collect()

    # Merge chunk CSVs into prod_state.csv cleanly for archival parity
    chunk_csvs = sorted(glob.glob(os.path.join(chunks_dir, "state_*.csv")))
    if chunk_csvs:
        with open(log_path, "w", encoding="utf-8") as fout:
            for ci, cp in enumerate(chunk_csvs):
                lines = open(cp, "r", encoding="utf-8").readlines()
                fout.writelines(lines if ci == 0 else lines[1:])

    log(f"  [ANALYSIS] Computing RMSD, RMSF, Rg, H-bonds, Water Bridges, & OBC2 MM-GBSA for {name}...")
    pdb_res_labels = load_pdb_res_labels(rec_pdb)
    summary = analyze_and_mmgbsa(name, cdir, solv_pdb, dry_top, rec_dry_top, lig_dry_top, dcd_path, off_mol, platform, plat_props, pdb_res_labels)
    summary["Elapsed_Minutes"] = round((time.time() - t0) / 60.0, 2)
    with open(summary_json, "w") as fs:
        json.dump(summary, fs, indent=2)
    log(f"  [DONE] {name} | Prot RMSD={summary['Prot_CA_RMSD_Mean_A']:.2f}±{summary['Prot_CA_RMSD_SD_A']:.2f} Å | "
        f"Lig RMSD={summary['Lig_Heavy_RMSD_Mean_A']:.2f}±{summary['Lig_Heavy_RMSD_SD_A']:.2f} Å | "
        f"MM-GBSA={summary['MMGBSA_dG_Mean_kcal_mol']:.2f}±{summary['MMGBSA_dG_SEM_kcal_mol']:.2f} kcal/mol")
    return summary


def compute_implicit_energies(omm_topology, traj_xyz_nm, ff_vac, platform, plat_props):
    sys_gb = ff_vac.createSystem(omm_topology, nonbondedMethod=app.NoCutoff, constraints=None)
    add_obc2_gbsa_force(sys_gb, omm_topology)
    integ = mm.VerletIntegrator(0.001 * u.picoseconds)
    ctx = mm.Context(sys_gb, integ, platform, plat_props)
    energies = []
    for xyz in traj_xyz_nm:
        ctx.setPositions(xyz * u.nanometers)
        e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)
        energies.append(float(e))
    del ctx, sys_gb
    return np.array(energies)


def analyze_and_mmgbsa(name, cdir, solv_pdb, dry_top, rec_dry_top, lig_dry_top, dcd_path, off_mol, platform, plat_props, pdb_res_labels=None):
    chunk_dcds = sorted(glob.glob(os.path.join(cdir, "chunks", "chunk_*.dcd")))
    if chunk_dcds:
        trajs = []
        if os.path.exists(dcd_path) and os.path.getsize(dcd_path) > 1024:
            trajs.append(md.load(dcd_path, top=solv_pdb))
        for cdcd in chunk_dcds:
            trajs.append(md.load(cdcd, top=solv_pdb))
        traj = md.join(trajs)
        del trajs
    else:
        traj = md.load(dcd_path, top=solv_pdb)
    top_md = traj.topology
    n_rec = rec_dry_top.getNumAtoms()
    n_dry = dry_top.getNumAtoms()
    prot_idx = np.arange(n_rec, dtype=int)
    lig_all_idx = np.arange(n_rec, n_dry, dtype=int)
    dry_idx = np.arange(n_dry, dtype=int)
    traj.image_molecules(inplace=True, anchor_molecules=[set(top_md.atom(i) for i in prot_idx)])

    ca_idx = top_md.select("protein and name CA")
    lig_heavy_idx = np.array([i for i in lig_all_idx if top_md.atom(i).element is not None and top_md.atom(i).element.symbol != "H"], dtype=int)

    traj.superpose(traj, frame=0, atom_indices=ca_idx)

    prot_rmsd_A = md.rmsd(traj, traj, frame=0, atom_indices=ca_idx) * 10.0
    lig_xyz = traj.xyz[:, lig_heavy_idx, :]
    lig_rmsd_A = np.sqrt(np.mean(np.sum((lig_xyz - lig_xyz[0:1]) ** 2, axis=2), axis=1)) * 10.0

    n_frames = traj.n_frames
    eq_start = max(1, n_frames // 2)
    rg_A = md.compute_rg(traj.atom_slice(prot_idx)) * 10.0
    rmsf_A = md.rmsf(traj[eq_start:], traj[eq_start:], frame=0, atom_indices=ca_idx) * 10.0

    time_ns = np.linspace(SAVE_INTERVAL_PS / 1000.0, n_frames * SAVE_INTERVAL_PS / 1000.0, n_frames)
    pd.DataFrame({
        "Time_ns": np.round(time_ns, 3),
        "Prot_CA_RMSD_A": np.round(prot_rmsd_A, 3),
        "Lig_Heavy_RMSD_A": np.round(lig_rmsd_A, 3),
        "Prot_Rg_A": np.round(rg_A, 3),
    }).to_csv(os.path.join(cdir, "rmsd_rg_timeseries.csv"), index=False)

    pdb_res_labels = pdb_res_labels or []
    ca_res = [fmt_res_label(top_md.atom(i).residue, pdb_res_labels) for i in ca_idx]
    pd.DataFrame({"Residue": ca_res, "RMSF_A": np.round(rmsf_A, 3)}).to_csv(os.path.join(cdir, "rmsf_per_residue.csv"), index=False)

    traj_dry = traj.atom_slice(dry_idx)
    traj_dry[-1].save_pdb(os.path.join(cdir, "final_dry_aligned.pdb"))
    hbonds = md.baker_hubbard(traj_dry[eq_start:], freq=0.05, periodic=False)
    md_dry_top = traj_dry.topology
    hb_rows = []
    for d_i, h_i, a_i in hbonds:
        r_d = md_dry_top.atom(d_i).residue
        r_a = md_dry_top.atom(a_i).residue
        if (r_d.name == "LIG") != (r_a.name == "LIG"):
            prot_res = r_a if r_d.name == "LIG" else r_d
            dists = md.compute_distances(traj_dry[eq_start:], [[d_i, a_i]])[:, 0]
            angs = md.compute_angles(traj_dry[eq_start:], [[d_i, h_i, a_i]])[:, 0]
            occ = float(np.mean((dists <= 0.35) & (angs >= (2.0 * np.pi / 3.0))) * 100.0)
            if occ >= 5.0:
                hb_rows.append({
                    "Residue": fmt_res_label(prot_res, pdb_res_labels),
                    "Donor": f"{fmt_res_label(r_d, pdb_res_labels)}:{md_dry_top.atom(d_i).name}",
                    "Acceptor": f"{fmt_res_label(r_a, pdb_res_labels)}:{md_dry_top.atom(a_i).name}",
                    "Occupancy_Pct": round(occ, 1),
                    "Mean_Dist_A": round(float(np.mean(dists[dists <= 0.35]) * 10.0) if np.any(dists <= 0.35) else float(np.mean(dists) * 10.0), 2),
                })
    hb_df = pd.DataFrame(hb_rows).sort_values("Occupancy_Pct", ascending=False) if hb_rows else pd.DataFrame(columns=["Residue", "Donor", "Acceptor", "Occupancy_Pct", "Mean_Dist_A"])
    hb_df.to_csv(os.path.join(cdir, "hbond_occupancy.csv"), index=False)

    water_o_idx = top_md.select("water and name O")
    lig_pol_idx = top_md.select("resname LIG and (element N or element O or element S)")
    prot_pol_idx = top_md.select("protein and (element N or element O or element S)")
    bridge_counts = {}
    sub_frames = traj[eq_start::max(1, 2 if (n_frames - eq_start) >= 20 else 1)]
    if len(water_o_idx) > 0 and len(lig_pol_idx) > 0:
        for fr in sub_frames:
            w_near_lig = md.compute_neighbors(fr, 0.35, query_indices=lig_pol_idx, haystack_indices=water_o_idx)[0]
            if len(w_near_lig) == 0:
                continue
            p_near_w = md.compute_neighbors(fr, 0.35, query_indices=w_near_lig, haystack_indices=prot_pol_idx)[0]
            seen_res_frame = set()
            for pi in p_near_w:
                r = top_md.atom(pi).residue
                seen_res_frame.add(fmt_res_label(r, pdb_res_labels))
            for rk in seen_res_frame:
                bridge_counts[rk] = bridge_counts.get(rk, 0) + 1
    wb_rows = [
        {"Residue": k, "Water_Bridge_Pct": round(100.0 * v / len(sub_frames), 1)}
        for k, v in sorted(bridge_counts.items(), key=lambda x: -x[1])
        if (100.0 * v / len(sub_frames)) >= 10.0
    ]
    pd.DataFrame(wb_rows).to_csv(os.path.join(cdir, "water_bridges.csv"), index=False)

    sample_idx = np.linspace(eq_start, n_frames - 1, min(MMGBSA_N_FRAMES, n_frames - eq_start), dtype=int)
    traj_sub = traj_dry[sample_idx]
    rec_only_idx = np.arange(n_rec, dtype=int)
    lig_only_idx = np.arange(n_rec, n_dry, dtype=int)

    smirnoff = SMIRNOFFTemplateGenerator(molecules=off_mol, forcefield="openff-2.2.0")
    ff_vac_com = app.ForceField("amber14/protein.ff14SB.xml")
    ff_vac_com.registerTemplateGenerator(smirnoff.generator)
    ff_vac_rec = app.ForceField("amber14/protein.ff14SB.xml")
    ff_vac_lig = app.ForceField()
    ff_vac_lig.registerTemplateGenerator(smirnoff.generator)

    e_com = compute_implicit_energies(dry_top, traj_sub.xyz, ff_vac_com, platform, plat_props)
    e_rec = compute_implicit_energies(rec_dry_top, traj_sub.atom_slice(rec_only_idx).xyz, ff_vac_rec, platform, plat_props)
    e_lig = compute_implicit_energies(lig_dry_top, traj_sub.atom_slice(lig_only_idx).xyz, ff_vac_lig, platform, plat_props)

    sasa_com = np.sum(md.shrake_rupley(traj_sub, mode="residue"), axis=1) * 100.0
    sasa_rec = np.sum(md.shrake_rupley(traj_sub.atom_slice(rec_only_idx), mode="residue"), axis=1) * 100.0
    sasa_lig = np.sum(md.shrake_rupley(traj_sub.atom_slice(lig_only_idx), mode="residue"), axis=1) * 100.0
    dg_sa = 0.00542 * (sasa_com - sasa_rec - sasa_lig)
    dg_tot = (e_com - e_rec - e_lig) + dg_sa

    pd.DataFrame({
        "Frame": sample_idx,
        "Time_ns": np.round(time_ns[sample_idx], 2),
        "E_Com_kcal": np.round(e_com, 3),
        "E_Rec_kcal": np.round(e_rec, 3),
        "E_Lig_kcal": np.round(e_lig, 3),
        "dG_SA_kcal": np.round(dg_sa, 3),
        "dG_MMGBSA_kcal": np.round(dg_tot, 3),
    }).to_csv(os.path.join(cdir, "mmgbsa_frames.csv"), index=False)

    top_hb_str = "; ".join(f"{r['Residue']} ({r['Occupancy_Pct']:.1f}%)" for _, r in hb_df.head(5).iterrows()) if len(hb_df) else "Hydrophobic/vdW dominant"
    top_wb_str = "; ".join(f"{r['Residue']} ({r['Water_Bridge_Pct']:.1f}%)" for r in wb_rows[:4]) if wb_rows else "None >= 10%"

    return {
        "Complex": name,
        "Production_ns": PRODUCTION_NS,
        "N_Frames": int(n_frames),
        "Prot_CA_RMSD_Mean_A": round(float(np.mean(prot_rmsd_A[eq_start:])), 3),
        "Prot_CA_RMSD_SD_A": round(float(np.std(prot_rmsd_A[eq_start:])), 3),
        "Lig_Heavy_RMSD_Mean_A": round(float(np.mean(lig_rmsd_A[eq_start:])), 3),
        "Lig_Heavy_RMSD_SD_A": round(float(np.std(lig_rmsd_A[eq_start:])), 3),
        "Prot_Rg_Mean_A": round(float(np.mean(rg_A[eq_start:])), 3),
        "Prot_Rg_SD_A": round(float(np.std(rg_A[eq_start:])), 3),
        "MMGBSA_dG_Mean_kcal_mol": round(float(np.mean(dg_tot)), 2),
        "MMGBSA_dG_SD_kcal_mol": round(float(np.std(dg_tot)), 2),
        "MMGBSA_dG_SEM_kcal_mol": round(float(np.std(dg_tot) / np.sqrt(len(dg_tot))), 2),
        "MMGBSA_dG_SA_Mean_kcal_mol": round(float(np.mean(dg_sa)), 2),
        "Top_Direct_HBonds_Eq": top_hb_str,
        "Top_Water_Bridges_Eq": top_wb_str,
    }


def main():
    full_manifest = json.load(open(os.path.join(INPUT_DIR, "md_manifest.json")))
    if ONLY_COMPLEX:
        tokens = [t.strip() for t in ONLY_COMPLEX.split(",") if t.strip()]
        manifest = [j for j in full_manifest if any(j["name"] == t or j["name"].startswith(t) for t in tokens)]
        log(f"[FILTER] ONLY_COMPLEX={ONLY_COMPLEX!r} -> running {len(manifest)}/{len(full_manifest)} complex(es): {[j['name'] for j in manifest]}")
    else:
        manifest = full_manifest

    platform, plat_props = get_platform()
    summaries = []
    for job in manifest:
        summaries.append(run_complex_md(job, platform, plat_props))

    # Always check whether all 8 complexes in full_manifest now exist on Drive
    all_summaries = []
    for job in full_manifest:
        sp = os.path.join(DRIVE_OUT, job["name"], "summary.json")
        if os.path.exists(sp):
            rec_pdb = os.path.join(INPUT_DIR, "receptors", job["receptor_pdb"])
            repair_completed_residue_labels(job["name"], os.path.join(DRIVE_OUT, job["name"]), rec_pdb)
            all_summaries.append(json.load(open(sp)))

    df_run = pd.DataFrame(summaries)
    log("\n" + "=" * 76)
    log(f"SESSION COMPLETED {len(summaries)} COMPLEX(ES) ({len(all_summaries)}/{len(full_manifest)} TOTAL ON DRIVE):")
    log(df_run[["Complex", "Prot_CA_RMSD_Mean_A", "Lig_Heavy_RMSD_Mean_A", "MMGBSA_dG_Mean_kcal_mol", "MMGBSA_dG_SEM_kcal_mol"]].to_string(index=False))

    if len(all_summaries) == len(full_manifest):
        df_all = pd.DataFrame(all_summaries)
        sum_csv = os.path.join(DRIVE_OUT, "md_summary_all_complexes.csv")
        df_all.to_csv(sum_csv, index=False)
        bundle_zip = os.path.join(DRIVE_OUT, "protocol_b_md_results_bundle.zip")
        with zipfile.ZipFile(bundle_zip, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.write(sum_csv, "md_summary_all_complexes.csv")
            for job in full_manifest:
                nm = job["name"]
                for fn in ("summary.json", "rmsd_rg_timeseries.csv", "rmsf_per_residue.csv", "hbond_occupancy.csv", "water_bridges.csv", "mmgbsa_frames.csv", "final_dry_aligned.pdb"):
                    fp = os.path.join(DRIVE_OUT, nm, fn)
                    if os.path.exists(fp):
                        zf.write(fp, f"{nm}/{fn}")
        log(f"\n[BUNDLE READY] ALL {len(full_manifest)} COMPLEXES COMPLETE -> {bundle_zip} ({os.path.getsize(bundle_zip)/1024:.1f} KB)")
    else:
        rem = [j["name"] for j in full_manifest if not os.path.exists(os.path.join(DRIVE_OUT, j["name"], "summary.json"))]
        log(f"[STATUS] {len(all_summaries)}/{len(full_manifest)} complete. Remaining: {rem}")


if __name__ == "__main__":
    main()
'''


def main():
    cells = [
        md_cell(
            "# Protocol B — 50 ns/Complex (400 ns Aggregate) Explicit-Water OpenMM MD & OBC2 MM-GBSA Suite\n"
            "**Targets:** Human DPP-4 (`PDB: 4A5S`, 728 residues, 0 gaps) & Hepatic CD36 Ectodomain (`PDB: 5LGD`, 400 residues, 0 gaps) at physiological pH 7.4.\n\n"
            "**Force Fields:** `AMBER ff14SB` + `TIP3P` water (`0.15 M NaCl`) + `OpenFF Sage 2.2.0` (`AM1-BCC` partial charges via AmberTools `sqm`, cached to Drive), `4 fs` HMR (`3.0 amu`, Hopkins et al. *JCTC* 2015), `310 K`, `1 bar`, Google Drive auto-checkpointing every `5 ns`, native OpenMM `OBC2 (igb=5, mbondi2)` MM-GBSA (`250` equilibrated frames).\n\n"
            "**8 Complexes:**\n"
            "- **DPP-4 (`4A5S`):** `N7F-901` (crystal control), `Teneligliptin`, `Naringenin`, `Berberine`\n"
            "- **CD36 (`5LGD`):** `PLM-512` (native palmitate), `Linagliptin`, `Teneligliptin`, `Naringenin`\n\n"
            "> **Colab Pro Tip:** Select **Runtime → Change runtime type → L4 GPU** (50 ns default fits 100 CU budget cleanly) or **A100 GPU**. If disconnected, simply re-run Cells 1–4—each complex automatically resumes from its last 5 ns Google Drive checkpoint."
        ),
        code_cell(
            "# Cell 1: Verify GPU, Mount Google Drive, and Install Isolated Micromamba OpenMM (CUDA 12.6 PTX-Safe!) + OpenFF + AmberTools\n"
            "import os, subprocess\n"
            "from google.colab import drive\n"
            "drive.mount('/content/drive')\n"
            "!nvidia-smi\n\n"
            "MM_BIN = '/content/micromamba/bin/micromamba'\n"
            "ENV_PY = '/content/micromamba/envs/md/bin/python'\n"
            "if not os.path.exists(ENV_PY):\n"
            "    !mkdir -p /content/micromamba/bin\n"
            "    !curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /content/micromamba/bin --strip-components=1 bin/micromamba\n"
            "    !{MM_BIN} create -y -r /content/micromamba -n md -c conda-forge \\\n"
            "        python=3.11 'openmm>=8.1' 'cuda-version=12.6' pdbfixer openff-toolkit openmmforcefields ambertools mdtraj rdkit pandas matplotlib scipy\n\n"
            "with open('/content/test_cuda_ptx.py', 'w') as _f:\n"
            "    _f.write(\n"
            "        'import openmm as mm\\n'\n"
            "        's = mm.System(); s.addParticle(12.0); s.addParticle(12.0)\\n'\n"
            "        'nb = mm.NonbondedForce(); nb.addParticle(0, 0.3, 0.1); nb.addParticle(0, 0.3, 0.1); s.addForce(nb)\\n'\n"
            "        'p = mm.Platform.getPlatformByName(\"CUDA\")\\n'\n"
            "        'c = mm.Context(s, mm.VerletIntegrator(0.001), p, {\"Precision\": \"mixed\"})\\n'\n"
            "        'c.setPositions([mm.Vec3(0, 0, 0), mm.Vec3(0.5, 0, 0)])\\n'\n"
            "        'print(\"[CUDA PTX SELF-TEST PASSED] OpenMM\", mm.__version__, c.getState(getEnergy=True).getPotentialEnergy())\\n'\n"
            "    )\n"
            "test_cmd = 'unset OPENMM_CUDA_COMPILER; LD_LIBRARY_PATH=/content/micromamba/envs/md/lib:$LD_LIBRARY_PATH /content/micromamba/envs/md/bin/python /content/test_cuda_ptx.py'\n"
            "res = subprocess.run(test_cmd, shell=True, capture_output=True, text=True)\n"
            "print(res.stdout or res.stderr)\n"
            "if res.returncode != 0:\n"
            "    print('[REPAIR] Re-pinning cuda-version=12.6 + cuda-nvrtc=12.6 in existing md env...')\n"
            "    !rm -rf ~/.openmm\n"
            "    !{MM_BIN} install -y -r /content/micromamba -n md -c conda-forge 'cuda-version=12.6' 'cuda-nvrtc=12.6.*' 'openmm>=8.1'\n"
            "    subprocess.run(test_cmd, shell=True, check=True)\n"
        ),
        code_cell(
            "# Cell 2: Load protocol_b_md_inputs.zip (from Google Drive or direct upload)\n"
            "import os, zipfile, shutil\n"
            "from google.colab import files\n\n"
            "DRIVE_IN_ZIP = '/content/drive/MyDrive/CD36_DPP4_ProtocolB_MD/protocol_b_md_inputs.zip'\n"
            "os.makedirs('/content/drive/MyDrive/CD36_DPP4_ProtocolB_MD', exist_ok=True)\n"
            "os.makedirs('/content/md_work/inputs', exist_ok=True)\n\n"
            "if not os.path.exists(DRIVE_IN_ZIP):\n"
            "    print('Upload protocol_b_md_inputs.zip from protocol_b_hepatic_cd36_dpp4/colab_md/ ...')\n"
            "    uploaded = files.upload()\n"
            "    fname = list(uploaded.keys())[0]\n"
            "    shutil.move(fname, DRIVE_IN_ZIP)\n\n"
            "with zipfile.ZipFile(DRIVE_IN_ZIP, 'r') as zf:\n"
            "    zf.extractall('/content/md_work/inputs')\n"
            "!ls -la /content/md_work/inputs /content/md_work/inputs/receptors /content/md_work/inputs/ligands"
        ),
        code_cell("%%writefile /content/run_md_suite.py\n" + RUN_SCRIPT_BODY),
        code_cell(
            "# Cell 4: Run / Resume 8-Complex 50 ns Explicit-Water MD + OBC2 MM-GBSA Suite\n"
            "# Optional parallel speedup on Colab Pro (2 L4 tabs): set ONLY_COMPLEX='4A5S' in Tab 1 & ONLY_COMPLEX='5LGD' in Tab 2!\n"
            "!unset OPENMM_CUDA_COMPILER && PATH=\"/content/micromamba/envs/md/bin:$PATH\" LD_LIBRARY_PATH=\"/content/micromamba/envs/md/lib:$LD_LIBRARY_PATH\" MPLBACKEND=Agg \\\n"
            "    PRODUCTION_NS=50.0 CHECKPOINT_INTERVAL_NS=5.0 SAVE_INTERVAL_PS=100.0 MMGBSA_N_FRAMES=250 ONLY_COMPLEX=\"\" \\\n"
            "    /content/micromamba/envs/md/bin/python -u /content/run_md_suite.py"
        ),
        code_cell(
            "# Cell 5: Inspect Summary Table & Download protocol_b_md_results_bundle.zip back to Local PC\n"
            "import pandas as pd\n"
            "from google.colab import files\n"
            "sum_csv = '/content/drive/MyDrive/CD36_DPP4_ProtocolB_MD/results/md_summary_all_complexes.csv'\n"
            "display(pd.read_csv(sum_csv))\n"
            "files.download('/content/drive/MyDrive/CD36_DPP4_ProtocolB_MD/results/protocol_b_md_results_bundle.zip')"
        ),
    ]

    nb = {
        "cells": cells,
        "metadata": {
            "accelerator": "GPU",
            "colab": {"gpuType": "L4", "provenance": []},
            "kernelspec": {"display_name": "Python 3", "name": "python3"},
            "language_info": {"name": "python", "version": "3.11"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    with open(NB_PATH, "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=2)
    print(f"[SUCCESS] Generated Colab Pro Notebook: {NB_PATH} ({os.path.getsize(NB_PATH)/1024:.1f} KB)")


if __name__ == "__main__":
    main()
