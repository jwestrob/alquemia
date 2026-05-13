#!/usr/bin/env python3
"""
Phase 2 PRE-CHECK: alchemical BVS on 8DQ2 X-ray (HansLanM).

Goal:
  Validate that the OpenMM-Merz-12-6-4 pipeline produces non-trivial
  Ln-discriminating geometry on a known-good X-ray structure (PDB 8DQ2,
  HansLanM, 3 La sites in chain A) BEFORE we scale to a Protenix fold of
  canonical M. extorquens LanM.

Pre-check pipeline per (site, target_Ln):
  1. Strip waters, keep chain A protein + 3 La ions.
  2. Add hydrogens (PDBFixer at pH 7).
  3. Build OpenMM system: amber14-all + amber14/tip3p (for the La residue
     template). VACUUM (no implicit solvent — see DESIGN NOTES below).
  4. Override target La's nonbonded params with target Ln's Merz 12-6-4
     values from fep/prep/ln_parameters.json (Rmin/2, eps, charge=+3).
  5. Add CustomNonbondedForce for the C4 ion-induced-dipole correction.
  6. Restraints (k = 50 kcal/mol/Å² = 20920 kJ/mol/nm²):
       - All protein backbone heavy atoms (N, CA, C, O): restrained.
       - Sidechain heavy atoms of residues NOT within 6 Å of target metal:
         restrained.
       - Sidechain heavy atoms of residues WITHIN 6 Å of target metal: free.
       - Other 2 La ions: restrained at X-ray positions.
  7. L-BFGS minimize: tolerance = 0.4 kJ/mol/nm (~0.01 kcal/mol/Å),
     maxIters = 5000.
  8. Compute BVS at relaxed geometry using on_scanner.bvs.compute_bvs()
     for the target Ln cation.
  9. Record: BVS, residual=|BVS-q|, n_donors, mean Ln-O dist, max Ln-O dist,
     RMSD of free atoms vs X-ray.

Output:
  alchemical_bvs/results/precheck_8dq2_residuals.tsv  — per (site, Ln) rows
  alchemical_bvs/results/precheck_8dq2_donors.tsv     — per (site, Ln, donor) rows

DESIGN NOTES:
  - VACUUM, NOT GBSA: we don't have amber14/obc2.xml in this OpenMM install,
    and amber99_obc uses integer atom-type indexing that doesn't compose with
    amber14/tip3p's La. Vacuum + restrained backbone is the cleanest minimal
    setup. The +3 over-attraction is uniform across the Ln panel, so
    differential ranking is preserved (LJ σ differences give the Ln-O distance
    spread we want to measure).

  - C4 INCLUDED: this is canonical 12-6-4. May ablate in a v2 if needed.

  - OTHER 2 La RESTRAINED, NOT SWAPPED: we're testing one site at a time;
    the other sites stay at La to keep the protein scaffold near native.

Author: Claude Sonnet 4.7 + Jacob Westrob
Date: 2026-05-03
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from io import StringIO
from multiprocessing import Pool, get_context

import numpy as np

# Limit per-worker OpenMM threading so process-level parallelism wins.
os.environ.setdefault("OPENMM_CPU_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

# Defer heavy imports
ALCH_DIR = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/"
                "lanthanide_binding/on_density_scanner/alchemical_bvs")
SCANNER_DIR = ALCH_DIR.parent
PDB_8DQ2 = SCANNER_DIR / "calibration/pdb_controls/8DQ2.pdb"
LN_PARAMS = ("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/"
             "lanthanide_binding/fep/prep/ln_parameters.json")
SITE_RECON = ALCH_DIR / "inputs/8dq2_chainA_sites.json"

LN_PANEL = ["La", "Nd", "Sm", "Eu", "Tb", "Dy", "Yb", "Lu"]

# Restraint stiffness: 50 kcal/mol/Å² = 50 * 4.184 * 100 = 20920 kJ/mol/nm²
K_RESTRAINT_KJ_NM2 = 50.0 * 4.184 * 100.0  # = 20920

# For BVS computation
sys.path.insert(0, str(SCANNER_DIR))


def load_ln_params() -> dict:
    """Load Merz 12-6-4 OPC parameters. Confirms presence of all panel Ln."""
    with open(LN_PARAMS) as f:
        params = json.load(f)
    missing = [ln for ln in LN_PANEL if ln not in params]
    if missing:
        raise RuntimeError(f"Missing Ln params: {missing}")
    return params


def load_8dq2_chainA_clean(pdb_path: Path):
    """Return prepared topology + positions: chain A protein + 3 La ions,
    no waters, hydrogens added. Avoids any PDB round-trip so original
    residue IDs (201, 202, 203 for La) are preserved.
    Returns (Topology, positions, la_atom_indices_by_resseq)."""
    from pdbfixer import PDBFixer
    from openmm.app import Modeller

    fixer = PDBFixer(str(pdb_path))

    # Keep chain A only.
    chains_to_remove = [c.index for c in fixer.topology.chains() if c.id != "A"]
    fixer.removeChains(chainIndices=chains_to_remove)

    # PDBFixer prep on the protein. addMissingHydrogens uses topology
    # heuristics, not FF templates — robust against the LINK-record
    # Ln-O "bonds" that confuse FF template matching.
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    # Strip waters via Modeller (PDBFixer has no water-only deletion).
    # Use Modeller AFTER addMissingHydrogens so we don't waste time
    # adding H to waters we'd just delete.
    modeller = Modeller(fixer.topology, fixer.positions)
    to_delete = [res for res in modeller.topology.residues() if res.name == "HOH"]
    modeller.delete(to_delete)

    # Identify La indices, preserving original residue IDs (201, 202, 203).
    la_atoms_by_seq = {}
    for residue in modeller.topology.residues():
        if residue.name == "LA":
            for atom in residue.atoms():
                la_atoms_by_seq[int(residue.id)] = atom.index

    # CRITICAL: strip Ln-protein bonds from the topology. The 8DQ2 PDB
    # contains LINK records describing each Ln-O coordination as a covalent
    # bond. OpenMM materializes these as topology bonds, and NonbondedForce
    # auto-creates 1-2 EXCLUSIONS for them — which kills both the LJ wall
    # and the Coulomb attraction between metal and donors. Donors then
    # collapse to ~1.5 Å with no opposing force (verified empirically).
    la_indices_set = set(la_atoms_by_seq.values())
    new_bonds = []
    n_stripped = 0
    for bond in modeller.topology.bonds():
        a1, a2 = bond.atom1, bond.atom2
        if a1.index in la_indices_set or a2.index in la_indices_set:
            n_stripped += 1
            continue
        new_bonds.append(bond)
    modeller.topology._bonds = new_bonds

    return modeller.topology, modeller.positions, la_atoms_by_seq


# Element-based OBC radii (Bondi 1964) and OBC scale factors.
# These are added manually because amber14/obc2.xml doesn't ship with this
# OpenMM install. Without GBSA, vacuum +3 Coulomb pulls the Ln across the
# protein (verified empirically: 7.8 Å metal migration in vacuum smoke test).
OBC_PARAMS_BY_ELEMENT = {
    "H": (0.12, 0.85),
    "C": (0.17, 0.72),
    "N": (0.155, 0.79),
    "O": (0.15, 0.85),
    "S": (0.18, 0.96),
    "P": (0.185, 0.86),
}
# Ln OBC radius: between Sr²⁺ and Cs⁺ Bondi-like estimates; use 0.18 nm
# uniformly across the Ln series (the metal is buried by 10 carboxylates,
# so its GBSA contribution is small via Born-radius screening anyway —
# this value is not a sensitive dial).
OBC_LN_RADIUS = 0.18
OBC_LN_SCALE = 0.85


def build_system_for_site(
    topology, positions, la_atom_indices, target_la_seqid: int,
    target_ln: str, ln_params: dict, free_residue_set: set,
):
    """Build an OpenMM system with GBSA OBC2 implicit solvent for a
    single (site, Ln) combo.

    Returns: (system, target_la_idx, restrained_atom_indices)
    """
    from openmm.app import ForceField, NoCutoff
    from openmm import (CustomNonbondedForce, NonbondedForce,
                        CustomExternalForce, GBSAOBCForce, unit)

    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    # ignoreExternalBonds=True is critical: the PDB has LINK records that
    # define Ln-O coordination as topology bonds, but neither the protein
    # FF nor the LA template knows about these. We tell ForceField to
    # match standard residue templates and ignore inter-residue bonds it
    # doesn't have parameters for.
    system = ff.createSystem(
        topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
        removeCMMotion=True,
        ignoreExternalBonds=True,
    )

    target_la_idx = la_atom_indices[target_la_seqid]
    other_la_indices = {idx for seq, idx in la_atom_indices.items()
                        if seq != target_la_seqid}

    # 1) Set target La's nonbonded params to target Ln (Merz 12-6-4 LJ).
    p_target = ln_params[target_ln]
    rmin_half_A = p_target["Rmin_half"]
    eps_kcal = p_target["epsilon"]
    sigma_nm = rmin_half_A * 2.0 * (2.0 ** (-1.0 / 6.0)) * 0.1
    eps_kj = eps_kcal * 4.184

    nb_force = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    nb_force.setParticleParameters(
        target_la_idx,
        3.0 * unit.elementary_charge,
        sigma_nm * unit.nanometer,
        eps_kj * unit.kilojoule_per_mole,
    )
    nb_force.setUseDispersionCorrection(False)

    # 2) Add C4 CustomNonbondedForce for ion-O ion-induced-dipole.
    # NOTE: OpenMM CustomNonbondedForce appends "1"/"2" directly (no
    # underscore) to per-particle parameter names. So if we name the
    # parameter "C4", it becomes "C41"/"C42" in the expression.
    c4_force = CustomNonbondedForce(
        "-C4prod / (r^4);"
        "C4prod = C41 * isO2 + C42 * isO1"
    )
    c4_force.addPerParticleParameter("C4")
    c4_force.addPerParticleParameter("isO")
    c4_force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)

    # C4 in kcal/mol·Å⁴ → kJ/mol·nm⁴: × 4.184 × 0.1^4 = × 4.184e-4
    c4_target_kjnm4 = p_target["C4"] * 4.184e-4
    # Other La ions retain La's C4
    c4_la_kjnm4 = ln_params["La"]["C4"] * 4.184e-4

    for atom in topology.atoms():
        if atom.index == target_la_idx:
            c4_force.addParticle([c4_target_kjnm4, 0.0])
        elif atom.index in other_la_indices:
            c4_force.addParticle([c4_la_kjnm4, 0.0])
        elif atom.element is not None and atom.element.symbol == "O":
            c4_force.addParticle([0.0, 1.0])
        else:
            c4_force.addParticle([0.0, 0.0])

    # Copy NonbondedForce exceptions to keep 1-2/1-3 exclusions consistent.
    for i in range(nb_force.getNumExceptions()):
        p1, p2, _, _, _ = nb_force.getExceptionParameters(i)
        c4_force.addExclusion(p1, p2)
    system.addForce(c4_force)

    # 2.5) Add GBSAOBCForce manually with element-based OBC radii.
    # Required to screen the +3 Coulomb pull at long range (without GBSA the
    # metal migrates across the protein — verified empirically).
    gbsa = GBSAOBCForce()
    gbsa.setSoluteDielectric(1.0)
    gbsa.setSolventDielectric(78.5)
    gbsa.setNonbondedMethod(GBSAOBCForce.NoCutoff)
    # Charges from NonbondedForce (these were set by the FF for protein,
    # and we just set +3 for the target Ln above).
    for atom in topology.atoms():
        idx = atom.index
        charge_q, _, _ = nb_force.getParticleParameters(idx)
        # Element-based OBC radius/scale; metal special-cased.
        if idx in la_atom_indices.values():
            radius_nm = OBC_LN_RADIUS
            scale = OBC_LN_SCALE
        elif atom.element is not None and atom.element.symbol in OBC_PARAMS_BY_ELEMENT:
            radius_nm, scale = OBC_PARAMS_BY_ELEMENT[atom.element.symbol]
        else:
            radius_nm, scale = 0.15, 0.85  # fallback
        gbsa.addParticle(charge_q, radius_nm * unit.nanometer, scale)
    system.addForce(gbsa)

    # 3) Restraints: backbone always restrained; sidechains restrained UNLESS
    #    in free_residue_set; other Ln ions also restrained at X-ray pos.
    restraint = CustomExternalForce(
        "k_r * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
    )
    restraint.addGlobalParameter(
        "k_r", K_RESTRAINT_KJ_NM2 * unit.kilojoule_per_mole / unit.nanometer**2
    )
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    BACKBONE = {"N", "CA", "C", "O"}
    restrained = []
    for atom in topology.atoms():
        res = atom.residue
        elem = atom.element.symbol if atom.element else None
        if elem == "H":
            continue
        # La ions: restrain non-target ions at X-ray pos; target is free.
        if atom.index in other_la_indices:
            pos = positions[atom.index]
            restraint.addParticle(atom.index, [
                pos[0].value_in_unit(unit.nanometer),
                pos[1].value_in_unit(unit.nanometer),
                pos[2].value_in_unit(unit.nanometer),
            ])
            restrained.append(atom.index)
            continue
        if atom.index == target_la_idx:
            # Target metal: ALSO restrained at X-ray position. GBSA implicit
            # solvent self-solvation pushes +3 ions out of buried sites into
            # bulk dielectric (verified: metal flew 137 Å w/o this). Pinning
            # the metal here lets the σ-driven differential geometry come
            # entirely from donor sidechain rearrangement around the metal —
            # which is what we're trying to measure with BVS at relaxed
            # geometry.
            pos = positions[atom.index]
            restraint.addParticle(atom.index, [
                pos[0].value_in_unit(unit.nanometer),
                pos[1].value_in_unit(unit.nanometer),
                pos[2].value_in_unit(unit.nanometer),
            ])
            restrained.append(atom.index)
            continue
        # Standard residues: backbone always restrained
        if atom.name in BACKBONE:
            pos = positions[atom.index]
            restraint.addParticle(atom.index, [
                pos[0].value_in_unit(unit.nanometer),
                pos[1].value_in_unit(unit.nanometer),
                pos[2].value_in_unit(unit.nanometer),
            ])
            restrained.append(atom.index)
            continue
        # Sidechain: free if residue is in free_residue_set, else restrained
        if (res.name, int(res.id)) in free_residue_set:
            continue
        pos = positions[atom.index]
        restraint.addParticle(atom.index, [
            pos[0].value_in_unit(unit.nanometer),
            pos[1].value_in_unit(unit.nanometer),
            pos[2].value_in_unit(unit.nanometer),
        ])
        restrained.append(atom.index)
    system.addForce(restraint)
    return system, target_la_idx, restrained


def minimize_and_analyze(
    system, topology, init_positions, target_la_idx, target_ln: str,
    ln_params: dict, free_residue_set: set,
):
    """Run L-BFGS minimization, return relaxed positions + diagnostics."""
    from openmm import LangevinMiddleIntegrator, Platform, unit
    from openmm.app import Simulation

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picosecond
    )
    # CPU platform: fast for system size (1668 atoms), no GPU required,
    # deterministic enough for minimization.
    try:
        platform = Platform.getPlatformByName("CPU")
    except Exception:
        platform = Platform.getPlatformByName("Reference")
    sim = Simulation(topology, system, integrator, platform)
    sim.context.setPositions(init_positions)

    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    t0 = time.time()
    sim.minimizeEnergy(
        tolerance=0.4 * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=5000,
    )
    dt = time.time() - t0
    state = sim.context.getState(getPositions=True, getEnergy=True)
    e1 = state.getPotentialEnergy()
    return state.getPositions(), e0, e1, dt


def compute_first_shell_donors(
    topology, positions, target_la_idx, free_residue_set: set,
    cutoff_A: float = 4.0,
):
    """Return Atom-like list of donors (O, N) within cutoff_A of metal."""
    sys.path.insert(0, str(SCANNER_DIR))
    from on_scanner.geometry import Atom
    from openmm import unit

    metal_pos = np.array(positions[target_la_idx].value_in_unit(unit.angstrom))
    donors = []
    for atom in topology.atoms():
        if atom.index == target_la_idx:
            continue
        elem = atom.element.symbol if atom.element else None
        if elem not in ("O", "N"):
            continue
        if atom.name in ("N", "C"):  # backbone amide N — usually too far for Ln; allow
            pass
        pos = np.array(positions[atom.index].value_in_unit(unit.angstrom))
        d = float(np.linalg.norm(pos - metal_pos))
        if d > cutoff_A:
            continue
        donors.append(Atom(
            pos=pos,
            element=elem,
            atom_name=atom.name,
            res_name=atom.residue.name,
            res_seq=int(atom.residue.id),
            chain=atom.residue.chain.id,
            kind="bb_O" if atom.name == "O" else "sc",
        ))
    return donors, metal_pos


def run_one_combo(site_id: int, target_ln: str, ln_params, recon, top, pos, la_idx_map):
    """Run minimization + BVS for one (site, Ln) combo."""
    site_key = f"LA{site_id}"
    site_info = recon[site_key]
    free_residue_set = {(rn, rs) for rn, rs in site_info["free_residues"]}

    # X-ray metal pos (for RMSD reference)
    xray_metal_pos = np.array(site_info["metal_pos_A"])

    # Build system
    system, target_la_idx, restrained = build_system_for_site(
        top, pos, la_idx_map, site_id, target_ln, ln_params, free_residue_set
    )

    # Minimize
    relaxed_positions, e0, e1, dt = minimize_and_analyze(
        system, top, pos, target_la_idx, target_ln, ln_params, free_residue_set
    )

    # First-shell donors at relaxed geometry
    donors, relaxed_metal_pos = compute_first_shell_donors(
        top, relaxed_positions, target_la_idx, free_residue_set, cutoff_A=4.0
    )

    # Compute BVS for the target Ln cation
    from on_scanner.bvs import compute_bvs
    target_cation = f"{target_ln}3+"
    bvs_result = compute_bvs(relaxed_metal_pos, donors, target_cation, cutoff=4.0)

    # Diagnostics
    metal_displacement_A = float(np.linalg.norm(relaxed_metal_pos - xray_metal_pos))

    return {
        "site_id": site_id,
        "target_ln": target_ln,
        "n_donors_first_shell": bvs_result.n_donors_used,
        "bvs": bvs_result.bvs,
        "bvs_residual": bvs_result.residual,
        "bvs_per_donor": bvs_result.bvs_per_donor,
        "donor_distances_A": [
            (atom_id, dist)
            for atom_id, dist, _ in bvs_result.per_donor_contributions
        ],
        "metal_displacement_from_xray_A": metal_displacement_A,
        "energy_init_kJ": e0.value_in_unit_system(__import__("openmm").unit.md_unit_system),
        "energy_final_kJ": e1.value_in_unit_system(__import__("openmm").unit.md_unit_system),
        "minimize_time_s": dt,
        "n_restrained_atoms": len(restrained),
    }


def _emit_result(status, result, rows, donor_rows, i, n_total):
    """Format + accumulate one combo's result. Mutates rows/donor_rows."""
    if status == "FAIL":
        print(f"[{i+1}/{n_total}] FAILED  site LA{result['site_id']} {result['target_ln']}: {result['error']}")
        rows.append({**result, "status": "FAILED"})
        return
    r = result
    print(f"[{i+1}/{n_total}] OK  site LA{r['site_id']} {r['target_ln']}: "
          f"BVS={r['bvs']:.3f} resid={r['bvs_residual']:.3f} "
          f"n_donors={r['n_donors_first_shell']} dmin={min(d for _,d in r['donor_distances_A']):.3f}Å "
          f"({r['minimize_time_s']:.1f}s)")
    for atom_id, dist in r["donor_distances_A"]:
        donor_rows.append({
            "site_id": r["site_id"], "target_ln": r["target_ln"],
            "donor_atom": atom_id, "distance_A": dist,
        })
    rr = {k: v for k, v in r.items() if k != "donor_distances_A"}
    rr["status"] = "OK"
    rows.append(rr)


def _worker(args):
    """Worker fn for multiprocessing. Reloads everything fresh per combo —
    avoids pickling OpenMM objects (Topology isn't reliably picklable).
    Cost is ~1s of overhead per combo, well below the minimization itself."""
    site_id, target_ln = args
    ln_params = load_ln_params()
    with open(SITE_RECON) as f:
        recon = json.load(f)
    top, pos, la_idx_map = load_8dq2_chainA_clean(PDB_8DQ2)
    try:
        return ("OK", run_one_combo(site_id, target_ln, ln_params, recon,
                                     top, pos, la_idx_map))
    except Exception as e:
        return ("FAIL", {"site_id": site_id, "target_ln": target_ln,
                          "error": f"{type(e).__name__}: {e}"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default=str(ALCH_DIR / "results"))
    parser.add_argument("--sites", nargs="*", default=None,
                        help="Restrict to specific site IDs (e.g., 201 202)")
    parser.add_argument("--lns", nargs="*", default=None,
                        help="Restrict to specific Ln (e.g., La Eu Yb)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers (default: $SLURM_CPUS_ON_NODE or os.cpu_count())")
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_workers = args.workers or int(
        os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1)
    )

    print("=" * 70)
    print("Pre-check: alchemical BVS on 8DQ2 chain A")
    print("=" * 70)
    print(f"Output dir: {out_dir}")
    print(f"Workers: {n_workers}")

    sites_to_run = [int(s) for s in args.sites] if args.sites else [201, 202, 203]
    lns_to_run = args.lns if args.lns else LN_PANEL
    combos = [(s, ln) for s in sites_to_run for ln in lns_to_run]
    print(f"Combos: {len(combos)} = {len(sites_to_run)} sites × {len(lns_to_run)} Ln")

    # Sanity-load on the parent process so we fail-fast on a missing input.
    _ = load_ln_params()
    with open(SITE_RECON) as f:
        _ = json.load(f)

    t_start = time.time()
    rows = []
    donor_rows = []
    if n_workers > 1 and len(combos) > 1:
        # 'spawn' avoids fork-related issues with OpenMM's CPU platform.
        ctx = get_context("spawn")
        with ctx.Pool(n_workers) as pool:
            for i, (status, result) in enumerate(
                pool.imap_unordered(_worker, combos)
            ):
                _emit_result(status, result, rows, donor_rows, i, len(combos))
    else:
        for i, combo in enumerate(combos):
            status, result = _worker(combo)
            _emit_result(status, result, rows, donor_rows, i, len(combos))
    print(f"\n[done in {time.time()-t_start:.1f}s wallclock]")

    # Write per-combo TSV
    import csv
    tsv_path = out_dir / "precheck_8dq2_residuals.tsv"
    with open(tsv_path, "w", newline="") as f:
        cols = ["site_id", "target_ln", "status", "n_donors_first_shell",
                "bvs", "bvs_residual", "bvs_per_donor",
                "metal_displacement_from_xray_A",
                "energy_init_kJ", "energy_final_kJ", "minimize_time_s",
                "n_restrained_atoms", "error"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t",
                            extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nWrote {tsv_path}")

    donor_path = out_dir / "precheck_8dq2_donors.tsv"
    with open(donor_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["site_id", "target_ln", "donor_atom",
                                           "distance_A"], delimiter="\t")
        w.writeheader()
        for r in donor_rows:
            w.writerow(r)
    print(f"Wrote {donor_path}")

    # Pretty summary: per-site Ln residual vs target (closest = predicted Ln)
    print("\n" + "=" * 70)
    print("SUMMARY — predicted Ln per site (smallest |BVS-q|)")
    print("=" * 70)
    print(f"{'Site':<8s}" + "".join(f"{ln:>10s}" for ln in lns_to_run) + "  predicted")
    for site_id in sites_to_run:
        row_str = f"LA{site_id:<6d}"
        site_rows = {r["target_ln"]: r["bvs_residual"] for r in rows
                     if r.get("site_id") == site_id and r.get("status") == "OK"}
        for ln in lns_to_run:
            v = site_rows.get(ln, float("nan"))
            row_str += f"{v:>10.3f}"
        if site_rows:
            best = min(site_rows, key=site_rows.get)
            row_str += f"  → {best}"
        print(row_str)


if __name__ == "__main__":
    main()
