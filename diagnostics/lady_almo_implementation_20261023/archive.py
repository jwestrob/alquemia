"""Exact archived quartet inputs, parsed without importing a quantum engine.

The metadata CLI is safe preparation. ``build_molecule`` is used only in a
root-submitted Slurm job. No CP2K-converted orbital or ECP values are imported.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("LADY_REPO_ROOT", HERE.parent.parent)).resolve()
OLD = ROOT / "diagnostics/lady_electronic_decomposition_20261023"
PROVENANCE = OLD / "inputs/PROVENANCE.json"
LEDGER = OLD / "cp2k_inputs_v1/ENERGY_LEDGER.json"
BSE = OLD / "cp2k_inputs_v1/assets/def2-TZVP_HCO_BSE_v1.json"
MANIFEST = ROOT / "diagnostics/lady_almo_implementation_20261023/INPUT_PROVENANCE.json"
SYSTEMS = ("aquo8_central_La", "aquo8_central_Dy",
           "hydrated_acetate_mono_central_La", "hydrated_acetate_mono_central_Dy")
Z = {"H": 1, "C": 6, "O": 8, "La": 57, "Dy": 66}
ANGULAR = {s: l for l, s in enumerate("spdfghik")}
HARTREE_TO_KCAL = 627.509474


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": sha256(path)}


def verify_pin(record):
    if sha256(record["path"]) != record["sha256"]:
        raise ValueError(f"Input hash changed: {record['path']}")
    return Path(record["path"])


def _number(token):
    result = float(token.replace("D", "E").replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError("Nonfinite archived coefficient")
    return result


def _section(text, keyword, metal):
    matches = re.findall(rf"(?im)^\s*{keyword}\s+{metal}\s*\n(.*?)^\s*end\s*$",
                         text, flags=re.S | re.M | re.I)
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {keyword} {metal} section")
    return [line.split() for line in matches[0].splitlines() if line.strip()]


def parse_orca_metal(text, metal):
    """Translate only the exact archived semilocal ECP and spherical basis.

    ORCA ECP rows are index, exponent, coefficient, radial n; potential terms
    are c*r**(n-2)*exp(-a*r*r). PySCF ECP stores n in the radial-power bucket.
    The highest (g) channel is local, encoded as l=-1 in PySCF. The original
    g term is zero, so no nonlocal-channel subtraction is required. A nonzero
    local channel is deliberately rejected rather than silently generalized.
    """
    rows = _section(text, "NewGTO", metal)
    basis, raw_basis = [], []
    i = 0
    while i < len(rows):
        symbol, count = rows[i]
        l, count = ANGULAR[symbol.lower()], int(count)
        terms = rows[i + 1:i + 1 + count]
        if len(terms) != count:
            raise ValueError("Incomplete orbital shell")
        shell = [l]
        for j, row in enumerate(terms, 1):
            if len(row) != 3 or int(row[0]) != j:
                raise ValueError("Unexpected raw ORCA orbital row")
            exponent, coefficient = map(_number, row[1:])
            if exponent <= 0:
                raise ValueError("Nonpositive orbital exponent")
            shell.append([exponent, coefficient])
        basis.append(shell)
        raw_basis.append({"l": l, "rows": terms})
        i += count + 1
    rows = _section(text, "NewECP", metal)
    if rows[0][0].lower() != "n_core" or rows[1][0].lower() != "lmax":
        raise ValueError("Unexpected ECP header")
    core, lmax = int(rows[0][1]), ANGULAR[rows[1][1].lower()]
    if core != {"La": 46, "Dy": 55}[metal] or lmax != 4:
        raise ValueError("Archived core count / local channel changed")
    channels, raw_ecp = {}, []
    i = 2
    while i < len(rows):
        symbol, count = rows[i]
        l, count = ANGULAR[symbol.lower()], int(count)
        if l in channels or l > lmax:
            raise ValueError("Duplicate or unexpected ECP channel")
        terms = rows[i + 1:i + 1 + count]
        if len(terms) != count:
            raise ValueError("Incomplete ECP channel")
        buckets = [[] for _ in range(7)]
        for j, row in enumerate(terms, 1):
            if len(row) != 4 or int(row[0]) != j:
                raise ValueError("Unexpected raw ORCA ECP row")
            exponent, coefficient = map(_number, row[1:3])
            n = int(row[3])
            if n != 2 or exponent <= 0:
                raise ValueError("Only archived n=2 ECP terms are admitted")
            if l == lmax and coefficient != 0:
                raise ValueError("Nonzero local ECP requires separate conversion review")
            buckets[n].append([exponent, coefficient])
        channels[l] = buckets
        raw_ecp.append({"l": l, "local": l == lmax, "rows": terms})
        i += count + 1
    if set(channels) != set(range(lmax + 1)):
        raise ValueError("Missing ECP channel")
    ecp = [core, [[-1, channels[lmax]]] + [[l, channels[l]] for l in range(lmax)]]
    if sum(2 * shell[0] + 1 for shell in basis) != 47:
        raise ValueError("Metal spherical basis count differs from archive")
    return basis, ecp, {"raw_basis": raw_basis, "raw_ecp": raw_ecp,
                        "ecp_core_electrons": core, "local_l": lmax,
                        "radial_convention": "c*r^(n-2)*exp(-a*r^2); n=2 retained"}


def parse_bse_light_atoms(data):
    if data["name"].lower() != "def2-tzvp" or str(data["version"]) != "1":
        raise ValueError("Unexpected light-atom basis source")
    basis = {}
    for symbol in ("H", "C", "O"):
        shells = []
        for shell in data["elements"][str(Z[symbol])]["electron_shells"]:
            if len(shell["angular_momentum"]) != 1 or len(shell["coefficients"]) != 1:
                raise ValueError("Unexpected contracted BSE shell representation")
            exp, coef = shell["exponents"], shell["coefficients"][0]
            if len(exp) != len(coef):
                raise ValueError("Mismatched BSE exponents and coefficients")
            shells.append([shell["angular_momentum"][0]] +
                          [[_number(a), _number(c)] for a, c in zip(exp, coef)])
        basis[symbol] = shells
    expected = {"H": 6, "C": 31, "O": 31}
    if {s: sum(2 * sh[0] + 1 for sh in shells) for s, shells in basis.items()} != expected:
        raise ValueError("Light-atom spherical basis counts differ from archive")
    return basis


def load_system(system_id, *, verify_manifest=True):
    """Return the exact quartet member as pure Python input and provenance."""
    if system_id not in SYSTEMS:
        raise ValueError(f"Outside declared quartet: {system_id}")
    if verify_manifest:
        manifest = json.loads(MANIFEST.read_text())
        for record in manifest["input_pins"]:
            verify_pin(record)
        verify_pin(manifest["loader_source"])
    provenance = json.loads(PROVENANCE.read_text())
    saved = provenance["systems"][system_id]
    xyz_path = verify_pin(saved["xyz"])
    xyz = xyz_path.read_text().splitlines()
    count = int(xyz[0])
    rows = [line.split() for line in xyz[2:] if line.strip()]
    if count != saved["atoms"] or len(rows) != count or any(len(r) != 4 for r in rows):
        raise ValueError("XYZ count/format mismatch")
    atoms = [(r[0], tuple(map(_number, r[1:]))) for r in rows]
    metal = saved["metal"]
    if atoms[0][0] != metal or any(a[0] not in ("H", "C", "O") for a in atoms[1:]):
        raise ValueError("Archived atom-order / identity contract violated")
    asset = verify_pin(provenance["files"]["metal_assets"][metal]["basis"]["copy"])
    basis_m, ecp, conversion = parse_orca_metal(asset.read_text(), metal)
    native = verify_pin(saved["original_native_input"])
    if parse_orca_metal(native.read_text(), metal) != (basis_m, ecp, conversion):
        raise ValueError("Raw metal asset does not match executed ORCA input")
    basis = parse_bse_light_atoms(json.loads(BSE.read_text()))
    basis[metal] = basis_m
    fragments = []
    for f in saved["fragments"]:
        if f["effective_multiplicity"] != 1:
            raise ValueError("Only archived effective singlets are admitted")
        fragments.append({"label": f["label"],
                          "atom_indices": [i - 1 for i in f["atom_indices_one_based"]],
                          "charge": f["charge"], "spin": 0,
                          "nelectron": f["explicit_electrons"]})
    if sorted(i for f in fragments for i in f["atom_indices"]) != list(range(count)):
        raise ValueError("Fragment partition is not exact")
    electrons = sum(Z[a[0]] for a in atoms) - ecp[0] - saved["charge"]
    if electrons != saved["explicit_electrons"] or electrons != sum(f["nelectron"] for f in fragments):
        raise ValueError("Explicit electron count does not close")
    if sum(f["charge"] for f in fragments) != saved["charge"]:
        raise ValueError("Fragment charges do not close")
    if saved["effective_multiplicity"] != 1:
        raise ValueError("Only archived effective singlets are admitted")
    energy = json.loads(LEDGER.read_text())["systems"][system_id]["energy_hartree"]
    return {"system_id": system_id, "atoms": atoms, "charge": saved["charge"],
            "spin": 0, "nelectron": electrons, "basis": basis, "ecp": {metal: ecp},
            "fragments": fragments, "metal": metal, "basis_conversion": conversion,
            "physical_multiplicity_metadata": saved["physical_multiplicity_metadata"],
            "physical_state_qualified": False,
            "expected_nao": 391 if count == 25 else 490,
            "archive_energies": {
                "electronic_hartree": float(energy["electronic_after_final_cosx_integration"]["token"]),
                "full_pbe0_d4_hartree": float(energy["full_pbe0_d4"]["token"]),
                "d4_printed_hartree": float(energy["d4_printed"]["token"]),
                "d4_effective_token_difference_hartree": float(energy["d4_effective_from_full_minus_electronic_token_difference"]),
                "note": "Only archived printed D4 is added once; effective token difference retained as bookkeeping, not a new dispersion calculation."},
            "provenance": {"xyz": saved["xyz"], "native_input": saved["original_native_input"],
                           "metal_asset": pin(asset), "light_basis": pin(BSE),
                           "energy_ledger": pin(LEDGER), "original_output": saved["original_output"]}}


def build_molecule(system_id, output, max_memory_mb, fragment=None):
    """Build a spherical PySCF Mole; numerical jobs must be Slurm allocations."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Molecular setup is restricted to root-submitted Slurm jobs")
    from pyscf import gto
    info = load_system(system_id)
    selected = list(range(len(info["atoms"])))
    charge, spin, nelectron = info["charge"], info["spin"], info["nelectron"]
    if fragment is not None:
        f = info["fragments"][fragment]
        selected = f["atom_indices"]
        charge, spin, nelectron = f["charge"], f["spin"], f["nelectron"]
    atoms = [info["atoms"][i] for i in selected]
    present = {a[0] for a in atoms}
    mol = gto.Mole()
    mol.atom = atoms
    mol.unit = "Angstrom"
    mol.charge, mol.spin = charge, spin
    mol.basis = {s: copy.deepcopy(info["basis"][s]) for s in present}
    mol.ecp = {s: copy.deepcopy(v) for s, v in info["ecp"].items() if s in present}
    mol.cart = False
    mol.symmetry = False
    mol.max_memory = max_memory_mb
    mol.verbose = 4
    mol.output = str(output)
    mol.build()
    expected_nao = sum(sum(2 * sh[0] + 1 for sh in info["basis"][a[0]]) for a in atoms)
    if mol.nelectron != nelectron or mol.spin != 0 or mol.nao_nr() != expected_nao:
        raise ValueError("PySCF runtime atom/electron/spin/AO count does not match archive")
    for i, atom in enumerate(atoms):
        expected_core = info["ecp"].get(atom[0], [0])[0]
        if mol.atom_nelec_core(i) != expected_core:
            raise ValueError("PySCF runtime core count mismatch")
    info["runtime_subset"] = {"fragment": fragment, "atom_indices": selected,
                              "charge": charge, "spin": spin, "nelectron": nelectron,
                              "nao": mol.nao_nr(), "cartesian": mol.cart}
    return mol, info


def write_provenance():
    """Write source-only declaration; refuse to replace existing provenance."""
    provenance = json.loads(PROVENANCE.read_text())
    paths = {PROVENANCE, LEDGER, BSE}
    for saved in provenance["systems"].values():
        for key in ("xyz", "original_native_input", "original_output"):
            paths.add(verify_pin(saved[key]))
    for metal in ("La", "Dy"):
        paths.add(verify_pin(provenance["files"]["metal_assets"][metal]["basis"]["copy"]))
    systems = {name: load_system(name, verify_manifest=False) for name in SYSTEMS}
    record = {"schema": "archived_quartet_pyscf_translation_v1",
              "status": "source_and_metadata_prepared_not_numerically_validated",
              "new_molecular_calls": 0, "loader_source": pin(__file__),
              "input_pins": [pin(p) for p in sorted(paths)], "systems": systems,
              "light_basis_evidence": "BSE def2-TZVP v1 named-library equivalence; original ORCA light-atom primitive coefficients were not printed.",
              "ecp_mapping": "Direct raw ORCA: g(lmax=4) zero channel -> local l=-1; s/p/d/f nonlocal channels; n=2 radial index retained. No CP2K converted asset used.",
              "basis_normalization": "PySCF spherical contracted Gaussian primitive/contraction normalization; coefficient tokens preserved, no decontraction.",
              "numerical_comparison_required": True}
    with MANIFEST.open("x") as stream:
        json.dump(record, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-provenance", action="store_true", required=True)
    parser.parse_args()
    record = write_provenance()
    print(json.dumps({"status": record["status"], "systems": list(record["systems"]),
                      "manifest": str(MANIFEST), "new_molecular_calls": 0}))
