#!/usr/bin/env python3
"""Run the pinned protonator without downloading definitions for PQQ or La.

PDBFixer 1.12 attempts to download CCD hydrogen definitions for every
nonstandard residue before protonating standard amino acids.  That is neither
needed nor wanted here: PQQ hydrogens belong exclusively to the pinned PQQ(3-)
carver specification, and La has no hydrogens.  This wrapper temporarily makes
that one private lookup return an empty mapping while retaining the canonical
protonator and OpenMM's standard-residue hydrogen definitions unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random
import sys
from typing import Any, Mapping, Sequence

import gemmi
from openmm import Platform, app
import pdbfixer as pdbfixer_module
from pdbfixer import PDBFixer as OriginalPDBFixer


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SCRIPTS = PROJECT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import protonate_cif as canonical  # noqa: E402


PROTOCOL_ID = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
POLICY_ID = "standard_residues_only_no_nonstandard_CCD_hydrogens_v1"
PROTONATION_SUBPROTOCOL_ID = (
    "pdbfixer_standard_only_rng20260914_openmm_cpu_threads1_v1"
)
RNG_SEED = 20260914
STANDARD_PROTEIN_RESIDUES = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
)
EXPECTED_NONSTANDARD = {"LA": 1, "PQQ": 1}
EXPECTED_SUBPROTOCOL = {
    "id": PROTONATION_SUBPROTOCOL_ID,
    "nonstandard_definition_policy_id": POLICY_ID,
    "python_random_seed": RNG_SEED,
    "openmm_platform": "CPU",
    "openmm_cpu_threads": 1,
    "forcefield": None,
    "allowed_nonstandard_residue_counts": EXPECTED_NONSTANDARD,
    "required_output_nonstandard_hydrogen_counts": {"LA": 0, "PQQ": 0},
}


class StandardOnlyProtonationError(RuntimeError):
    """The narrow standard-residue-only protonation contract cannot be met."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise StandardOnlyProtonationError(f"cannot read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise StandardOnlyProtonationError(f"JSON root is not an object: {path}")
    return value


def verify_file_record(record: Any, label: str) -> Path:
    if not isinstance(record, dict):
        raise StandardOnlyProtonationError(f"{label} is not a file record")
    path = Path(str(record.get("path", ""))).resolve()
    expected = record.get("sha256")
    if not path.is_file() or not isinstance(expected, str):
        raise StandardOnlyProtonationError(f"{label} file/hash is missing")
    if sha256_file(path) != expected:
        raise StandardOnlyProtonationError(f"{label} SHA-256 mismatch")
    return path


def write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def inspect_normalized_input(path: Path) -> dict[str, int]:
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise StandardOnlyProtonationError("normalized input does not contain exactly one model")
    residue_counts: dict[str, int] = {}
    pqq_hydrogens = 0
    for chain in structure[0]:
        for residue in chain:
            name = residue.name.upper().strip()
            residue_counts[name] = residue_counts.get(name, 0) + 1
            if name == "PQQ":
                pqq_hydrogens += sum(
                    atom.element.name.upper() in {"H", "D"} for atom in residue
                )
    unexpected = sorted(
        name
        for name in residue_counts
        if name not in STANDARD_PROTEIN_RESIDUES and name not in EXPECTED_NONSTANDARD
    )
    observed_nonstandard = {
        name: residue_counts.get(name, 0) for name in EXPECTED_NONSTANDARD
    }
    if unexpected:
        raise StandardOnlyProtonationError(
            f"normalized input contains unsupported nonstandard residues: {unexpected}"
        )
    if observed_nonstandard != EXPECTED_NONSTANDARD:
        raise StandardOnlyProtonationError(
            f"expected exactly one LA and one PQQ, observed {observed_nonstandard}"
        )
    if pqq_hydrogens:
        raise StandardOnlyProtonationError(
            "normalized PQQ already has hydrogens; the carver must supply them"
        )
    return observed_nonstandard


def nonstandard_hydrogen_counts(path: Path) -> dict[str, int]:
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise StandardOnlyProtonationError("protonated output does not contain exactly one model")
    counts = {name: 0 for name in EXPECTED_NONSTANDARD}
    residue_counts = {name: 0 for name in EXPECTED_NONSTANDARD}
    for chain in structure[0]:
        for residue in chain:
            name = residue.name.upper().strip()
            if name not in EXPECTED_NONSTANDARD:
                continue
            residue_counts[name] += 1
            counts[name] += sum(
                atom.element.name.upper() in {"H", "D"} for atom in residue
            )
    if residue_counts != EXPECTED_NONSTANDARD:
        raise StandardOnlyProtonationError(
            f"protonated output changed PQQ/LA residue counts: {residue_counts}"
        )
    return counts


def protonate_standard_only(
    source: Path,
    output: Path,
    *,
    pins_path: Path,
    ph: float = 7.0,
) -> Path:
    pins_path = pins_path.resolve()
    pins = read_object(pins_path)
    if pins.get("protocol_id") != PROTOCOL_ID:
        raise StandardOnlyProtonationError("implementation pins carry another protocol ID")
    if pins.get("protonation_subprotocol") != EXPECTED_SUBPROTOCOL:
        raise StandardOnlyProtonationError("protonation subprotocol pins changed")
    canonical_record = pins.get("canonical_helpers", {}).get("protonate_cif")
    canonical_path = verify_file_record(canonical_record, "canonical protonator")
    if canonical_path != Path(canonical.__file__).resolve():
        raise StandardOnlyProtonationError("imported protonator differs from the frozen helper")
    own_record = pins.get("experiment_helpers", {}).get("protonate_standard_only")
    own_path = verify_file_record(own_record, "standard-only protonation wrapper")
    if own_path != Path(__file__).resolve():
        raise StandardOnlyProtonationError("executing wrapper differs from its frozen path")
    observed_nonstandard = inspect_normalized_input(source)

    cpu_platform = Platform.getPlatformByName("CPU")
    calls = {"constructor": 0, "add_missing_hydrogens": 0, "ccd_suppression": 0}

    class StandardOnlyPDBFixer(OriginalPDBFixer):
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            requested = kwargs.get("platform")
            if requested is not None and requested.getName() != "CPU":
                raise StandardOnlyProtonationError("only the explicit CPU platform is allowed")
            kwargs["platform"] = cpu_platform
            calls["constructor"] += 1
            super().__init__(*args, **kwargs)

        def _downloadNonstandardDefinitions(self) -> dict[str, Any]:
            # Preserve the useful side effect of PDBFixer's original method.
            # Modeller.addHydrogens() also performs this load, but doing it here
            # makes the only changed behavior the nonstandard CCD download loop.
            calls["ccd_suppression"] += 1
            app.Modeller._loadStandardHydrogenDefinitions()
            return {}

        def addMissingHydrogens(self, pH: float = 7.0, forcefield: Any = None) -> Any:
            if forcefield is not None:
                raise StandardOnlyProtonationError(
                    "fixed protonation subprotocol requires forcefield=None"
                )
            if cpu_platform.getPropertyDefaultValue("Threads") != "1":
                raise StandardOnlyProtonationError(
                    "OpenMM CPU platform is not pinned to one thread"
                )
            calls["add_missing_hydrogens"] += 1
            return super().addMissingHydrogens(pH=pH, forcefield=None)

    original_class = pdbfixer_module.PDBFixer
    rng_state = random.getstate()
    original_cpu_threads = cpu_platform.getPropertyDefaultValue("Threads")
    random.seed(RNG_SEED)
    cpu_platform.setPropertyDefaultValue("Threads", "1")
    pdbfixer_module.PDBFixer = StandardOnlyPDBFixer
    try:
        canonical.protonate(
            source,
            output,
            ph=ph,
            add_missing_residues=False,
        )
    finally:
        pdbfixer_module.PDBFixer = original_class
        cpu_platform.setPropertyDefaultValue("Threads", original_cpu_threads)
        random.setstate(rng_state)

    if calls != {"constructor": 1, "add_missing_hydrogens": 1, "ccd_suppression": 1}:
        raise StandardOnlyProtonationError(
            f"unexpected PDBFixer wrapper call ledger: {calls}"
        )
    output_nonstandard_h = nonstandard_hydrogen_counts(output)
    if output_nonstandard_h != {"LA": 0, "PQQ": 0}:
        raise StandardOnlyProtonationError(
            f"PDBFixer added a nonstandard-residue hydrogen: {output_nonstandard_h}"
        )

    manifest_path = canonical.protonation_manifest_path(output)
    manifest = read_object(manifest_path)
    manifest["canonical_protonator"] = canonical_record
    manifest["experiment_wrapper"] = own_record
    manifest["experiment_protonation_protocol_id"] = PROTONATION_SUBPROTOCOL_ID
    manifest["nonstandard_definition_policy"] = {
        "policy_id": POLICY_ID,
        "suppressed_method": "PDBFixer._downloadNonstandardDefinitions",
        "replacement_return_value": {},
        "standard_OpenMM_hydrogen_definitions_loaded": True,
        "observed_nonstandard_residue_counts": observed_nonstandard,
        "output_nonstandard_hydrogen_counts": output_nonstandard_h,
        "PQQ_hydrogen_policy": "none_added_by_PDBFixer_all_supplied_by_pinned_carver",
        "La_hydrogen_policy": "none",
        "forcefield": None,
        "python_random_seed": RNG_SEED,
        "python_random_state_restored_after_call": True,
        "openmm_platform": "CPU",
        "openmm_CPU_Threads": 1,
        "openmm_CPU_default_threads_restored_after_call": True,
        "wrapper_call_ledger": calls,
        "scientific_chemistry_changed": False,
    }
    write_json_atomic(manifest_path, manifest)
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--implementation-pins", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = protonate_standard_only(
            args.source,
            args.output,
            pins_path=args.implementation_pins,
            ph=args.ph,
        )
    except (StandardOnlyProtonationError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
