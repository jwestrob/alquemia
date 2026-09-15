#!/usr/bin/env python3
"""Prepare the frozen 25-protein fixed-core PQQ-MDH calibration panel.

This is an explicit preparation command, not an import-time action.  It
normalizes, protonates, carves, and validates immutable La/Ca input pairs.  It
never launches ORCA and refuses to reuse or overwrite a partial preparation.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Mapping, Sequence

import gemmi


HERE = Path(__file__).resolve().parent
PROTOCOL_ID = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
PREPARATION_SCHEMA = "alchemical_bvs.pqq_fixed_core_preparation.v1"
NORMALIZATION_SCHEMA = "alchemical_bvs.pqq_fixed_core_normalization.v1"
HEAVY_CHECK_SCHEMA = "alchemical_bvs.heavy_coordinate_check.v1"
EXPECTED_COUNTS = {"La": 11, "Ca": 14}
EXPECTED_ORCA_DIRECTIVE = "! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3"
PDB_SERIALIZATION_TOLERANCE_A = 0.001


class PreparationError(RuntimeError):
    """A frozen input, runtime, or prepared artifact violates the contract."""


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
        raise PreparationError(f"cannot read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PreparationError(f"JSON root is not an object: {path}")
    return value


def write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def file_record(path: Path) -> dict[str, str]:
    path = path.resolve()
    if not path.is_file():
        raise PreparationError(f"required file is missing: {path}")
    return {"path": str(path), "sha256": sha256_file(path)}


def verify_file_record(record: Any, label: str) -> Path:
    if not isinstance(record, dict):
        raise PreparationError(f"{label} record is not an object")
    raw_path = record.get("path")
    expected = record.get("sha256")
    if not isinstance(raw_path, str) or not isinstance(expected, str):
        raise PreparationError(f"{label} lacks path/SHA-256")
    path = Path(raw_path).resolve()
    if not path.is_file():
        raise PreparationError(f"{label} is missing: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise PreparationError(
            f"{label} SHA-256 mismatch: expected {expected}, observed {observed}"
        )
    return path


def load_and_verify_pins(path: Path) -> dict[str, Any]:
    pins = read_object(path)
    if pins.get("schema_version") != "alchemical_bvs.pqq_fixed_core_pins.v1":
        raise PreparationError("unsupported implementation-pins schema")
    if pins.get("protocol_id") != PROTOCOL_ID:
        raise PreparationError("implementation pins carry another protocol ID")
    for key in (
        "preregistration",
        "core_map",
        "core_map_audit",
        "reserved_holdout_spec",
        "reserved_holdout_audit",
        "preparation_amendment",
    ):
        verify_file_record(pins.get(key), f"frozen {key}")
    for group in ("canonical_helpers", "experiment_helpers"):
        records = pins.get(group)
        if not isinstance(records, dict) or not records:
            raise PreparationError(f"implementation pins lack {group}")
        for name, record in records.items():
            verify_file_record(record, f"{group}.{name}")
    own_record = pins["experiment_helpers"].get("prepare")
    if (
        not isinstance(own_record, dict)
        or Path(str(own_record.get("path", ""))).resolve() != Path(__file__).resolve()
    ):
        raise PreparationError("preparation driver is not self-bound in implementation pins")
    verify_file_record(pins.get("aquo_reference"), "pinned aquo reference")
    orca_runtime = pins.get("orca_runtime")
    if not isinstance(orca_runtime, dict):
        raise PreparationError("implementation pins lack ORCA runtime")
    verify_file_record(orca_runtime.get("executable"), "pinned ORCA executable")

    runtime = pins.get("preparation_runtime")
    if not isinstance(runtime, dict):
        raise PreparationError("implementation pins lack preparation runtime")
    python_record = runtime.get("python")
    if not isinstance(python_record, dict):
        raise PreparationError("preparation runtime lacks Python record")
    expected_python = Path(str(python_record.get("path", ""))).resolve()
    if Path(sys.executable).resolve() != expected_python:
        raise PreparationError(
            f"run with pinned interpreter {expected_python}, not {Path(sys.executable).resolve()}"
        )
    observed_python = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    if observed_python != python_record.get("version"):
        raise PreparationError("Python version differs from frozen runtime")
    expected_packages = runtime.get("packages")
    observed_packages = {
        "gemmi": getattr(gemmi, "__version__", None),
        "openmm": importlib.metadata.version("openmm"),
        "pdbfixer": importlib.metadata.version("pdbfixer"),
    }
    if observed_packages != expected_packages:
        raise PreparationError(
            f"package versions {observed_packages} differ from {expected_packages}"
        )
    subprotocol = pins.get("protonation_subprotocol")
    if (
        not isinstance(subprotocol, dict)
        or subprotocol.get("id")
        != "pdbfixer_standard_only_rng20260914_openmm_cpu_threads1_v1"
        or subprotocol.get("python_random_seed") != 20260914
        or subprotocol.get("openmm_platform") != "CPU"
        or subprotocol.get("openmm_cpu_threads") != 1
        or subprotocol.get("forcefield") is not None
        or subprotocol.get("allowed_nonstandard_residue_counts")
        != {"LA": 1, "PQQ": 1}
        or subprotocol.get("required_output_nonstandard_hydrogen_counts")
        != {"LA": 0, "PQQ": 0}
    ):
        raise PreparationError("protonation subprotocol pins are incomplete or changed")
    return pins


def load_core_map(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ids = [row.get("id", "") for row in rows]
    counts = {label: sum(row.get("class") == label for row in rows) for label in EXPECTED_COUNTS}
    if len(rows) != 25 or len(set(ids)) != 25 or counts != EXPECTED_COUNTS:
        raise PreparationError(
            f"core map is not the frozen 25-member panel: n={len(rows)}, counts={counts}"
        )
    return rows


def selector_parts(value: str, *, with_atom: bool = False) -> tuple[Any, ...]:
    pattern = (
        r"([^:]+):([A-Za-z0-9_]+?)(-?[0-9]+):([A-Za-z0-9_]+)"
        if with_atom
        else r"([^:]+):([A-Za-z0-9_]+?)(-?[0-9]+)"
    )
    match = re.fullmatch(pattern, value)
    if match is None:
        raise PreparationError(f"invalid selector {value!r}")
    if with_atom:
        return match.group(1), match.group(2).upper(), int(match.group(3)), match.group(4).upper()
    return match.group(1), match.group(2).upper(), int(match.group(3))


def insertion_code(residue: gemmi.Residue) -> str:
    value = str(residue.seqid.icode).strip()
    return "" if value in {"", "\x00", ".", "?"} else value


def altloc(atom: gemmi.Atom) -> str:
    value = str(atom.altloc).strip()
    return "" if value in {"", "\x00", ".", "?"} else value


def heavy_atom_map(
    path: Path, row: Mapping[str, str], *, normalize_source_names: bool
) -> dict[tuple[str, int, str, str, str, str], tuple[float, float, float]]:
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise PreparationError(f"{path} does not contain exactly one model")
    metal_chain, _, metal_num, _ = selector_parts(
        row["metal_selector_source"], with_atom=True
    )
    pqq_chain, _, pqq_num = selector_parts(row["pqq_selector_source"])
    result: dict[tuple[str, int, str, str, str, str], tuple[float, float, float]] = {}
    for chain in structure[0]:
        for residue in chain:
            resname = residue.name.upper().strip()
            if normalize_source_names:
                if chain.name == metal_chain and residue.seqid.num == metal_num:
                    resname = "LA"
                elif chain.name == pqq_chain and residue.seqid.num == pqq_num:
                    resname = "PQQ"
            for atom in residue:
                element = atom.element.name.upper()
                if element in {"H", "D"}:
                    continue
                if altloc(atom):
                    raise PreparationError(f"unexpected alternate conformer in {path}")
                key = (
                    chain.name,
                    residue.seqid.num,
                    insertion_code(residue),
                    resname,
                    atom.name.upper().strip(),
                    element,
                )
                if key in result:
                    raise PreparationError(f"duplicate heavy-atom identity {key} in {path}")
                result[key] = (atom.pos.x, atom.pos.y, atom.pos.z)
    return result


def coordinate_digest(
    atoms: Mapping[tuple[str, int, str, str, str, str], tuple[float, float, float]]
) -> str:
    serial = [
        [*identity, *[round(float(value), 6) for value in atoms[identity]]]
        for identity in sorted(atoms)
    ]
    return hashlib.sha256(
        json.dumps(serial, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def compare_heavy_atoms(
    reference: Mapping[tuple[str, int, str, str, str, str], tuple[float, float, float]],
    query: Mapping[tuple[str, int, str, str, str, str], tuple[float, float, float]],
    *,
    label: str,
) -> dict[str, Any]:
    if set(reference) != set(query):
        missing = sorted(set(reference) - set(query))
        added = sorted(set(query) - set(reference))
        raise PreparationError(
            f"{label} changed heavy-atom identities: missing={missing[:5]}, added={added[:5]}"
        )
    displacements = {
        key: math.dist(reference[key], query[key]) for key in reference
    }
    worst_key, maximum = max(displacements.items(), key=lambda item: item[1])
    if maximum > PDB_SERIALIZATION_TOLERANCE_A:
        raise PreparationError(
            f"{label} moved {worst_key} by {maximum:.8f} A, beyond "
            f"{PDB_SERIALIZATION_TOLERANCE_A:.8f} A"
        )
    return {
        "label": label,
        "atom_count": len(reference),
        "identity_sets_equal": True,
        "reference_coordinate_digest": coordinate_digest(reference),
        "query_coordinate_digest": coordinate_digest(query),
        "maximum_displacement_A": maximum,
        "maximum_displacement_atom": list(worst_key),
        "allowed_PDB_serialization_tolerance_A": PDB_SERIALIZATION_TOLERANCE_A,
        "passes": True,
    }


def run_checked(command: Sequence[str], *, label: str) -> None:
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise PreparationError(f"{label} failed with exit code {completed.returncode}")


def parse_xyz(path: Path) -> tuple[str, list[str]]:
    lines = path.read_text().splitlines()
    if len(lines) < 3 or not lines[0].isdigit() or int(lines[0]) != len(lines) - 2:
        raise PreparationError(f"malformed XYZ file: {path}")
    return lines[1], lines[2:]


def validate_prepared_pair(
    manifest_path: Path,
    row: Mapping[str, str],
    pins: Mapping[str, Any],
    pins_path: Path,
) -> dict[str, Any]:
    manifest = read_object(manifest_path)
    if manifest.get("protocol_id") != PROTOCOL_ID or manifest.get("panel_id") != row["id"]:
        raise PreparationError(f"wrong protocol/panel identity in {manifest_path}")
    carver_record = pins["experiment_helpers"]["fixed_core_carver"]
    if (
        manifest.get("carver", {}).get("path") != carver_record["path"]
        or manifest.get("carver", {}).get("implementation_sha256")
        != carver_record["sha256"]
        or manifest.get("experiment_provenance", {}).get("implementation_pins")
        != {"path": str(pins_path.resolve()), "sha256": sha256_file(pins_path)}
    ):
        raise PreparationError(f"{row['id']} carve implementation provenance fails")
    tasks = manifest.get("tasks")
    if not isinstance(tasks, list) or [task.get("task_id") for task in tasks] != ["La", "Ca"]:
        raise PreparationError(f"{row['id']} does not have exactly La and Ca tasks")
    acidic = row["include_plus2_as_ln_specific_asp"].lower() == "true"
    expected_charges = {"La": -2, "Ca": -3} if acidic else {"La": -1, "Ca": -2}
    ledger = manifest.get("charge_ledger")
    if not isinstance(ledger, dict) or {
        "La": ledger.get("La_total"), "Ca": ledger.get("Ca_total")
    } != expected_charges or ledger.get("expected_total_charges") != expected_charges:
        raise PreparationError(f"{row['id']} has the wrong frozen charge ledger")
    fragments = manifest.get("qm_fragments")
    if not isinstance(fragments, list) or len(fragments) != (6 if acidic else 5):
        raise PreparationError(f"{row['id']} has the wrong fixed fragment count")
    for fragment in fragments:
        records = fragment.get("atom_records") if isinstance(fragment, dict) else None
        if not isinstance(records, list) or len(records) != fragment.get("atom_count"):
            raise PreparationError(f"{row['id']} lacks complete ordered atom provenance")
        for record in records:
            if not isinstance(record, dict) or not isinstance(record.get("xyz_A"), list):
                raise PreparationError(f"{row['id']} has malformed atom provenance")
    coordination = manifest.get("coordination")
    if (
        not isinstance(coordination, dict)
        or coordination.get("coordination_number", 0) < 6
        or "source_typed_direct_donors_ordered_by_distance" not in coordination
        or manifest.get("fixed_core", {}).get("unexpected_direct_donors") != []
    ):
        raise PreparationError(f"{row['id']} fails the frozen typed-CN contract")
    fixed_core = manifest.get("fixed_core")
    if (
        not isinstance(fixed_core, dict)
        or fixed_core.get("water_policy") != "dry_exclude_all_source_and_synthetic_waters"
        or fixed_core.get("point_charge_embedding") is not False
    ):
        raise PreparationError(f"{row['id']} changed the dry/no-embedding protocol")

    root = manifest_path.parent
    xyz_paths: dict[str, Path] = {}
    for task in tasks:
        task_id = task["task_id"]
        for key in ("input", "xyz"):
            record = task.get(key)
            if not isinstance(record, dict):
                raise PreparationError(f"{row['id']} task {task_id} lacks {key}")
            artifact = (root / str(record.get("path", ""))).resolve()
            if artifact.parent != root or sha256_file(artifact) != record.get("sha256"):
                raise PreparationError(f"{row['id']} task {task_id} {key} provenance fails")
            if key == "input":
                text = artifact.read_text()
                if EXPECTED_ORCA_DIRECTIVE not in text or "%pal" in text.lower():
                    raise PreparationError(f"{row['id']} task {task_id} ORCA template changed")
            else:
                xyz_paths[task_id] = artifact
    _, la_atoms = parse_xyz(xyz_paths["La"])
    _, ca_atoms = parse_xyz(xyz_paths["Ca"])
    if len(la_atoms) != len(ca_atoms) or not la_atoms[0].split()[0].upper() == "LA":
        raise PreparationError(f"{row['id']} La/Ca XYZ atom counts or La identity fail")
    if ca_atoms[0].split()[0].upper() != "CA" or la_atoms[1:] != ca_atoms[1:]:
        raise PreparationError(
            f"{row['id']} La/Ca nonmetal coordinate payloads are not byte-identical"
        )
    return {
        "manifest": file_record(manifest_path),
        "coordination_number": coordination["coordination_number"],
        "direct_nitrogen_count": coordination["direct_nitrogen_count"],
        "fragment_count": len(fragments),
        "atom_count_per_arm": len(la_atoms),
        "charges": expected_charges,
        "nonmetal_coordinates_byte_identical": True,
        "tasks": {task["task_id"]: task for task in tasks},
        "carver_sha256": manifest.get("carver", {}).get("implementation_sha256"),
        "runner_sha256": pins["canonical_helpers"]["run_orca_task_manifest"]["sha256"],
    }


def prepare(out_root: Path, pins_path: Path) -> Path:
    pins_path = pins_path.resolve()
    pins = load_and_verify_pins(pins_path)
    protonation_subprotocol = pins["protonation_subprotocol"]
    core_map_path = verify_file_record(pins["core_map"], "frozen core map")
    rows = load_core_map(core_map_path)
    out_root = out_root.resolve()
    out_root.mkdir(parents=True, exist_ok=False)
    python = str(Path(pins["preparation_runtime"]["python"]["path"]).resolve())
    normalizer = str(verify_file_record(
        pins["canonical_helpers"]["normalize_af3_cif"], "normalizer"
    ))
    protonator = str(verify_file_record(
        pins["experiment_helpers"]["protonate_standard_only"],
        "standard-only protonation wrapper",
    ))
    carver = str(verify_file_record(
        pins["experiment_helpers"]["fixed_core_carver"], "fixed-core carver"
    ))
    records: list[dict[str, Any]] = []

    for index, row in enumerate(rows, start=1):
        target_id = row["id"]
        safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", target_id)
        stem = f"pmdh_fc_{safe_id}"
        target_dir = out_root / f"{index:02d}_{safe_id}"
        target_dir.mkdir()
        source = Path(row["source_cif"]).resolve()
        if sha256_file(source) != row["source_sha256"]:
            raise PreparationError(f"source CIF hash mismatch for {target_id}")
        normalized = target_dir / f"{stem}_normalized.pdb"
        protonated = target_dir / f"{stem}_protonated.pdb"
        run_checked([python, normalizer, str(source), str(normalized)], label=f"normalize {target_id}")

        source_atoms = heavy_atom_map(source, row, normalize_source_names=True)
        normalized_atoms = heavy_atom_map(normalized, row, normalize_source_names=False)
        normalized_residues = {identity[:4] for identity in normalized_atoms}
        if sum(identity[3] == "LA" for identity in normalized_residues) != 1 or sum(
            identity[3] == "PQQ" for identity in normalized_residues
        ) != 1:
            raise PreparationError(
                f"normalization did not produce exactly one LA and one PQQ for {target_id}"
            )
        source_to_normalized = compare_heavy_atoms(
            source_atoms, normalized_atoms, label="source_CIF_to_normalized_PDB"
        )
        normalization_manifest = target_dir / f"{stem}_normalization_manifest.json"
        normalization_record = {
            "schema_version": NORMALIZATION_SCHEMA,
            "protocol_id": PROTOCOL_ID,
            "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "source": {"path": str(source), "sha256": row["source_sha256"]},
            "output": file_record(normalized),
            "normalizer": pins["canonical_helpers"]["normalize_af3_cif"],
            "renamed_residues": {"lanthanum": 1, "pqq": 1},
            "source_heavy_coordinate_check": source_to_normalized,
        }
        write_json_atomic(normalization_manifest, normalization_record)

        run_checked(
            [
                python,
                protonator,
                str(normalized),
                str(protonated),
                "--ph", "7.0",
                "--implementation-pins", str(pins_path),
            ],
            label=f"protonate {target_id}",
        )
        protonation_manifest = target_dir / f"{stem}_protonation_manifest.json"
        protonation = read_object(protonation_manifest)
        nonstandard_policy = protonation.get("nonstandard_definition_policy")
        if (
            protonation.get("protocol_id") != "pdbfixer_standard_residue_protonation_v2"
            or protonation.get("experiment_protonation_protocol_id")
            != protonation_subprotocol["id"]
            or protonation.get("canonical_protonator")
            != pins["canonical_helpers"]["protonate_cif"]
            or protonation.get("experiment_wrapper")
            != pins["experiment_helpers"]["protonate_standard_only"]
            or protonation.get("ph") != 7.0
            or protonation.get("add_missing_residues") is not False
            or protonation.get("repaired_missing_atom_count") != 0
            or protonation.get("repaired_missing_terminal_atom_count") != 0
            or protonation.get("cofactor_protonation") != "not_assigned"
            or protonation.get("software") != pins["preparation_runtime"]["packages"]
            or not isinstance(nonstandard_policy, dict)
            or nonstandard_policy.get("policy_id")
            != protonation_subprotocol["nonstandard_definition_policy_id"]
            or nonstandard_policy.get("observed_nonstandard_residue_counts")
            != protonation_subprotocol["allowed_nonstandard_residue_counts"]
            or nonstandard_policy.get("output_nonstandard_hydrogen_counts")
            != protonation_subprotocol[
                "required_output_nonstandard_hydrogen_counts"
            ]
            or nonstandard_policy.get("forcefield")
            != protonation_subprotocol["forcefield"]
            or nonstandard_policy.get("python_random_seed")
            != protonation_subprotocol["python_random_seed"]
            or nonstandard_policy.get("openmm_platform")
            != protonation_subprotocol["openmm_platform"]
            or nonstandard_policy.get("openmm_CPU_Threads")
            != protonation_subprotocol["openmm_cpu_threads"]
            or nonstandard_policy.get("scientific_chemistry_changed") is not False
        ):
            raise PreparationError(f"protonation contract failed for {target_id}")
        protonated_atoms = heavy_atom_map(protonated, row, normalize_source_names=False)
        normalized_to_protonated = compare_heavy_atoms(
            normalized_atoms, protonated_atoms, label="normalized_to_protonated_PDB"
        )
        source_to_protonated = compare_heavy_atoms(
            source_atoms, protonated_atoms, label="source_CIF_to_protonated_PDB"
        )
        heavy_check_path = target_dir / f"{stem}_heavy_coordinate_check.json"
        heavy_check = {
            "schema_version": HEAVY_CHECK_SCHEMA,
            "protocol_id": PROTOCOL_ID,
            "source": {"path": str(source), "sha256": row["source_sha256"]},
            "normalized": file_record(normalized),
            "protonated": file_record(protonated),
            "checks": [source_to_normalized, normalized_to_protonated, source_to_protonated],
            "all_source_heavy_atoms_retained": True,
            "no_heavy_atoms_added": True,
            "coordinate_changes_limited_to_PDB_serialization_tolerance": True,
            "passes": True,
        }
        write_json_atomic(heavy_check_path, heavy_check)

        run_checked(
            [
                python,
                carver,
                str(protonated),
                str(target_dir),
                "--stem", stem,
                "--core-map", str(core_map_path),
                "--target-id", target_id,
                "--implementation-pins", str(pins_path),
            ],
            label=f"carve {target_id}",
        )
        carve_manifest = target_dir / f"{stem}_carve_manifest.json"
        pair_record = validate_prepared_pair(carve_manifest, row, pins, pins_path)
        records.append(
            {
                "panel_index": index,
                "panel_id": target_id,
                "class": row["class"],
                "acidic_Dplus2": row["include_plus2_as_ln_specific_asp"].lower() == "true",
                "source_cif": {"path": str(source), "sha256": row["source_sha256"]},
                "normalized": file_record(normalized),
                "normalization_manifest": file_record(normalization_manifest),
                "protonated": file_record(protonated),
                "protonation_manifest": file_record(protonation_manifest),
                "heavy_coordinate_check": file_record(heavy_check_path),
                **pair_record,
            }
        )

    preparation_path = out_root / "preparation.json"
    preparation = {
        "schema_version": PREPARATION_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "ready_for_orca",
        "implementation_pins": file_record(pins_path),
        "preregistration": pins["preregistration"],
        "core_map": pins["core_map"],
        "reserved_holdouts": {
            "spec": pins["reserved_holdout_spec"],
            "audit": pins["reserved_holdout_audit"],
            "consumed": False,
        },
        "preparation_runtime": pins["preparation_runtime"],
        "protonation_subprotocol": {
            **protonation_subprotocol,
            "wrapper": pins["experiment_helpers"]["protonate_standard_only"],
            "amendment": pins["preparation_amendment"],
        },
        "panel_counts": EXPECTED_COUNTS,
        "target_count": len(records),
        "task_count": 2 * len(records),
        "all_source_heavy_coordinates_preserved": True,
        "all_nonmetal_arm_coordinates_byte_identical": True,
        "water_policy": "dry_no_source_or_synthetic_water",
        "point_charge_embedding": False,
        "targets": records,
    }
    write_json_atomic(preparation_path, preparation)
    return preparation_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "prepared",
        help="fresh output directory (must not already exist)",
    )
    parser.add_argument(
        "--implementation-pins", type=Path, default=HERE / "implementation_pins.json"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = prepare(args.output, args.implementation_pins)
    except (PreparationError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
