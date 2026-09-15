#!/usr/bin/env python3
"""Prepare the preregistered Q46444/1KB0 fixed_core_v3 external test.

This wrapper supplies a separately pinned selector authority to the immutable
generic crystal-transfer implementation.  It prepares two ORCA inputs and does
not execute ORCA.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import importlib.util
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
CALIBRATION_DIR = PROJECT / "diagnostics/pqq_pmdh_fixed_core_calibration_20260914"
GENERIC_PREPARER = CALIBRATION_DIR / "reserved_crystal_holdout/prepare_holdouts.py"
DEFAULT_PINS = HERE / "implementation_pins.json"
DEFAULT_OUTPUT = HERE / "prepared"
PROTOCOL_ID = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
PINS_SCHEMA = "alchemical_bvs.pqq_q46444_1kb0_external_pins.v1"
PREPARATION_SCHEMA = "alchemical_bvs.pqq_q46444_1kb0_external_preparation.v1"


def _load_generic_preparer():
    spec = importlib.util.spec_from_file_location(
        "pqq_fixed_core_v3_generic_crystal_preparer", GENERIC_PREPARER
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load immutable generic preparer {GENERIC_PREPARER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


holdout = _load_generic_preparer()


class ExternalPreparationError(RuntimeError):
    """The external selector, pins, or prepared payload is invalid."""


def read_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ExternalPreparationError(f"cannot read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ExternalPreparationError(f"JSON root is not an object: {path}")
    return value


def write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def verify_record(record: Any, label: str) -> Path:
    try:
        return holdout.verify_file_record(record, label)
    except holdout.HoldoutPreparationError as exc:
        raise ExternalPreparationError(str(exc)) from exc


def verify_external_pins(
    pins_path: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, str], dict[str, Any]]:
    pins_path = pins_path.resolve()
    pins = read_object(pins_path)
    if pins.get("schema_version") != PINS_SCHEMA or pins.get("protocol_id") != PROTOCOL_ID:
        raise ExternalPreparationError("external pins schema or protocol changed")

    calibration_pins_path = verify_record(
        pins.get("calibration_implementation_pins"), "calibration implementation pins"
    )
    calibration_pins = holdout.fixed.verify_pins(calibration_pins_path)
    calibration_result = verify_record(
        pins.get("calibration_result"), "locked calibration result"
    )
    for key, label in (
        ("preregistration", "external preregistration"),
        ("target_spec", "external target spec"),
        ("selector_audit", "external selector audit"),
    ):
        verify_record(pins.get(key), label)

    implementation = pins.get("implementation")
    if not isinstance(implementation, dict) or set(implementation) != {
        "prepare_external",
        "run_external",
        "run_sbatch",
    }:
        raise ExternalPreparationError("external implementation file set changed")
    verified_implementation = {
        name: verify_record(record, f"external implementation {name}")
        for name, record in implementation.items()
    }
    if verified_implementation["prepare_external"] != Path(__file__).resolve():
        raise ExternalPreparationError("executing preparer differs from external pins")

    inherited = pins.get("inherited")
    if not isinstance(inherited, dict) or set(inherited) != {
        "generic_crystal_preparer",
        "fixed_core_carver",
        "protonate_standard_only",
        "run_orca_task_manifest",
        "render_orca_runtime_input",
    }:
        raise ExternalPreparationError("inherited implementation set changed")
    for name, record in inherited.items():
        verify_record(record, f"inherited implementation {name}")
    expected_inherited = {
        "fixed_core_carver": calibration_pins["experiment_helpers"]["fixed_core_carver"],
        "protonate_standard_only": calibration_pins["experiment_helpers"]["protonate_standard_only"],
        "run_orca_task_manifest": calibration_pins["canonical_helpers"]["run_orca_task_manifest"],
        "render_orca_runtime_input": calibration_pins["canonical_helpers"]["render_orca_runtime_input"],
    }
    for name, expected in expected_inherited.items():
        if inherited.get(name) != expected:
            raise ExternalPreparationError(f"inherited {name} differs from calibration pins")
    if Path(inherited["generic_crystal_preparer"]["path"]).resolve() != Path(
        holdout.__file__
    ).resolve():
        raise ExternalPreparationError("imported generic crystal preparer differs from pins")

    expected_policy = {
        "model": "1",
        "chain": "A",
        "alternate_conformer_policy": holdout.ALTLOC_POLICY_ID,
        "water_policy": "exclude_all_source_and_synthetic_waters",
        "noncore_heterogen_policy": "exclude_all_without_replacement_including_TFB1810",
        "geometry_relaxation": False,
        "point_charge_embedding": False,
    }
    if pins.get("selection_policy") != expected_policy:
        raise ExternalPreparationError("external selection policy changed")
    expected_rule = {
        "expected_band": "Ca-supported",
        "pass_if_S_kcal_mol_at_most": 14.857129202922806,
        "indeterminate_open_interval": [14.857129202922806, 23.460061205609236],
        "fail_if_S_kcal_mol_at_least": 23.460061205609236,
        "invalid_pair": "failure",
        "threshold_refit_allowed": False,
    }
    if pins.get("frozen_decision_rule") != expected_rule:
        raise ExternalPreparationError("external frozen decision rule changed")

    target_path = verify_record(pins["target_spec"], "external target spec")
    with target_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ExternalPreparationError("external target spec must contain exactly one row")
    row = rows[0]
    exact = {
        "holdout_role": "external_primary",
        "pdb_id": "1KB0",
        "biological_class": "Ca/QH-ADH",
        "model": "1",
        "selected_chain": "A",
        "metal_selector": "A:CA801/CA",
        "pqq_selector": "A:PQQ1800",
        "fixed_typed_CN_3p1": "7",
        "glu_selector": "A:GLU185",
        "asn_selector": "A:ASN263",
        "catalytic_asp_selector": "A:ASP308",
        "plus2_selector": "A:THR310",
        "include_plus2": "false",
        "cationic_partner_selector": "A:LYS335",
        "source_water_count": "1016",
        "waters_within_3p6_A": "0",
        "noncore_direct_ligand": "none",
    }
    for field, expected in exact.items():
        if row.get(field) != expected:
            raise ExternalPreparationError(
                f"external target field {field} changed: {row.get(field)!r}"
            )
    if "TFB1810 OXT at 3.338 A" not in row.get("dry_policy", ""):
        raise ExternalPreparationError("external dry policy no longer names TFB1810/OXT")
    holdout.verify_runtime(calibration_pins)
    gate = holdout.verify_calibration_gate(
        calibration_result, calibration_pins_path, calibration_pins
    )
    return pins, calibration_pins, row, gate


def prepare_external(*, output: Path, pins_path: Path) -> Path:
    pins_path = pins_path.resolve()
    pins, calibration_pins, row, gate = verify_external_pins(pins_path)
    calibration_pins_path = Path(
        pins["calibration_implementation_pins"]["path"]
    ).resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    incomplete = output / "EXTERNAL_PREPARATION_INCOMPLETE"
    holdout.write_text_new(
        incomplete,
        "This directory is incomplete and must not be submitted to ORCA.\n",
    )

    inherited_context = {
        "holdout_spec": pins["target_spec"],
        "holdout_selector_audit": pins["selector_audit"],
        "holdout_implementation": {
            "prepare_holdouts": pins["inherited"]["generic_crystal_preparer"]
        },
    }
    target = holdout.prepare_one(
        row,
        1,
        output,
        gate,
        pins_path,
        inherited_context,
        calibration_pins_path,
        calibration_pins,
    )
    manifest_path = verify_record(target["manifest"], "prepared 1KB0 carve manifest")
    manifest = read_object(manifest_path)
    normalized_path = verify_record(target["normalization_manifest"], "normalization manifest")
    normalization = read_object(normalized_path)
    excluded = normalization.get("selection", {}).get("excluded_selected_chain_residues", [])
    excluded_selectors = {
        item.get("selector") for item in excluded if isinstance(item, dict)
    }
    if "A:TFB1810" not in excluded_selectors:
        raise ExternalPreparationError("normalization did not explicitly exclude A:TFB1810")
    if (
        manifest.get("panel_id") != "1KB0"
        or manifest.get("coordination", {}).get("coordination_number") != 7
        or manifest.get("charge_ledger", {}).get("expected_total_charges")
        != {"La": -1, "Ca": -2}
        or manifest.get("fixed_core", {}).get("synthetic_water_count") != 0
        or manifest.get("fixed_core", {}).get("replacement_ligand_count") != 0
        or manifest.get("fixed_core", {}).get("point_charge_embedding") is not False
        or manifest.get("fixed_core", {}).get("geometry_relaxation") is not False
        or manifest.get("paired_arm_invariant", {}).get(
            "nonmetal_coordinates_byte_identical"
        )
        is not True
        or manifest.get("fixed_core", {}).get("excluded_noncore_direct_ligands") != []
    ):
        raise ExternalPreparationError("prepared 1KB0 manifest violates preregistration")
    if list(output.rglob("*.out")):
        raise ExternalPreparationError("preparation unexpectedly created ORCA output")

    preparation = {
        "schema_version": PREPARATION_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "ready_for_orca_after_passing_calibration_gate",
        "experiment_id": "Q46444_1KB0_external_Ca_structural_transfer",
        "expected_band": "Ca-supported",
        "calibration_gate": gate,
        "external_pins": holdout.file_record(pins_path),
        "wrapper": pins["implementation"]["prepare_external"],
        "inherited_generic_crystal_preparer": pins["inherited"][
            "generic_crystal_preparer"
        ],
        "selector_authority": pins["target_spec"],
        "selector_audit": pins["selector_audit"],
        "target_count": 1,
        "task_count": 2,
        "water_policy": "dry_exclude_all_source_and_synthetic_waters",
        "noncore_heterogen_policy": "exclude_all_without_replacement_including_TFB1810",
        "TFB1810_OXT_source_distance_A": 3.338,
        "all_retained_source_heavy_coordinates_preserved": True,
        "all_nonmetal_arm_coordinates_byte_identical": True,
        "orca_executed": False,
        "target": target,
    }
    preparation_path = output / "external_preparation.json"
    write_json_atomic(preparation_path, preparation)
    incomplete.unlink()
    return preparation_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--pins", type=Path, default=DEFAULT_PINS)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = prepare_external(output=args.output, pins_path=args.pins)
    except (
        ExternalPreparationError,
        holdout.HoldoutPreparationError,
        holdout.fixed.FixedCoreError,
        holdout.base.PQQCarveError,
        OSError,
        ValueError,
    ) as exc:
        output = args.output.resolve()
        if output.is_dir() and (output / "EXTERNAL_PREPARATION_INCOMPLETE").is_file():
            write_json_atomic(
                output / "EXTERNAL_PREPARATION_ERROR.json",
                {
                    "schema_version": "alchemical_bvs.pqq_q46444_1kb0_external_preparation_error.v1",
                    "protocol_id": PROTOCOL_ID,
                    "error_type": type(exc).__name__,
                    "message": str(exc),
                    "valid_for_orca": False,
                },
            )
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

