#!/usr/bin/env python3
"""Apply only the preregistered fixed-core PQQ-MDH calibration rule."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
PROTOCOL_ID = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
RESULT_SCHEMA = "alchemical_bvs.pqq_fixed_core_calibration_result.v1"
EXECUTION_SCHEMA = "alchemical_bvs.orca_task_execution.v1"
HA2KCAL = 627.509474
MIN_GAP_KCAL = 5.0
ENERGY_RE = re.compile(
    r"FINAL SINGLE POINT ENERGY\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)


class ScoreError(RuntimeError):
    """The complete, frozen calibration panel cannot be scored safely."""


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
        raise ScoreError(f"cannot read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ScoreError(f"JSON root is not an object: {path}")
    return value


def verify_file_record(record: Any, label: str) -> Path:
    if not isinstance(record, dict):
        raise ScoreError(f"{label} is not a file record")
    path = Path(str(record.get("path", ""))).resolve()
    expected = record.get("sha256")
    if not path.is_file() or not isinstance(expected, str):
        raise ScoreError(f"{label} file/hash is missing")
    if sha256_file(path) != expected:
        raise ScoreError(f"{label} SHA-256 mismatch")
    return path


def write_new(path: Path, text: str) -> None:
    with path.open("x") as handle:
        handle.write(text)


def load_core_map(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    mapping = {row["id"]: row for row in rows}
    if len(rows) != 25 or len(mapping) != 25:
        raise ScoreError("frozen core map is not the unique 25-member panel")
    if sum(row["class"] == "La" for row in rows) != 11 or sum(
        row["class"] == "Ca" for row in rows
    ) != 14:
        raise ScoreError("frozen class counts changed")
    return mapping


def resolve_under(root: Path, value: str, label: str) -> Path:
    path = (root / value).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ScoreError(f"{label} escapes prepared target directory") from exc
    return path


def parse_energy(output_path: Path) -> float:
    text = output_path.read_text(errors="strict")
    values = ENERGY_RE.findall(text)
    if len(values) != 1:
        raise ScoreError(
            f"{output_path} contains {len(values)} FINAL SINGLE POINT ENERGY records"
        )
    if "ORCA TERMINATED NORMALLY" not in text:
        raise ScoreError(f"{output_path} lacks normal ORCA termination")
    if "SCF CONVERGED AFTER" not in text or "SCF NOT CONVERGED" in text:
        raise ScoreError(f"{output_path} lacks unambiguous SCF convergence")
    return float(values[0])


def score_task(
    *,
    task_id: str,
    task: Mapping[str, Any],
    manifest_path: Path,
    manifest_hash: str,
    pins: Mapping[str, Any],
) -> tuple[float, dict[str, Any]]:
    root = manifest_path.parent
    input_record = task.get("input")
    xyz_record = task.get("xyz")
    if not isinstance(input_record, dict) or not isinstance(xyz_record, dict):
        raise ScoreError(f"{manifest_path} task {task_id} lacks input/XYZ records")
    input_path = resolve_under(root, str(input_record.get("path", "")), f"{task_id} input")
    xyz_path = resolve_under(root, str(xyz_record.get("path", "")), f"{task_id} XYZ")
    output_path = resolve_under(root, str(task.get("output_path", "")), f"{task_id} output")
    for artifact, record, label in (
        (input_path, input_record, "input"),
        (xyz_path, xyz_record, "XYZ"),
    ):
        if not artifact.is_file() or sha256_file(artifact) != record.get("sha256"):
            raise ScoreError(f"{manifest_path} task {task_id} {label} hash mismatch")
    if not output_path.is_file():
        raise ScoreError(f"missing ORCA output {output_path}")
    execution_path = output_path.with_name(output_path.name + ".execution.json")
    execution = read_object(execution_path)
    runner = pins["canonical_helpers"]["run_orca_task_manifest"]
    renderer = pins["canonical_helpers"]["render_orca_runtime_input"]
    orca = pins["orca_runtime"]["executable"]
    valid = bool(
        execution.get("schema_version") == EXECUTION_SCHEMA
        and execution.get("task_id") == task_id
        and execution.get("returncode") == 0
        and execution.get("normal_termination") is True
        and execution.get("scf_converged") is True
        and execution.get("orca_version") == pins["orca_runtime"]["version"]
        and execution.get("manifest") == {
            "path": str(manifest_path), "sha256": manifest_hash
        }
        and execution.get("task_runner") == runner
        and execution.get("runtime_renderer") == renderer
        and execution.get("orca_executable") == orca
        and execution.get("artifacts", {}).get("template_input", {}).get("sha256")
        == input_record.get("sha256")
        and execution.get("artifacts", {}).get("xyz", {}).get("sha256")
        == xyz_record.get("sha256")
        and execution.get("artifacts", {}).get("output", {}).get("sha256")
        == sha256_file(output_path)
    )
    if not valid:
        raise ScoreError(f"execution provenance failed for {output_path}")
    energy = parse_energy(output_path)
    return energy, {
        "energy_hartree": energy,
        "input": {"path": str(input_path), "sha256": sha256_file(input_path)},
        "xyz": {"path": str(xyz_path), "sha256": sha256_file(xyz_path)},
        "output": {"path": str(output_path), "sha256": sha256_file(output_path)},
        "execution": {"path": str(execution_path), "sha256": sha256_file(execution_path)},
    }


def auc_larger_is_la(rows: Sequence[Mapping[str, Any]]) -> float:
    la = [float(row["R_kcal_mol"]) for row in rows if row["class"] == "La"]
    ca = [float(row["R_kcal_mol"]) for row in rows if row["class"] == "Ca"]
    wins = sum(1.0 if left > right else 0.5 if left == right else 0.0 for left in la for right in ca)
    return wins / (len(la) * len(ca))


def leave_one_out(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for held in rows:
        training = [row for row in rows if row["panel_id"] != held["panel_id"]]
        ca_scores = [float(row["R_kcal_mol"]) for row in training if row["class"] == "Ca"]
        la_scores = [float(row["R_kcal_mol"]) for row in training if row["class"] == "La"]
        upper_ca = max(ca_scores)
        lower_la = min(la_scores)
        training_strict = lower_la > upper_ca
        threshold = (upper_ca + lower_la) / 2.0 if training_strict else None
        value = float(held["R_kcal_mol"])
        if threshold is None or value == threshold:
            predicted = "failure"
            correct = False
        else:
            predicted = "La" if value > threshold else "Ca"
            correct = predicted == held["class"]
        results.append(
            {
                "held_out_panel_id": held["panel_id"],
                "true_class": held["class"],
                "R_kcal_mol": value,
                "training_U_max_Ca_kcal_mol": upper_ca,
                "training_L_min_La_kcal_mol": lower_la,
                "training_strict_separation": training_strict,
                "training_midpoint_kcal_mol": threshold,
                "prediction": predicted,
                "correct": correct,
            }
        )
    return results


def render_markdown(result: Mapping[str, Any]) -> str:
    calibration = result["calibration"]
    lines = [
        "# Fixed-core PQQ-MDH calibration result",
        "",
        f"**Protocol:** `{PROTOCOL_ID}`  ",
        f"**Verdict:** **{calibration['verdict']}**",
        "",
        "## Scientific interpretation",
        "",
        result["interpretation"]["incremental_discrimination_statement"],
        "",
        "## Frozen gates",
        "",
        "| Gate | Result |",
        "|---|---:|",
        f"| Valid unchanged pairs | {calibration['valid_pair_count']}/25 |",
        f"| Strict separation (L > U) | {calibration['strict_separation']} |",
        f"| Gap G >= 5.0 kcal/mol | {calibration['gap_floor_pass']} |",
        f"| Strict leave-one-out midpoint | {calibration['loo_correct_count']}/25 |",
        f"| AUROC | {calibration['auroc']:.6f} |",
        "",
        f"U (largest Ca R): {calibration['U_max_Ca_kcal_mol']:.6f} kcal/mol  ",
        f"L (smallest La R): {calibration['L_min_La_kcal_mol']:.6f} kcal/mol  ",
        f"G = L-U: {calibration['gap_kcal_mol']:.6f} kcal/mol  ",
        f"Released midpoint T_R (only if every gate passes): {calibration['threshold_R_kcal_mol']}",
        "",
        f"Centered U_S (largest Ca S): {calibration['U_max_Ca_S_kcal_mol']:.6f} kcal/mol  ",
        f"Centered L_S (smallest La S): {calibration['L_min_La_S_kcal_mol']:.6f} kcal/mol  ",
        f"Released centered midpoint T_S: {calibration['threshold_S_kcal_mol']}  ",
        "The gap and every classification decision are identical in R and S; S is only a constant-shifted display scale.",
        "",
        "## Scores",
        "",
        "| Panel ID | Class | D+2 Asp | E(La), Eh | E(Ca), Eh | R, kcal/mol | S (aquo gauge), kcal/mol |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in result["scores"]:
        lines.append(
            f"| {row['panel_id']} | {row['class']} | {row['acidic_Dplus2']} | "
            f"{row['energies_hartree']['La']:.12f} | {row['energies_hartree']['Ca']:.12f} | "
            f"{row['R_kcal_mol']:.6f} | {row['S_aquo_gauge_kcal_mol']:.6f} |"
        )
    lines.extend(
        [
            "",
            "R is the primary preregistered contrast. S differs only by the pinned constant aquo gauge; it was not used for any gate.",
            "",
        ]
    )
    return "\n".join(lines)


def score(preparation_path: Path, pins_path: Path, output_dir: Path) -> tuple[Path, Path]:
    preparation_path = preparation_path.resolve()
    pins_path = pins_path.resolve()
    output_dir = output_dir.resolve()
    preparation = read_object(preparation_path)
    pins = read_object(pins_path)
    if pins.get("protocol_id") != PROTOCOL_ID:
        raise ScoreError("implementation pins carry another protocol ID")
    for group in ("canonical_helpers", "experiment_helpers"):
        records = pins.get(group)
        if not isinstance(records, dict):
            raise ScoreError(f"implementation pins lack {group}")
        for name, record in records.items():
            verify_file_record(record, f"{group}.{name}")
    own_record = pins["experiment_helpers"].get("score")
    if (
        not isinstance(own_record, dict)
        or Path(str(own_record.get("path", ""))).resolve() != Path(__file__).resolve()
    ):
        raise ScoreError("scorer is not self-bound in implementation pins")
    preregistration = verify_file_record(pins.get("preregistration"), "preregistration")
    core_map_path = verify_file_record(pins.get("core_map"), "core map")
    verify_file_record(pins.get("reserved_holdout_spec"), "reserved holdout spec")
    verify_file_record(pins.get("reserved_holdout_audit"), "reserved holdout audit")
    verify_file_record(pins.get("preparation_amendment"), "preparation amendment")
    aquo_path = verify_file_record(pins.get("aquo_reference"), "aquo reference")
    aquo_payload = read_object(aquo_path)
    if (
        aquo_payload.get("reference_id") != pins["aquo_reference"].get("reference_id")
        or aquo_payload.get("delta_E_aquo_hartree")
        != pins["aquo_reference"].get("delta_E_aquo_hartree")
    ):
        raise ScoreError("pinned aquo identity or energy difference changed")
    verify_file_record(pins.get("orca_runtime", {}).get("executable"), "ORCA executable")
    if (
        preparation.get("schema_version") != "alchemical_bvs.pqq_fixed_core_preparation.v1"
        or preparation.get("protocol_id") != PROTOCOL_ID
        or preparation.get("status") != "ready_for_orca"
        or preparation.get("target_count") != 25
        or preparation.get("task_count") != 50
        or preparation.get("protonation_subprotocol", {}).get("id")
        != "pdbfixer_standard_only_rng20260914_openmm_cpu_threads1_v1"
    ):
        raise ScoreError("preparation is not the complete frozen panel")
    expected_pin_record = preparation.get("implementation_pins")
    if expected_pin_record != {"path": str(pins_path), "sha256": sha256_file(pins_path)}:
        raise ScoreError("preparation and supplied implementation pins differ")

    core_map = load_core_map(core_map_path)
    targets = preparation.get("targets")
    if not isinstance(targets, list) or len(targets) != 25:
        raise ScoreError("preparation target records are incomplete")
    scores: list[dict[str, Any]] = []
    for target in targets:
        if not isinstance(target, dict):
            raise ScoreError("malformed preparation target record")
        panel_id = target.get("panel_id")
        if not isinstance(panel_id, str) or panel_id not in core_map:
            raise ScoreError(f"unknown prepared panel ID {panel_id!r}")
        row = core_map[panel_id]
        if target.get("class") != row["class"]:
            raise ScoreError(f"prepared class changed for {panel_id}")
        acidic = row["include_plus2_as_ln_specific_asp"].lower() == "true"
        if target.get("acidic_Dplus2") is not acidic:
            raise ScoreError(f"prepared D+2 identity changed for {panel_id}")
        manifest_path = verify_file_record(target.get("manifest"), f"{panel_id} carve manifest")
        manifest_hash = sha256_file(manifest_path)
        manifest = read_object(manifest_path)
        if manifest.get("protocol_id") != PROTOCOL_ID or manifest.get("panel_id") != panel_id:
            raise ScoreError(f"wrong carve protocol/panel identity for {panel_id}")
        carver_pin = pins["experiment_helpers"]["fixed_core_carver"]
        if (
            manifest.get("carver", {}).get("path") != carver_pin["path"]
            or manifest.get("carver", {}).get("implementation_sha256")
            != carver_pin["sha256"]
            or manifest.get("experiment_provenance", {}).get("implementation_pins")
            != {"path": str(pins_path), "sha256": sha256_file(pins_path)}
            or manifest.get("fixed_core", {}).get("water_policy")
            != "dry_exclude_all_source_and_synthetic_waters"
            or manifest.get("fixed_core", {}).get("point_charge_embedding") is not False
            or manifest.get("protonation_manifest", {}).get(
                "experiment_protonation_protocol_id"
            )
            != "pdbfixer_standard_only_rng20260914_openmm_cpu_threads1_v1"
        ):
            raise ScoreError(f"frozen carve provenance/policy changed for {panel_id}")
        heavy_check_path = verify_file_record(
            target.get("heavy_coordinate_check"), f"{panel_id} heavy-coordinate check"
        )
        if manifest.get("heavy_coordinate_check", {}).get("sha256") != sha256_file(
            heavy_check_path
        ):
            raise ScoreError(f"heavy-coordinate record changed for {panel_id}")
        expected_charges = {"La": -2, "Ca": -3} if acidic else {"La": -1, "Ca": -2}
        charge_ledger = manifest.get("charge_ledger")
        if (
            not isinstance(charge_ledger, dict)
            or charge_ledger.get("expected_total_charges") != expected_charges
            or charge_ledger.get("La_total") != expected_charges["La"]
            or charge_ledger.get("Ca_total") != expected_charges["Ca"]
            or target.get("charges") != expected_charges
        ):
            raise ScoreError(f"frozen charge ledger changed for {panel_id}")
        tasks_raw = manifest.get("tasks")
        if not isinstance(tasks_raw, list):
            raise ScoreError(f"{panel_id} manifest lacks tasks")
        tasks = {task.get("task_id"): task for task in tasks_raw if isinstance(task, dict)}
        if set(tasks) != {"La", "Ca"}:
            raise ScoreError(f"{panel_id} manifest does not have exactly La/Ca")
        artifacts: dict[str, Any] = {}
        energies: dict[str, float] = {}
        for task_id in ("La", "Ca"):
            energy, task_record = score_task(
                task_id=task_id,
                task=tasks[task_id],
                manifest_path=manifest_path,
                manifest_hash=manifest_hash,
                pins=pins,
            )
            energies[task_id] = energy
            artifacts[task_id] = task_record
        raw = (energies["Ca"] - energies["La"]) * HA2KCAL
        aquo_delta = float(pins["aquo_reference"]["delta_E_aquo_hartree"])
        gauge = (energies["Ca"] - energies["La"] - aquo_delta) * HA2KCAL
        scores.append(
            {
                "panel_index": target["panel_index"],
                "panel_id": panel_id,
                "class": row["class"],
                "acidic_Dplus2": acidic,
                "expected_charges": expected_charges,
                "energies_hartree": energies,
                "R_kcal_mol": raw,
                "S_aquo_gauge_kcal_mol": gauge,
                "artifacts": artifacts,
            }
        )
    scores.sort(key=lambda item: item["panel_index"])
    if len({row["panel_id"] for row in scores}) != 25:
        raise ScoreError("scored panel IDs are not unique")

    ca_scores = [row["R_kcal_mol"] for row in scores if row["class"] == "Ca"]
    la_scores = [row["R_kcal_mol"] for row in scores if row["class"] == "La"]
    upper_ca = max(ca_scores)
    lower_la = min(la_scores)
    strict = lower_la > upper_ca
    gap = lower_la - upper_ca
    candidate_threshold = (upper_ca + lower_la) / 2.0 if strict else None
    loo = leave_one_out(scores)
    loo_correct = sum(record["correct"] is True for record in loo)
    auc = auc_larger_is_la(scores)
    gates = {
        "all_25_valid": len(scores) == 25,
        "strict_separation": strict,
        "gap_at_least_5_kcal_mol": strict and gap >= MIN_GAP_KCAL,
        "strict_LOO_midpoint_25_of_25": loo_correct == 25,
    }
    calibratable = all(gates.values())
    threshold = candidate_threshold if calibratable else None
    aquo_A = float(pins["aquo_reference"]["delta_E_aquo_hartree"]) * HA2KCAL
    upper_ca_S = upper_ca - aquo_A
    lower_la_S = lower_la - aquo_A
    threshold_S = threshold - aquo_A if threshold is not None else None

    motif_correct = sum(
        (("La" if row["acidic_Dplus2"] else "Ca") == row["class"])
        for row in scores
    )
    if motif_correct != 25:
        raise ScoreError("frozen motif-only baseline no longer classifies 25/25")
    statement = (
        "The D+2 Asp motif-only baseline classifies all 25/25 controls, and that motif "
        "also determines the fixed-core atom count and net charge. Therefore this panel "
        "cannot support any claim that DFT adds discrimination beyond the motif/charge "
        "confound, even if every preregistered calibration gate passes."
    )
    result = {
        "schema_version": RESULT_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "scored_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "preregistration": {"path": str(preregistration), "sha256": sha256_file(preregistration)},
        "preparation": {"path": str(preparation_path), "sha256": sha256_file(preparation_path)},
        "implementation_pins": {"path": str(pins_path), "sha256": sha256_file(pins_path)},
        "primary_contrast": "R=(E_Ca_site-E_La_site)*627.509474_kcal_per_mol",
        "aquo_reporting_gauge": {
            **pins["aquo_reference"],
            "A_kcal_mol": aquo_A,
            "classification_use": False,
            "registry_compatibility_claimed": False,
        },
        "scores": scores,
        "calibration": {
            "verdict": "CALIBRATABLE" if calibratable else "NOT CALIBRATABLE",
            "valid_pair_count": len(scores),
            "U_max_Ca_kcal_mol": upper_ca,
            "L_min_La_kcal_mol": lower_la,
            "strict_separation": strict,
            "gap_kcal_mol": gap,
            "minimum_preregistered_gap_kcal_mol": MIN_GAP_KCAL,
            "gap_floor_pass": gates["gap_at_least_5_kcal_mol"],
            "threshold_R_kcal_mol": threshold,
            "U_max_Ca_S_kcal_mol": upper_ca_S,
            "L_min_La_S_kcal_mol": lower_la_S,
            "threshold_S_kcal_mol": threshold_S,
            "R_and_S_gates_identical": True,
            "threshold_released_only_if_all_gates_pass": True,
            "released_supported_bands": (
                {
                    "Ca": {"R_kcal_mol_max": upper_ca, "S_kcal_mol_max": upper_ca_S},
                    "indeterminate": {
                        "R_kcal_mol_open_interval": [upper_ca, lower_la],
                        "S_kcal_mol_open_interval": [upper_ca_S, lower_la_S],
                    },
                    "La": {"R_kcal_mol_min": lower_la, "S_kcal_mol_min": lower_la_S},
                }
                if calibratable
                else None
            ),
            "auroc": auc,
            "loo_correct_count": loo_correct,
            "loo_records": loo,
            "gates": gates,
            "posthoc_threshold_search_performed": False,
        },
        "interpretation": {
            "motif_only_rule": "predict La iff catalytic-Asp+2 is acidic Asp",
            "motif_only_correct_count": motif_correct,
            "motif_only_total": 25,
            "motif_core_charge_perfectly_confounded_with_class": True,
            "incremental_discrimination_claim_supported": False,
            "incremental_discrimination_statement": statement,
            "reserved_crystal_pair_consumed": False,
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "result.json"
    markdown_path = output_dir / "RESULT.md"
    if json_path.exists() or markdown_path.exists():
        raise ScoreError("refusing to overwrite an existing scientific result")
    write_new(json_path, json.dumps(result, indent=2, sort_keys=True) + "\n")
    write_new(markdown_path, render_markdown(result))
    return json_path, markdown_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--preparation", type=Path, default=HERE / "prepared" / "preparation.json"
    )
    parser.add_argument(
        "--implementation-pins", type=Path, default=HERE / "implementation_pins.json"
    )
    parser.add_argument("--output-dir", type=Path, default=HERE)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        json_path, markdown_path = score(
            args.preparation, args.implementation_pins, args.output_dir
        )
    except (ScoreError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json_path)
    print(markdown_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
