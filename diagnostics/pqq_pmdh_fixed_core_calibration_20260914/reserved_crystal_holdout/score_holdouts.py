#!/usr/bin/env python3
"""Score the fixed-core crystal holdouts against frozen calibration bands."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import prepare_holdouts as holdout
import run_holdouts as runner

sys.path.insert(0, str(holdout.CALIBRATION_DIR))
import score as calibration_score  # noqa: E402


HERE = Path(__file__).resolve().parent
PROTOCOL_ID = holdout.PROTOCOL_ID
RESULT_SCHEMA = "alchemical_bvs.pqq_fixed_core_holdout_result.v1"
HA2KCAL = 627.509474
PRIMARY_IDS = runner.PRIMARY_IDS
FULL_IDS = runner.FULL_IDS


class HoldoutScoreError(RuntimeError):
    """The completed holdout set cannot be scored under the frozen rule."""


def read_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise HoldoutScoreError(f"cannot read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise HoldoutScoreError(f"JSON root is not an object: {path}")
    return value


def write_new(path: Path, text: str) -> None:
    with path.open("x") as handle:
        handle.write(text)


def verify_execution_receipt(
    receipt_path: Path,
    preparation_path: Path,
    pins_path: Path,
    pins: Mapping[str, Any],
    calibration_pins: Mapping[str, Any],
    expected_ids: Sequence[str],
) -> dict[str, Any]:
    receipt_path = receipt_path.resolve()
    receipt = read_object(receipt_path)
    expected_ids = tuple(expected_ids)
    if expected_ids not in (PRIMARY_IDS, FULL_IDS):
        raise HoldoutScoreError("scorer received an invalid expected target set")
    target_count = len(expected_ids)
    leg_count = 2 * target_count
    targets = receipt.get("targets")
    if (
        not isinstance(targets, list)
        or any(not isinstance(item, dict) for item in targets)
        or [item.get("pdb_id") for item in targets] != list(expected_ids)
        or any(item.get("status") not in {"complete", "failed"} for item in targets)
    ):
        raise HoldoutScoreError("execution receipt has an invalid target ledger")
    failed_ids = sorted(
        item["pdb_id"] for item in targets if item.get("status") == "failed"
    )
    expected_status = "failed" if failed_ids else "complete"
    expected_secondary_included = expected_ids == FULL_IDS
    parallelism = receipt.get("parallelism")
    allocation = receipt.get("allocation")
    if not isinstance(parallelism, dict) or not isinstance(allocation, dict):
        raise HoldoutScoreError("execution receipt lacks allocation/parallelism ledgers")
    total_cpus = allocation.get("slurm_cpus_on_node")
    ranks_per_leg = parallelism.get("mpi_ranks_per_leg")
    if (
        isinstance(total_cpus, bool)
        or not isinstance(total_cpus, int)
        or isinstance(ranks_per_leg, bool)
        or not isinstance(ranks_per_leg, int)
    ):
        raise HoldoutScoreError("execution receipt has malformed CPU/rank counts")
    expected_ranks_per_leg = min(16, total_cpus // leg_count)
    if (
        receipt.get("schema_version") != runner.EXECUTION_SCHEMA
        or receipt.get("protocol_id") != PROTOCOL_ID
        or receipt.get("status") != expected_status
        or receipt.get("preparation") != holdout.file_record(preparation_path)
        or receipt.get("holdout_pins") != holdout.file_record(pins_path)
        or receipt.get("failed_targets") != failed_ids
        or receipt.get("target_ids") != list(expected_ids)
        or receipt.get("target_count") != target_count
        or receipt.get("task_count") != leg_count
        or receipt.get("secondary_included") is not expected_secondary_included
        or parallelism.get("simultaneous_targets") != target_count
        or parallelism.get("simultaneous_legs_per_target") != 2
        or parallelism.get("simultaneous_ORCA_legs") != leg_count
        or ranks_per_leg != expected_ranks_per_leg
        or not 1 <= ranks_per_leg <= 16
        or parallelism.get("assigned_mpi_ranks") != leg_count * ranks_per_leg
        or parallelism.get("unassigned_CPUs") != total_cpus - leg_count * ranks_per_leg
        or parallelism.get("omp_threads_per_rank") != 1
        or parallelism.get("aggregate_ORCA_maxcore_upper_bound_MB")
        != leg_count * ranks_per_leg * 8000
        or receipt.get("runner") != pins["holdout_implementation"]["run_holdouts"]
        or receipt.get("manifested_orca_runner")
        != pins["inherited_helpers"]["run_orca_task_manifest"]
        or receipt.get("orca_executable")
        != calibration_pins["orca_runtime"]["executable"]
    ):
        raise HoldoutScoreError(
            f"execution receipt is not a valid attempted {leg_count}-leg holdout run"
        )
    return receipt


def frozen_s_bands(
    preparation: Mapping[str, Any],
    calibration_pins: Mapping[str, Any],
    calibration_pins_path: Path,
) -> dict[str, float]:
    gate = preparation.get("calibration_gate")
    if not isinstance(gate, dict):
        raise HoldoutScoreError("holdout preparation lacks calibration gate")
    calibration_result_path = holdout.verify_file_record(
        gate.get("result"), "locked calibration result"
    )
    verified_gate = holdout.verify_calibration_gate(
        calibration_result_path, calibration_pins_path, calibration_pins
    )
    for field in (
        "U_max_Ca_kcal_mol",
        "L_min_La_kcal_mol",
        "gap_kcal_mol",
        "threshold_R_kcal_mol",
        "U_max_Ca_S_kcal_mol",
        "L_min_La_S_kcal_mol",
        "threshold_S_kcal_mol",
        "aquo_A_kcal_mol",
    ):
        if gate.get(field) != verified_gate.get(field):
            raise HoldoutScoreError(f"holdout preparation changed frozen gate field {field}")
    calibration_result = read_object(calibration_result_path)
    calibration = calibration_result.get("calibration")
    if not isinstance(calibration, dict):
        raise HoldoutScoreError("locked calibration result lacks released bands")
    try:
        upper_s = float(calibration["U_max_Ca_S_kcal_mol"])
        lower_s = float(calibration["L_min_La_S_kcal_mol"])
        threshold_s = float(calibration["threshold_S_kcal_mol"])
        upper_r = float(calibration["U_max_Ca_kcal_mol"])
        lower_r = float(calibration["L_min_La_kcal_mol"])
        threshold_r = float(calibration["threshold_R_kcal_mol"])
        aquo_delta_hartree = float(calibration_pins["aquo_reference"]["delta_E_aquo_hartree"])
    except (KeyError, TypeError, ValueError) as exc:
        raise HoldoutScoreError("released S bands or aquo gauge are malformed") from exc
    values = (
        upper_s,
        lower_s,
        threshold_s,
        upper_r,
        lower_r,
        threshold_r,
        aquo_delta_hartree,
    )
    if not all(math.isfinite(value) for value in values):
        raise HoldoutScoreError("released S bands or aquo gauge are nonfinite")
    aquo_a = aquo_delta_hartree * HA2KCAL
    if (
        abs(upper_s - (upper_r - aquo_a)) > 1.0e-9
        or abs(lower_s - (lower_r - aquo_a)) > 1.0e-9
        or abs(threshold_s - (threshold_r - aquo_a)) > 1.0e-9
        or not upper_s < lower_s
    ):
        raise HoldoutScoreError("released S bands are inconsistent with frozen R/aquo values")
    bands = calibration.get("released_supported_bands")
    if (
        not isinstance(bands, dict)
        or bands.get("Ca", {}).get("S_kcal_mol_max") != upper_s
        or bands.get("La", {}).get("S_kcal_mol_min") != lower_s
        or bands.get("indeterminate", {}).get("S_kcal_mol_open_interval")
        != [upper_s, lower_s]
    ):
        raise HoldoutScoreError("released S-band record changed")
    return {
        "U_max_Ca_S_kcal_mol": upper_s,
        "L_min_La_S_kcal_mol": lower_s,
        "threshold_S_kcal_mol": threshold_s,
        "aquo_A_kcal_mol": aquo_a,
        "calibration_result_path": str(calibration_result_path),
        "calibration_result_sha256": holdout.sha256_file(calibration_result_path),
    }


def classify_s(value: float, *, upper: float, lower: float) -> str:
    if value <= upper:
        return "Ca-supported"
    if value >= lower:
        return "Ln-supported"
    return "indeterminate-calibration-gap"


def render_markdown(result: Mapping[str, Any]) -> str:
    bands = result["frozen_S_bands"]
    lines = [
        "# Reserved PQQ-MDH crystal transfer result",
        "",
        f"**Protocol:** `{PROTOCOL_ID}`  ",
        f"**Primary verdict:** **{result['primary_transfer_verdict']}**",
        "",
        "The primary rule was frozen before these energies were read: 1H4I must have "
        f"S <= {bands['U_max_Ca_S_kcal_mol']:.9f} kcal/mol and 4MAE must have "
        f"S >= {bands['L_min_La_S_kcal_mol']:.9f} kcal/mol. The open interval is "
        "indeterminate and fails the primary transfer test.",
        "",
        "| PDB | Role | Expected | S, kcal/mol | Frozen-band call | Primary pass |",
        "|---|---|---|---:|---|---:|",
    ]
    for row in result["scores"]:
        score = row["S_aquo_gauge_kcal_mol"]
        rendered_score = "NA" if score is None else f"{score:.9f}"
        lines.append(
            f"| {row['pdb_id']} | {row['holdout_role']} | {row['expected_band']} | "
            f"{rendered_score} | {row['frozen_band_call']} | "
            f"{row['primary_pass']} |"
        )
        if row["unscorable_reasons"]:
            lines.append(
                f"|  |  | reason |  | {'; '.join(row['unscorable_reasons'])} |  |"
            )
    lines.append("")
    if result["secondary_6OC6"]["status"] == "secondary-not-run":
        lines.append(
            "6OC6 was not run; it remains a nonindependent optional secondary geometry check."
        )
    else:
        lines.append("6OC6 is a secondary geometry check only and did not affect the verdict.")
    lines.extend(
        [
            "4MAE was evaluated after explicit removal of coordinating 15P603/OXT "
            "without replacement; it is a dry fixed-coordinate transfer test with a ligand vacancy.",
            "No threshold or band was fit, shifted, or widened using holdout results.",
            "",
        ]
    )
    return "\n".join(lines)


def score_holdouts(
    *,
    preparation_path: Path,
    execution_receipt_path: Path,
    pins_path: Path,
    output_dir: Path,
) -> tuple[Path, Path, str]:
    preparation_path = preparation_path.resolve()
    pins_path = pins_path.resolve()
    output_dir = output_dir.resolve()
    pins, calibration_pins, _ = holdout.verify_holdout_pins(pins_path)
    own_record = pins["holdout_implementation"]["score_holdouts"]
    if Path(own_record["path"]).resolve() != Path(__file__).resolve():
        raise HoldoutScoreError("executing scorer differs from holdout pins")
    calibration_score_record = pins["inherited_helpers"]["calibration_score"]
    if Path(calibration_score_record["path"]).resolve() != Path(
        calibration_score.__file__
    ).resolve():
        raise HoldoutScoreError("imported execution/energy verifier differs from pins")
    preparation = read_object(preparation_path)
    manifests = runner.verify_preparation(
        preparation_path, pins_path, pins, calibration_pins
    )
    expected_ids = tuple(pdb_id for pdb_id, _ in manifests)
    if expected_ids not in (PRIMARY_IDS, FULL_IDS):
        raise HoldoutScoreError("scoring requires the atomic primary pair, optionally plus 6OC6")
    receipt = verify_execution_receipt(
        execution_receipt_path,
        preparation_path,
        pins_path,
        pins,
        calibration_pins,
        expected_ids,
    )
    calibration_pins_path = Path(
        pins["calibration_implementation_pins"]["path"]
    ).resolve()
    bands = frozen_s_bands(preparation, calibration_pins, calibration_pins_path)
    upper = bands["U_max_Ca_S_kcal_mol"]
    lower = bands["L_min_La_S_kcal_mol"]
    aquo_a = bands["aquo_A_kcal_mol"]
    target_records = {item["pdb_id"]: item for item in preparation["targets"]}
    receipt_targets = {item["pdb_id"]: item for item in receipt["targets"]}
    scores: list[dict[str, Any]] = []
    for pdb_id, manifest_path in manifests:
        manifest = read_object(manifest_path)
        tasks_raw = manifest.get("tasks")
        if not isinstance(tasks_raw, list):
            raise HoldoutScoreError(f"{pdb_id} manifest lacks tasks")
        tasks = {task.get("task_id"): task for task in tasks_raw if isinstance(task, dict)}
        if set(tasks) != {"La", "Ca"}:
            raise HoldoutScoreError(f"{pdb_id} manifest does not have exactly La/Ca")
        manifest_hash = holdout.sha256_file(manifest_path)
        energies: dict[str, float] = {}
        artifacts: dict[str, Any] = {}
        unscorable_reasons: list[str] = []
        receipt_target = receipt_targets[pdb_id]
        if receipt_target["status"] != "complete":
            unscorable_reasons.append(
                f"target runner failed with returncode {receipt_target.get('returncode')!r}"
            )
        for task_id in ("La", "Ca"):
            try:
                energy, record = calibration_score.score_task(
                    task_id=task_id,
                    task=tasks[task_id],
                    manifest_path=manifest_path,
                    manifest_hash=manifest_hash,
                    pins=calibration_pins,
                )
            except (calibration_score.ScoreError, OSError, ValueError) as exc:
                unscorable_reasons.append(f"{task_id}: {type(exc).__name__}: {exc}")
            else:
                energies[task_id] = energy
                artifacts[task_id] = record
        raw_r: float | None = None
        score_s: float | None = None
        call = "unscorable"
        if not unscorable_reasons and set(energies) == {"La", "Ca"}:
            raw_r = (energies["Ca"] - energies["La"]) * HA2KCAL
            score_s = raw_r - aquo_a
            if math.isfinite(raw_r) and math.isfinite(score_s):
                call = classify_s(score_s, upper=upper, lower=lower)
            else:
                unscorable_reasons.append("paired R/S score is nonfinite")
                raw_r = None
                score_s = None
        holdout_role = target_records[pdb_id]["holdout_role"]
        expected_band = (
            "Ca-supported"
            if pdb_id == "1H4I"
            else "Ln-supported"
            if pdb_id == "4MAE"
            else "secondary-only"
        )
        primary_pass: bool | None = (
            call == expected_band if holdout_role == "primary" else None
        )
        scores.append(
            {
                "pdb_id": pdb_id,
                "holdout_role": holdout_role,
                "expected_band": expected_band,
                "energies_hartree": energies,
                "R_kcal_mol": raw_r,
                "S_aquo_gauge_kcal_mol": score_s,
                "frozen_band_call": call,
                "primary_pass": primary_pass,
                "score_status": "unscorable" if unscorable_reasons else "scored",
                "unscorable_reasons": unscorable_reasons,
                "artifacts": artifacts,
            }
        )
    scores.sort(key=lambda item: FULL_IDS.index(item["pdb_id"]))
    primary = [item for item in scores if item["holdout_role"] == "primary"]
    if [item["pdb_id"] for item in primary] != ["1H4I", "4MAE"]:
        raise HoldoutScoreError("primary holdout identity/order changed")
    primary_pass = all(item["primary_pass"] is True for item in primary)
    secondary_score = next(
        (item for item in scores if item["pdb_id"] == runner.SECONDARY_ID),
        None,
    )
    secondary_record = (
        {
            "status": "secondary-not-run",
            "included_in_run": False,
            "included_in_primary_verdict": False,
            "sequence_already_represented_by": "C5B120",
            "frozen_band_call": None,
        }
        if secondary_score is None
        else {
            "status": secondary_score["score_status"],
            "included_in_run": True,
            "included_in_primary_verdict": False,
            "sequence_already_represented_by": "C5B120",
            "frozen_band_call": secondary_score["frozen_band_call"],
        }
    )
    result = {
        "schema_version": RESULT_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "scored_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "primary_transfer_verdict": "PASS" if primary_pass else "FAIL",
        "primary_pass_rule": {
            "1H4I": "S <= U_max_Ca_S_kcal_mol",
            "4MAE": "S >= L_min_La_S_kcal_mol",
            "open_gap": "indeterminate_and_primary_failure",
            "wrong_band": "primary_failure",
            "unscorable": "primary_failure",
        },
        "frozen_S_bands": bands,
        "threshold_refit_performed": False,
        "holdout_results_used_to_modify_bands": False,
        "preparation": holdout.file_record(preparation_path),
        "execution_receipt": holdout.file_record(execution_receipt_path),
        "holdout_pins": holdout.file_record(pins_path),
        "calibration_implementation_pins": holdout.file_record(calibration_pins_path),
        "executed_target_ids": list(expected_ids),
        "scores": scores,
        "secondary_6OC6": secondary_record,
        "interpretation_limits": {
            "4MAE_dry_15P_vacancy": True,
            "1H4I_mature_sequence_overlaps_P16027_calibration": True,
            "prior_protocol_results_exist_for_1H4I_and_4MAE": True,
            "blind_sequence_independent_validation": False,
        },
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    json_path = output_dir / "holdout_result.json"
    markdown_path = output_dir / "HOLDOUT_RESULT.md"
    if json_path.exists() or markdown_path.exists():
        raise HoldoutScoreError("refusing to overwrite an existing holdout result")
    write_new(json_path, json.dumps(result, indent=2, sort_keys=True) + "\n")
    write_new(markdown_path, render_markdown(result))
    return json_path, markdown_path, result["primary_transfer_verdict"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation", type=Path, required=True)
    parser.add_argument("--execution-receipt", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--holdout-pins",
        type=Path,
        default=HERE / "holdout_implementation_pins.json",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        json_path, markdown_path, verdict = score_holdouts(
            preparation_path=args.preparation,
            execution_receipt_path=args.execution_receipt,
            pins_path=args.holdout_pins,
            output_dir=args.output_dir,
        )
    except (
        HoldoutScoreError,
        holdout.HoldoutPreparationError,
        runner.HoldoutRunError,
        calibration_score.ScoreError,
        OSError,
        ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    for path in (json_path, markdown_path):
        print(path)
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
