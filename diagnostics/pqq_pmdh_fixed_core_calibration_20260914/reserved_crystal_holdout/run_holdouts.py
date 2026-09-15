#!/usr/bin/env python3
"""Run an atomic PQQ-MDH primary pair, with optional secondary, under SLURM."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Mapping, Sequence

import prepare_holdouts as holdout


HERE = Path(__file__).resolve().parent
PROTOCOL_ID = holdout.PROTOCOL_ID
EXECUTION_SCHEMA = "alchemical_bvs.pqq_fixed_core_holdout_execution.v1"
PRIMARY_IDS = ("1H4I", "4MAE")
SECONDARY_ID = "6OC6"
FULL_IDS = (*PRIMARY_IDS, SECONDARY_ID)


class HoldoutRunError(RuntimeError):
    """The prepared holdout set or SLURM allocation is invalid."""


def read_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise HoldoutRunError(f"cannot read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise HoldoutRunError(f"JSON root is not an object: {path}")
    return value


def write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def allocation() -> dict[str, Any]:
    job_id = os.environ.get("SLURM_JOB_ID")
    raw_cpus = os.environ.get("SLURM_CPUS_ON_NODE", "")
    cpu_match = re.fullmatch(r"([0-9]+)", raw_cpus)
    raw_nodes = os.environ.get("SLURM_JOB_NUM_NODES", "")
    if not job_id:
        raise HoldoutRunError("refusing ORCA execution outside SLURM")
    if cpu_match is None or int(cpu_match.group(1)) < 4:
        raise HoldoutRunError(f"invalid SLURM_CPUS_ON_NODE={raw_cpus!r}")
    if raw_nodes != "1":
        raise HoldoutRunError(f"holdout runner requires one node, not {raw_nodes!r}")
    if os.environ.get("SLURM_JOB_PARTITION") != "high-memory":
        raise HoldoutRunError("holdout release runner requires the high-memory partition")
    return {
        "slurm_job_id": job_id,
        "slurm_job_partition": os.environ["SLURM_JOB_PARTITION"],
        "slurm_job_num_nodes": 1,
        "slurm_cpus_on_node": int(cpu_match.group(1)),
        "host": os.uname().nodename,
    }


def verify_preparation(
    preparation_path: Path,
    pins_path: Path,
    pins: Mapping[str, Any],
    calibration_pins: Mapping[str, Any],
) -> list[tuple[str, Path]]:
    preparation_path = preparation_path.resolve()
    preparation = read_object(preparation_path)
    secondary_included = preparation.get("secondary_included")
    if secondary_included is True:
        expected_ids = FULL_IDS
        expected_scope = "atomic_primary_pair_plus_nonindependent_secondary"
    elif secondary_included is False:
        expected_ids = PRIMARY_IDS
        expected_scope = "atomic_primary_pair"
    else:
        raise HoldoutRunError("preparation has an invalid secondary-inclusion flag")
    expected_target_count = len(expected_ids)
    expected_task_count = 2 * expected_target_count
    if (
        preparation.get("schema_version") != holdout.HOLDOUT_PREPARATION_SCHEMA
        or preparation.get("protocol_id") != PROTOCOL_ID
        or preparation.get("status") != "ready_for_orca_after_passing_calibration_gate"
        or preparation.get("release_scope") != expected_scope
        or preparation.get("release_runnable") is not True
        or preparation.get("primary_holdouts") != ["1H4I", "4MAE"]
        or preparation.get("secondary_holdout") != "6OC6"
        or "secondary_omission_disposition" not in preparation
        or preparation.get("secondary_omission_disposition")
        != (None if secondary_included else "secondary-not-run")
        or preparation.get("target_count") != expected_target_count
        or preparation.get("task_count") != expected_task_count
        or preparation.get("all_nonmetal_arm_coordinates_byte_identical") is not True
        or preparation.get("water_policy") != "dry_exclude_all_source_and_synthetic_waters"
        or preparation.get("orca_executed") is not False
    ):
        raise HoldoutRunError("preparation is not a reviewed atomic holdout release")
    if preparation.get("holdout_pins") != holdout.file_record(pins_path):
        raise HoldoutRunError("preparation is not bound to supplied holdout pins")
    calibration_pins_path = Path(
        pins["calibration_implementation_pins"]["path"]
    ).resolve()
    calibration_result = holdout.verify_file_record(
        pins.get("calibration_result"), "locked calibration result"
    )
    verified_gate = holdout.verify_calibration_gate(
        calibration_result, calibration_pins_path, calibration_pins
    )
    if preparation.get("calibration_gate") != verified_gate:
        raise HoldoutRunError("preparation calibration gate differs from locked result")
    targets = preparation.get("targets")
    if (
        not isinstance(targets, list)
        or any(not isinstance(item, dict) for item in targets)
        or [item.get("pdb_id") for item in targets] != list(expected_ids)
    ):
        raise HoldoutRunError(
            f"prepared target identity/order differs from {','.join(expected_ids)}"
        )
    manifests: list[tuple[str, Path]] = []
    for target in targets:
        if not isinstance(target, dict):
            raise HoldoutRunError("malformed holdout target")
        pdb_id = target["pdb_id"]
        manifest_path = holdout.verify_file_record(target.get("manifest"), f"{pdb_id} manifest")
        manifest = read_object(manifest_path)
        expected_role = "secondary" if pdb_id == "6OC6" else "primary"
        expected_class = "Ca/MxaF" if pdb_id == "1H4I" else "Ln/XoxF"
        expected_cn = 6 if pdb_id == "1H4I" else 9
        expected_charges = (
            {"La": -1, "Ca": -2}
            if pdb_id == "1H4I"
            else {"La": -2, "Ca": -3}
        )
        tasks = manifest.get("tasks")
        fixed_core = manifest.get("fixed_core")
        charge_ledger = manifest.get("charge_ledger")
        pqq = manifest.get("pqq")
        noncore = fixed_core.get("excluded_noncore_direct_ligands") if isinstance(fixed_core, dict) else None
        if (
            manifest.get("schema_version") != holdout.HOLDOUT_CARVE_SCHEMA
            or manifest.get("protocol_id") != PROTOCOL_ID
            or manifest.get("panel_id") != pdb_id
            or manifest.get("status") != "ready_for_orca_after_passing_calibration_gate"
            or manifest.get("holdout", {}).get("role") != expected_role
            or manifest.get("holdout", {}).get("biological_class") != expected_class
            or manifest.get("calibration_gate") != verified_gate
            or manifest.get("experiment_provenance", {}).get("holdout_pins")
            != holdout.file_record(pins_path)
            or manifest.get("fixed_core", {}).get("water_policy")
            != "dry_exclude_all_source_and_synthetic_waters"
            or not isinstance(fixed_core, dict)
            or fixed_core.get("synthetic_water_count") != 0
            or fixed_core.get("replacement_ligand_count") != 0
            or fixed_core.get("point_charge_embedding") is not False
            or fixed_core.get("geometry_relaxation") is not False
            or manifest.get("coordination", {}).get("coordination_number") != expected_cn
            or not isinstance(charge_ledger, dict)
            or charge_ledger.get("expected_total_charges") != expected_charges
            or charge_ledger.get("La_total") != expected_charges["La"]
            or charge_ledger.get("Ca_total") != expected_charges["Ca"]
            or not isinstance(pqq, dict)
            or pqq.get("schema_id") != "pdb_ccd_pqq_v1"
            or pqq.get("microstate_id") != "pqq_ox_3minus_v1"
            or pqq.get("formal_charge") != -3
            or manifest.get("paired_arm_invariant", {}).get(
                "nonmetal_coordinates_byte_identical"
            )
            is not True
            or not isinstance(tasks, list)
            or [task.get("task_id") for task in tasks] != ["La", "Ca"]
        ):
            raise HoldoutRunError(f"{pdb_id} manifest violates the holdout contract")
        if pdb_id == "4MAE":
            if (
                not isinstance(noncore, list)
                or len(noncore) != 1
                or noncore[0].get("selector") != "A:15P603/OXT"
                or noncore[0].get("disposition") != "excluded_without_replacement"
            ):
                raise HoldoutRunError("4MAE manifest lacks explicit dry 15P vacancy")
        elif noncore != []:
            raise HoldoutRunError(f"{pdb_id} has an unexpected noncore direct ligand")
        policy = manifest.get("execution_policy")
        if (
            not isinstance(policy, dict)
            or policy.get("id") != holdout.EXECUTION_POLICY_ID
            or policy.get("task_runner") != pins["inherited_helpers"]["run_orca_task_manifest"]
            or policy.get("runtime_renderer")
            != pins["inherited_helpers"]["render_orca_runtime_input"]
        ):
            raise HoldoutRunError(f"{pdb_id} execution helper pins changed")
        manifests.append((pdb_id, manifest_path))
    return manifests


def run_target(
    *,
    pdb_id: str,
    manifest: Path,
    runner: Path,
    runner_hash: str,
    orca: Path,
    ranks_per_leg: int,
) -> dict[str, Any]:
    log_path = manifest.parent / "holdout_orca_runner.log"
    command = [
        str(Path(sys.executable).resolve()),
        str(runner),
        str(manifest),
        "--orca",
        str(orca),
        "--workers",
        "2",
        "--nprocs",
        str(ranks_per_leg),
        "--expected-runner-sha256",
        runner_hash,
    ]
    with log_path.open("ab") as handle:
        handle.write(
            f"\n[{dt.datetime.now(dt.timezone.utc).isoformat()}] {' '.join(command)}\n".encode()
        )
        completed = subprocess.run(
            command,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return {
        "pdb_id": pdb_id,
        "status": "complete" if completed.returncode == 0 else "failed",
        "returncode": completed.returncode,
        "manifest": holdout.file_record(manifest),
        "log": holdout.file_record(log_path),
    }


def run_holdouts(preparation_path: Path, pins_path: Path) -> Path:
    allocation_record = allocation()
    pins_path = pins_path.resolve()
    pins, calibration_pins, _ = holdout.verify_holdout_pins(pins_path)
    own_record = pins["holdout_implementation"]["run_holdouts"]
    if Path(own_record["path"]).resolve() != Path(__file__).resolve():
        raise HoldoutRunError("executing runner differs from holdout pins")
    manifests = verify_preparation(
        preparation_path, pins_path, pins, calibration_pins
    )
    runner_record = pins["inherited_helpers"]["run_orca_task_manifest"]
    runner = holdout.verify_file_record(runner_record, "manifested ORCA runner")
    orca_record = calibration_pins["orca_runtime"]["executable"]
    orca = holdout.verify_file_record(orca_record, "ORCA executable")
    total_cpus = int(allocation_record["slurm_cpus_on_node"])
    target_ids = tuple(pdb_id for pdb_id, _ in manifests)
    if target_ids not in (PRIMARY_IDS, FULL_IDS):
        raise HoldoutRunError("verified preparation produced an invalid target set")
    target_count = len(target_ids)
    leg_count = 2 * target_count
    # Match the calibration runner's empirically chosen ceiling.  Arbitrary
    # larger %pal values are legal, but these modest fixed-core single points
    # are not expected to scale usefully to 37--57 MPI ranks per leg.
    ranks_per_leg = min(16, total_cpus // leg_count)
    if ranks_per_leg < 1:
        raise HoldoutRunError(
            f"allocation cannot support {leg_count} simultaneous ORCA legs"
        )

    started = dt.datetime.now(dt.timezone.utc).isoformat()
    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=target_count) as executor:
        futures = {
            executor.submit(
                run_target,
                pdb_id=pdb_id,
                manifest=manifest,
                runner=runner,
                runner_hash=runner_record["sha256"],
                orca=orca,
                ranks_per_leg=ranks_per_leg,
            ): pdb_id
            for pdb_id, manifest in manifests
        }
        for future in as_completed(futures):
            pdb_id = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "pdb_id": pdb_id,
                    "status": "failed",
                    "returncode": None,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            results.append(result)
            print(f"{pdb_id}: {result['status']}", flush=True)

    failures = sorted(item["pdb_id"] for item in results if item["status"] != "complete")
    receipt_path = preparation_path.resolve().parent / (
        f"holdout_execution_{allocation_record['slurm_job_id']}.json"
    )
    if receipt_path.exists():
        raise HoldoutRunError(f"refusing to overwrite execution receipt {receipt_path}")
    receipt = {
        "schema_version": EXECUTION_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "started_at_utc": started,
        "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "complete" if not failures else "failed",
        "preparation": holdout.file_record(preparation_path),
        "holdout_pins": holdout.file_record(pins_path),
        "runner": own_record,
        "manifested_orca_runner": runner_record,
        "orca_executable": orca_record,
        "allocation": allocation_record,
        "target_ids": list(target_ids),
        "target_count": target_count,
        "task_count": leg_count,
        "secondary_included": SECONDARY_ID in target_ids,
        "parallelism": {
            "simultaneous_targets": target_count,
            "simultaneous_legs_per_target": 2,
            "simultaneous_ORCA_legs": leg_count,
            "mpi_ranks_per_leg": ranks_per_leg,
            "assigned_mpi_ranks": leg_count * ranks_per_leg,
            "unassigned_CPUs": total_cpus - leg_count * ranks_per_leg,
            "omp_threads_per_rank": 1,
            "rank_cap_rationale": (
                "match calibration cap16; avoid communication-dominated PAL37/PAL57 "
                "for modest fixed-core r2SCAN-3c single points"
            ),
            "aggregate_ORCA_maxcore_upper_bound_MB": leg_count * ranks_per_leg * 8000,
        },
        "failed_targets": failures,
        "targets": sorted(results, key=lambda item: target_ids.index(item["pdb_id"])),
    }
    write_json_atomic(receipt_path, receipt)
    if failures:
        raise HoldoutRunError(f"holdout ORCA failures: {failures}")
    return receipt_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation", type=Path, required=True)
    parser.add_argument(
        "--holdout-pins",
        type=Path,
        default=HERE / "holdout_implementation_pins.json",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        receipt = run_holdouts(args.preparation, args.holdout_pins)
    except (HoldoutRunError, holdout.HoldoutPreparationError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
