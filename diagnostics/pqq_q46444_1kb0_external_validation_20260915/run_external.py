#!/usr/bin/env python3
"""Run the audited Q46444/1KB0 La/Ca pair inside one SLURM allocation."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Mapping, Sequence

import prepare_external as preparation_module


HERE = Path(__file__).resolve().parent
PROTOCOL_ID = preparation_module.PROTOCOL_ID
EXECUTION_SCHEMA = "alchemical_bvs.pqq_q46444_1kb0_external_execution.v1"


class ExternalRunError(RuntimeError):
    """The external preparation or SLURM execution environment is invalid."""


def read_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ExternalRunError(f"JSON root is not an object: {path}")
    return value


def write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def allocation() -> dict[str, Any]:
    job_id = os.environ.get("SLURM_JOB_ID")
    raw_cpus = os.environ.get("SLURM_CPUS_ON_NODE", "")
    raw_nodes = os.environ.get("SLURM_JOB_NUM_NODES", "")
    match = re.fullmatch(r"([0-9]+)", raw_cpus)
    partition = os.environ.get("SLURM_JOB_PARTITION", "")
    if not job_id:
        raise ExternalRunError("refusing ORCA execution outside SLURM")
    if raw_nodes != "1":
        raise ExternalRunError(f"external runner requires one node, not {raw_nodes!r}")
    if match is None or int(match.group(1)) < 2:
        raise ExternalRunError(f"invalid SLURM_CPUS_ON_NODE={raw_cpus!r}")
    if partition in {"gpu", "gpu_h200", "test"}:
        raise ExternalRunError(f"refusing CPU-only ORCA work on partition {partition!r}")
    return {
        "slurm_job_id": job_id,
        "slurm_job_partition": partition,
        "slurm_job_num_nodes": 1,
        "slurm_cpus_on_node": int(match.group(1)),
        "host": os.uname().nodename,
    }


def verify_preparation(
    preparation_path: Path, pins_path: Path, pins: Mapping[str, Any], gate: Mapping[str, Any]
) -> Path:
    prep = read_object(preparation_path)
    if (
        prep.get("schema_version") != preparation_module.PREPARATION_SCHEMA
        or prep.get("protocol_id") != PROTOCOL_ID
        or prep.get("status") != "ready_for_orca_after_passing_calibration_gate"
        or prep.get("experiment_id") != "Q46444_1KB0_external_Ca_structural_transfer"
        or prep.get("expected_band") != "Ca-supported"
        or prep.get("calibration_gate") != gate
        or prep.get("external_pins")
        != preparation_module.holdout.file_record(pins_path)
        or prep.get("target_count") != 1
        or prep.get("task_count") != 2
        or prep.get("all_nonmetal_arm_coordinates_byte_identical") is not True
        or prep.get("orca_executed") is not False
    ):
        raise ExternalRunError("preparation is not the preregistered external release")
    target = prep.get("target")
    if not isinstance(target, dict) or target.get("pdb_id") != "1KB0":
        raise ExternalRunError("prepared external target is not 1KB0")
    manifest_path = preparation_module.verify_record(target.get("manifest"), "1KB0 manifest")
    manifest = read_object(manifest_path)
    if (
        manifest.get("protocol_id") != PROTOCOL_ID
        or manifest.get("panel_id") != "1KB0"
        or manifest.get("status") != "ready_for_orca_after_passing_calibration_gate"
        or manifest.get("calibration_gate") != gate
        or manifest.get("coordination", {}).get("coordination_number") != 7
        or manifest.get("charge_ledger", {}).get("expected_total_charges")
        != {"La": -1, "Ca": -2}
        or manifest.get("paired_arm_invariant", {}).get(
            "nonmetal_coordinates_byte_identical"
        )
        is not True
        or [task.get("task_id") for task in manifest.get("tasks", [])] != ["La", "Ca"]
    ):
        raise ExternalRunError("1KB0 manifest violates the external contract")
    policy = manifest.get("execution_policy", {})
    if (
        policy.get("task_runner") != pins["inherited"]["run_orca_task_manifest"]
        or policy.get("runtime_renderer") != pins["inherited"]["render_orca_runtime_input"]
    ):
        raise ExternalRunError("1KB0 manifest execution helpers differ from pins")
    return manifest_path


def run_external(preparation_path: Path, pins_path: Path) -> Path:
    allocation_record = allocation()
    pins_path = pins_path.resolve()
    pins, calibration_pins, _, gate = preparation_module.verify_external_pins(pins_path)
    own_record = pins["implementation"]["run_external"]
    if Path(own_record["path"]).resolve() != Path(__file__).resolve():
        raise ExternalRunError("executing runner differs from external pins")
    manifest = verify_preparation(preparation_path.resolve(), pins_path, pins, gate)
    runner_record = pins["inherited"]["run_orca_task_manifest"]
    runner = preparation_module.verify_record(runner_record, "manifested ORCA runner")
    orca_record = calibration_pins["orca_runtime"]["executable"]
    orca = preparation_module.verify_record(orca_record, "ORCA executable")
    ranks_per_leg = min(16, allocation_record["slurm_cpus_on_node"] // 2)
    if ranks_per_leg < 1:
        raise ExternalRunError("allocation cannot support the paired ORCA arms")

    log_path = manifest.parent / "external_orca_runner.log"
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
        runner_record["sha256"],
    ]
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    with log_path.open("ab") as handle:
        handle.write(f"\n[{started}] {' '.join(command)}\n".encode())
        completed = subprocess.run(
            command, stdout=handle, stderr=subprocess.STDOUT, check=False
        )

    receipt_path = preparation_path.resolve().parent / (
        f"external_execution_{allocation_record['slurm_job_id']}.json"
    )
    if receipt_path.exists():
        raise ExternalRunError(f"refusing to overwrite receipt {receipt_path}")
    status = "complete" if completed.returncode == 0 else "failed"
    receipt = {
        "schema_version": EXECUTION_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "started_at_utc": started,
        "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": status,
        "preparation": preparation_module.holdout.file_record(preparation_path.resolve()),
        "external_pins": preparation_module.holdout.file_record(pins_path),
        "runner": own_record,
        "manifested_orca_runner": runner_record,
        "orca_executable": orca_record,
        "allocation": allocation_record,
        "target_id": "1KB0",
        "task_count": 2,
        "parallelism": {
            "simultaneous_ORCA_legs": 2,
            "mpi_ranks_per_leg": ranks_per_leg,
            "assigned_mpi_ranks": 2 * ranks_per_leg,
            "unassigned_CPUs": allocation_record["slurm_cpus_on_node"] - 2 * ranks_per_leg,
            "omp_threads_per_rank": 1,
        },
        "returncode": completed.returncode,
        "manifest": preparation_module.holdout.file_record(manifest),
        "log": preparation_module.holdout.file_record(log_path),
    }
    write_json_atomic(receipt_path, receipt)
    if completed.returncode != 0:
        raise ExternalRunError(f"1KB0 ORCA pair failed with status {completed.returncode}")
    return receipt_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation", type=Path, required=True)
    parser.add_argument("--pins", type=Path, default=HERE / "implementation_pins.json")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = run_external(args.preparation, args.pins)
    except (
        ExternalRunError,
        preparation_module.ExternalPreparationError,
        preparation_module.holdout.HoldoutPreparationError,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

