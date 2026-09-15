#!/usr/bin/env python3
"""Run every prepared fixed-core La/Ca pair inside one SLURM allocation."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent
PROTOCOL_ID = "pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3"
EXECUTION_SCHEMA = "alchemical_bvs.pqq_fixed_core_panel_execution.v1"


class PanelRunError(RuntimeError):
    """The prepared panel or allocation violates the execution contract."""


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
        raise PanelRunError(f"cannot read JSON object {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PanelRunError(f"JSON root is not an object: {path}")
    return value


def write_json_atomic(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def verify_file_record(record: Any, label: str) -> Path:
    if not isinstance(record, dict):
        raise PanelRunError(f"{label} is not a file record")
    path = Path(str(record.get("path", ""))).resolve()
    expected = record.get("sha256")
    if not path.is_file() or not isinstance(expected, str):
        raise PanelRunError(f"{label} file/hash is missing")
    observed = sha256_file(path)
    if observed != expected:
        raise PanelRunError(f"{label} SHA-256 mismatch")
    return path


def allocation_cpus() -> int:
    if not os.environ.get("SLURM_JOB_ID"):
        raise PanelRunError("refusing ORCA execution outside a SLURM allocation")
    raw = os.environ.get("SLURM_CPUS_ON_NODE", "")
    match = re.match(r"^(\d+)", raw)
    if match is None or int(match.group(1)) < 2:
        raise PanelRunError(f"invalid SLURM_CPUS_ON_NODE={raw!r}")
    return int(match.group(1))


def choose_parallelism(total_cpus: int, target_count: int) -> tuple[int, int]:
    # Each target runner executes its La and Ca legs concurrently.  Favor up to
    # 16 MPI ranks per leg, then add a target slot when that uses the otherwise
    # stranded remainder of a full-node allocation.
    target_slots = min(target_count, max(1, math.ceil(total_cpus / 32)))
    nprocs = min(16, total_cpus // (2 * target_slots))
    if nprocs < 1:
        raise PanelRunError("allocation cannot support a two-leg target runner")
    return target_slots, nprocs


def run_one_target(
    *,
    panel_id: str,
    manifest_path: Path,
    runner: Path,
    runner_hash: str,
    orca: Path,
    nprocs: int,
) -> dict[str, Any]:
    log_path = manifest_path.parent / "orca_runner.log"
    command = [
        str(Path(sys.executable).resolve()),
        str(runner),
        str(manifest_path),
        "--orca", str(orca),
        "--workers", "2",
        "--nprocs", str(nprocs),
        "--expected-runner-sha256", runner_hash,
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
        "panel_id": panel_id,
        "returncode": completed.returncode,
        "status": "complete" if completed.returncode == 0 else "failed",
        "manifest": {"path": str(manifest_path), "sha256": sha256_file(manifest_path)},
        "log": {"path": str(log_path), "sha256": sha256_file(log_path)},
    }


def run_panel(preparation_path: Path, pins_path: Path) -> Path:
    total_cpus = allocation_cpus()
    preparation_path = preparation_path.resolve()
    pins_path = pins_path.resolve()
    preparation = read_object(preparation_path)
    pins = read_object(pins_path)
    if (
        preparation.get("schema_version") != "alchemical_bvs.pqq_fixed_core_preparation.v1"
        or preparation.get("protocol_id") != PROTOCOL_ID
        or preparation.get("status") != "ready_for_orca"
        or preparation.get("target_count") != 25
        or preparation.get("task_count") != 50
        or preparation.get("protonation_subprotocol", {}).get("id")
        != "pdbfixer_standard_only_rng20260914_openmm_cpu_threads1_v1"
    ):
        raise PanelRunError("preparation manifest is not the complete frozen panel")
    if pins.get("protocol_id") != PROTOCOL_ID:
        raise PanelRunError("implementation pins carry another protocol ID")
    experiment_helpers = pins.get("experiment_helpers")
    if not isinstance(experiment_helpers, dict):
        raise PanelRunError("implementation pins lack experiment helpers")
    for name, record in experiment_helpers.items():
        verify_file_record(record, f"experiment helper {name}")
    own_record = experiment_helpers.get("run_panel")
    if (
        not isinstance(own_record, dict)
        or Path(str(own_record.get("path", ""))).resolve() != Path(__file__).resolve()
    ):
        raise PanelRunError("panel runner is not self-bound in implementation pins")
    pins_record = preparation.get("implementation_pins")
    if not isinstance(pins_record, dict) or (
        Path(str(pins_record.get("path", ""))).resolve() != pins_path
        or pins_record.get("sha256") != sha256_file(pins_path)
    ):
        raise PanelRunError("preparation is not bound to the supplied implementation pins")
    runner_record = pins.get("canonical_helpers", {}).get("run_orca_task_manifest")
    runner = verify_file_record(runner_record, "ORCA task runner")
    runner_hash = str(runner_record["sha256"])
    orca_record = pins.get("orca_runtime", {}).get("executable")
    orca = verify_file_record(orca_record, "ORCA executable")
    targets = preparation.get("targets")
    if not isinstance(targets, list) or len(targets) != 25:
        raise PanelRunError("preparation target list is incomplete")
    jobs: list[tuple[str, Path]] = []
    for target in targets:
        if not isinstance(target, dict) or not isinstance(target.get("panel_id"), str):
            raise PanelRunError("malformed preparation target")
        manifest = verify_file_record(target.get("manifest"), f"{target.get('panel_id')} carve manifest")
        jobs.append((target["panel_id"], manifest))

    target_slots, nprocs = choose_parallelism(total_cpus, len(jobs))
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=target_slots) as executor:
        futures = {
            executor.submit(
                run_one_target,
                panel_id=panel_id,
                manifest_path=manifest,
                runner=runner,
                runner_hash=runner_hash,
                orca=orca,
                nprocs=nprocs,
            ): panel_id
            for panel_id, manifest in jobs
        }
        for future in as_completed(futures):
            panel_id = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "panel_id": panel_id,
                    "returncode": None,
                    "status": "failed",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            results.append(result)
            print(f"{result['panel_id']}: {result['status']}", flush=True)

    failed = sorted(result["panel_id"] for result in results if result["status"] != "complete")
    job_id = str(os.environ["SLURM_JOB_ID"])
    receipt_path = preparation_path.parent / f"panel_execution_{job_id}.json"
    if receipt_path.exists():
        raise PanelRunError(f"refusing to overwrite execution receipt {receipt_path}")
    receipt = {
        "schema_version": EXECUTION_SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "started_at_utc": started,
        "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "complete" if not failed else "failed",
        "preparation": {"path": str(preparation_path), "sha256": sha256_file(preparation_path)},
        "implementation_pins": {"path": str(pins_path), "sha256": sha256_file(pins_path)},
        "allocation": {
            "slurm_job_id": job_id,
            "slurm_job_partition": os.environ.get("SLURM_JOB_PARTITION"),
            "host": os.uname().nodename,
            "cpus_on_node": total_cpus,
        },
        "parallelism": {
            "concurrent_targets": target_slots,
            "concurrent_legs_per_target": 2,
            "mpi_ranks_per_leg": nprocs,
            "maximum_concurrent_mpi_ranks": target_slots * 2 * nprocs,
        },
        "failed_panel_ids": failed,
        "targets": sorted(results, key=lambda item: item["panel_id"]),
    }
    write_json_atomic(receipt_path, receipt)
    if failed:
        raise PanelRunError(f"ORCA failures in targets: {failed}")
    return receipt_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--preparation", type=Path, default=HERE / "prepared" / "preparation.json"
    )
    parser.add_argument(
        "--implementation-pins", type=Path, default=HERE / "implementation_pins.json"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        receipt = run_panel(args.preparation, args.implementation_pins)
    except (PanelRunError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
