#!/usr/bin/env python3
"""Run a manifested ORCA task panel inside an existing SLURM allocation.

The runner refuses login-node execution, verifies every prepared input/XYZ hash,
adds allocation-specific ``%pal`` controls without mutating scientific inputs,
and writes a hash-linked execution record for every attempt.  It submits no job
and requests no resources; invoke it from a project-compliant allocation.

Supported manifests are the site-mechanics and closed-shell charge/radius task
lists, both generic and PQQ carve outputs, and the native CN8 aquo
``calculations`` mapping.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

from render_orca_runtime_input import render_runtime_input
from result_protocol import sha256_file


EXECUTION_SCHEMA = "alchemical_bvs.orca_task_execution.v1"
NORMAL_TERMINATION = "ORCA TERMINATED NORMALLY"
SCF_CONVERGED = "SCF CONVERGED AFTER"
SCF_FAILED = "SCF NOT CONVERGED"
PROGRAM_VERSION_RE = re.compile(r"Program Version\s+([0-9][0-9A-Za-z._-]*)")


class TaskRunnerError(ValueError):
    """A manifest or execution context violates the runner contract."""


def _utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def _write_json_atomic(path: Path, value: Any) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def execution_record_path(output_path: Path) -> Path:
    return output_path.with_name(output_path.name + ".execution.json")


def _runner_identity() -> dict[str, str]:
    path = Path(__file__).resolve()
    digest = sha256_file(path)
    if digest is None:
        raise TaskRunnerError(f"cannot hash manifested ORCA runner: {path}")
    return {"path": str(path), "sha256": digest}


def _runtime_renderer_identity() -> dict[str, str]:
    path = Path(__file__).resolve().with_name("render_orca_runtime_input.py")
    digest = sha256_file(path)
    if digest is None:
        raise TaskRunnerError(f"cannot hash ORCA runtime renderer: {path}")
    return {"path": str(path), "sha256": digest}


def _verify_runner_policy(
    manifest: dict[str, Any], expected_runner_sha256: str | None
) -> dict[str, str]:
    """Bind execution to the runner implementation selected at preparation."""
    actual = _runner_identity()
    policy = manifest.get("execution_policy")
    implementation = manifest.get("implementation")
    declared_records: list[tuple[str, Any]] = []
    if isinstance(policy, dict) and "task_runner" in policy:
        declared_records.append(
            ("execution_policy.task_runner", policy["task_runner"])
        )
    if isinstance(implementation, dict) and "execution_runner" in implementation:
        declared_records.append(
            ("implementation.execution_runner", implementation["execution_runner"])
        )
    declared_hashes: set[str] = set()
    for label, declared in declared_records:
        if not isinstance(declared, dict):
            raise TaskRunnerError(f"manifest {label} is not an object")
        declared_hash = declared.get("sha256")
        if (
            not isinstance(declared_hash, str)
            or re.fullmatch(r"[0-9a-f]{64}", declared_hash) is None
        ):
            raise TaskRunnerError(f"manifest {label} SHA-256 is invalid")
        declared_hashes.add(declared_hash)
        if declared_hash != actual["sha256"]:
            raise TaskRunnerError(
                "manifest was prepared for a different ORCA task-runner implementation"
            )
    if len(declared_hashes) > 1:
        raise TaskRunnerError("manifest declares conflicting task-runner SHA-256 values")
    if expected_runner_sha256 is not None:
        if re.fullmatch(r"[0-9a-f]{64}", expected_runner_sha256) is None:
            raise TaskRunnerError("submit-script task-runner SHA-256 is invalid")
        if expected_runner_sha256 != actual["sha256"]:
            raise TaskRunnerError(
                "submit-script task-runner SHA-256 does not match the executing runner"
            )
        if declared_hashes and expected_runner_sha256 not in declared_hashes:
            raise TaskRunnerError(
                "submit-script and manifest task-runner SHA-256 values disagree"
            )
    return actual


def _verify_runtime_renderer_policy(
    manifest: dict[str, Any], actual: dict[str, str]
) -> None:
    policy = manifest.get("execution_policy")
    if not isinstance(policy, dict) or "runtime_renderer" not in policy:
        return
    declared = policy["runtime_renderer"]
    if not isinstance(declared, dict):
        raise TaskRunnerError("manifest execution_policy.runtime_renderer is not an object")
    declared_hash = declared.get("sha256")
    if (
        not isinstance(declared_hash, str)
        or re.fullmatch(r"[0-9a-f]{64}", declared_hash) is None
    ):
        raise TaskRunnerError(
            "manifest execution_policy.runtime_renderer SHA-256 is invalid"
        )
    if declared_hash != actual["sha256"]:
        raise TaskRunnerError(
            "manifest was prepared for a different ORCA runtime-renderer implementation"
        )


def _resolve_under(root: Path, value: str, label: str) -> Path:
    path = Path(value)
    resolved = (path if path.is_absolute() else root / path).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise TaskRunnerError(f"{label} escapes manifest directory: {value}") from exc
    return resolved


def _normalise_task_list(manifest: dict[str, Any], root: Path) -> list[dict[str, Any]]:
    raw_tasks = manifest.get("tasks")
    tasks: list[dict[str, Any]] = []
    if isinstance(raw_tasks, list):
        for raw in raw_tasks:
            if not isinstance(raw, dict):
                raise TaskRunnerError("task list contains a non-object record")
            task_id = raw.get("task_id") or raw.get("label")
            input_record = raw.get("input")
            if isinstance(input_record, dict):
                input_value = input_record.get("path")
                input_hash = input_record.get("sha256")
                xyz_record = raw.get("xyz")
                xyz_value = xyz_record.get("path") if isinstance(xyz_record, dict) else None
                xyz_hash = xyz_record.get("sha256") if isinstance(xyz_record, dict) else None
            else:
                input_value = raw.get("input_path")
                input_hash = raw.get("input_sha256")
                xyz_value = raw.get("xyz_path")
                xyz_hash = raw.get("xyz_sha256")
            output_value = raw.get("output_path")
            if not all(isinstance(value, str) and value for value in (task_id, input_value, input_hash, xyz_value, xyz_hash, output_value)):
                raise TaskRunnerError(f"task has incomplete paths/hashes: {task_id!r}")
            engrad_value = raw.get("engrad_path")
            tasks.append(
                {
                    "task_id": task_id,
                    "input": _resolve_under(root, input_value, f"{task_id} input"),
                    "input_sha256": input_hash,
                    "xyz": _resolve_under(root, xyz_value, f"{task_id} XYZ"),
                    "xyz_sha256": xyz_hash,
                    "output": _resolve_under(root, output_value, f"{task_id} output"),
                    "engrad": (
                        _resolve_under(root, engrad_value, f"{task_id} engrad")
                        if isinstance(engrad_value, str) and engrad_value
                        else None
                    ),
                }
            )
        return tasks

    calculations = manifest.get("calculations")
    if not isinstance(calculations, dict):
        outputs = manifest.get("outputs")
        stem = manifest.get("stem")
        if not isinstance(outputs, dict) or not isinstance(stem, str) or not stem:
            raise TaskRunnerError(
                "manifest has no supported task list, calculations mapping, or carve outputs"
            )
        for label in ("La", "Ca", "apo", "water"):
            nested = outputs.get(label)
            if isinstance(nested, dict):
                input_record = nested.get("orca_input")
                xyz_record = nested.get("xyz")
            else:
                input_record = outputs.get(f"{label}_input")
                xyz_record = outputs.get(f"{label}_xyz")
            if not isinstance(input_record, dict) or not isinstance(xyz_record, dict):
                raise TaskRunnerError(f"carve manifest lacks {label} input/XYZ records")
            input_value = input_record.get("path")
            input_hash = input_record.get("sha256")
            xyz_value = xyz_record.get("path")
            xyz_hash = xyz_record.get("sha256")
            if not all(
                isinstance(value, str) and value
                for value in (input_value, input_hash, xyz_value, xyz_hash)
            ):
                raise TaskRunnerError(f"carve manifest has incomplete {label} provenance")
            tasks.append(
                {
                    "task_id": label,
                    "input": _resolve_under(root, input_value, f"{label} input"),
                    "input_sha256": input_hash,
                    "xyz": _resolve_under(root, xyz_value, f"{label} XYZ"),
                    "xyz_sha256": xyz_hash,
                    "output": _resolve_under(
                        root, f"sp_{stem}_{label}.out", f"{label} output"
                    ),
                    "engrad": None,
                    "nprocs_override": 1 if label == "water" else None,
                }
            )
        return tasks
    for task_id, raw in calculations.items():
        if not isinstance(task_id, str) or not isinstance(raw, dict):
            raise TaskRunnerError("calculations mapping contains an invalid record")
        values = {
            "input": raw.get("input"),
            "input_sha256": raw.get("input_sha256"),
            "xyz": raw.get("xyz"),
            "xyz_sha256": raw.get("xyz_sha256"),
            "output": raw.get("expected_output"),
        }
        if not all(isinstance(value, str) and value for value in values.values()):
            raise TaskRunnerError(f"calculation has incomplete paths/hashes: {task_id}")
        tasks.append(
            {
                "task_id": task_id,
                "input": _resolve_under(root, values["input"], f"{task_id} input"),
                "input_sha256": values["input_sha256"],
                "xyz": _resolve_under(root, values["xyz"], f"{task_id} XYZ"),
                "xyz_sha256": values["xyz_sha256"],
                "output": _resolve_under(root, values["output"], f"{task_id} output"),
                "engrad": None,
                "nprocs_override": None,
            }
        )
    return tasks


def load_manifest_tasks(manifest_path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    manifest_path = manifest_path.resolve()
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise TaskRunnerError(f"cannot read task manifest: {manifest_path}") from exc
    if not isinstance(manifest, dict):
        raise TaskRunnerError("task manifest is not a JSON object")
    tasks = _normalise_task_list(manifest, manifest_path.parent)
    task_ids = [task["task_id"] for task in tasks]
    if len(set(task_ids)) != len(task_ids):
        raise TaskRunnerError("manifest task IDs are not unique")
    if not tasks:
        raise TaskRunnerError("manifest has no tasks")
    return manifest, tasks


def _parse_allocation() -> int:
    if not os.environ.get("SLURM_JOB_ID"):
        raise TaskRunnerError(
            "no SLURM_JOB_ID: refusing to execute ORCA outside an allocation"
        )
    raw = os.environ.get("SLURM_CPUS_ON_NODE", "")
    match = re.match(r"^(\d+)", raw)
    if match is None or int(match.group(1)) < 1:
        raise TaskRunnerError(f"invalid SLURM_CPUS_ON_NODE={raw!r}")
    return int(match.group(1))


def _resolve_orca(value: Path | None) -> Path:
    candidate: str | None
    if value is not None:
        candidate = str(value)
    elif os.environ.get("ORCA_PATH"):
        candidate = str(Path(os.environ["ORCA_PATH"]) / "orca")
    else:
        candidate = shutil.which("orca")
    if not candidate:
        raise TaskRunnerError("ORCA executable not found; set ORCA_PATH or --orca")
    path = Path(candidate).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise TaskRunnerError(f"ORCA executable is not an executable file: {path}")
    return path


def _verify_prepared_task(task: dict[str, Any]) -> None:
    if sha256_file(task["input"]) != task["input_sha256"]:
        raise TaskRunnerError(f"input SHA-256 mismatch: {task['task_id']}")
    if sha256_file(task["xyz"]) != task["xyz_sha256"]:
        raise TaskRunnerError(f"XYZ SHA-256 mismatch: {task['task_id']}")


def _completed_attempt_is_valid(
    record_path: Path,
    output_path: Path,
    *,
    manifest_sha256: str,
    task: dict[str, Any],
    runner_identity: dict[str, str],
    runtime_renderer_identity: dict[str, str],
) -> bool:
    try:
        record = json.loads(record_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    base_valid = bool(
        record.get("schema_version") == EXECUTION_SCHEMA
        and record.get("returncode") == 0
        and record.get("normal_termination") is True
        and record.get("scf_converged") is True
        and isinstance(record.get("orca_version"), str)
        and bool(record.get("orca_version"))
        and record.get("manifest", {}).get("sha256") == manifest_sha256
        and record.get("artifacts", {}).get("template_input", {}).get("sha256")
        == task["input_sha256"]
        and record.get("artifacts", {}).get("xyz", {}).get("sha256")
        == task["xyz_sha256"]
        and record.get("task_runner") == runner_identity
        and record.get("runtime_renderer") == runtime_renderer_identity
        and record.get("artifacts", {}).get("output", {}).get("sha256")
        == sha256_file(output_path)
    )
    if not base_valid:
        return False
    expected_engrad = task.get("engrad")
    if expected_engrad is None:
        return True
    observed_engrad_hash = sha256_file(expected_engrad)
    return bool(
        observed_engrad_hash
        and record.get("artifacts", {}).get("engrad", {}).get("sha256")
        == observed_engrad_hash
    )


def _run_one(
    task: dict[str, Any],
    *,
    manifest_path: Path,
    manifest_sha256: str,
    orca: Path,
    nprocs: int,
    runner_identity: dict[str, str],
    runtime_renderer_identity: dict[str, str],
) -> dict[str, Any]:
    _verify_prepared_task(task)
    output_path: Path = task["output"]
    record_path = execution_record_path(output_path)
    if output_path.exists() or record_path.exists():
        if output_path.is_file() and record_path.is_file() and _completed_attempt_is_valid(
            record_path,
            output_path,
            manifest_sha256=manifest_sha256,
            task=task,
            runner_identity=runner_identity,
            runtime_renderer_identity=runtime_renderer_identity,
        ):
            return {"task_id": task["task_id"], "status": "already_complete"}
        raise TaskRunnerError(
            f"refusing to overwrite an incomplete/unverified attempt for {task['task_id']}; "
            "use a fresh prepared directory"
        )

    runtime_path = task["input"].with_name(task["input"].stem + ".runtime.inp")
    runtime_sidecar = runtime_path.with_suffix(".json")
    if runtime_path.exists() or runtime_sidecar.exists():
        raise TaskRunnerError(
            f"stale runtime input exists for {task['task_id']}; use a fresh directory"
        )
    if task["engrad"] is not None and task["engrad"].exists():
        raise TaskRunnerError(f"stale engrad exists for {task['task_id']}")

    runtime_record = render_runtime_input(
        task["input"], runtime_path, nprocs=nprocs, sidecar_path=runtime_sidecar
    )
    # ORCA forwards one quoted argument string to mpirun (ORCA 6.1, section 2.5).
    command = [str(orca), runtime_path.name, "--bind-to none"]
    started = _utc_now()
    returncode: int | None = None
    error: str | None = None
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "OMPI_MCA_hwloc_base_binding_policy": "none",
            "ORCA_TMPDIR": str(runtime_path.parent),
            "TMPDIR": str(runtime_path.parent),
        }
    )
    try:
        with output_path.open("xb") as output_handle:
            completed = subprocess.run(
                command,
                cwd=runtime_path.parent,
                stdout=output_handle,
                stderr=subprocess.STDOUT,
                env=environment,
                check=False,
            )
        returncode = int(completed.returncode)
    except Exception as exc:  # the execution record must survive launch failures
        error = f"{type(exc).__name__}: {exc}"

    runtime_engrad = runtime_path.with_suffix(".engrad")
    if task["engrad"] is not None and runtime_engrad.is_file():
        runtime_engrad.replace(task["engrad"])
    output_text = output_path.read_text(errors="ignore") if output_path.is_file() else ""
    normal = NORMAL_TERMINATION in output_text
    scf_converged = SCF_CONVERGED in output_text and SCF_FAILED not in output_text
    versions = sorted(set(PROGRAM_VERSION_RE.findall(output_text)))
    orca_version = versions[0] if len(versions) == 1 else None
    artifacts: dict[str, Any] = {
        "template_input": {
            "path": str(task["input"]),
            "sha256": sha256_file(task["input"]),
        },
        "xyz": {"path": str(task["xyz"]), "sha256": sha256_file(task["xyz"])},
        "runtime_input": {
            "path": str(runtime_path),
            "sha256": sha256_file(runtime_path),
        },
        "runtime_sidecar": {
            "path": str(runtime_sidecar),
            "sha256": sha256_file(runtime_sidecar),
        },
        "output": {"path": str(output_path), "sha256": sha256_file(output_path)},
    }
    if task["engrad"] is not None:
        artifacts["engrad"] = {
            "path": str(task["engrad"]),
            "sha256": sha256_file(task["engrad"]),
        }
    expected_engrad_complete = bool(
        task["engrad"] is None
        or artifacts.get("engrad", {}).get("sha256")
    )
    record = {
        "schema_version": EXECUTION_SCHEMA,
        "task_id": task["task_id"],
        "started_at_utc": started,
        "finished_at_utc": _utc_now(),
        "manifest": {"path": str(manifest_path), "sha256": manifest_sha256},
        "orca_executable": {
            "path": str(orca),
            "sha256": sha256_file(orca),
        },
        "task_runner": runner_identity,
        "runtime_renderer": runtime_renderer_identity,
        "command": command,
        "parallelism": {
            "nprocs": nprocs,
            "omp_threads_per_rank": 1,
            "binding": "none",
        },
        "allocation": runtime_record["execution_context"],
        "returncode": returncode,
        "launch_error": error,
        "orca_version": orca_version,
        "normal_termination": normal,
        "scf_converged": scf_converged,
        "expected_engrad_complete": expected_engrad_complete,
        "artifacts": artifacts,
    }
    _write_json_atomic(record_path, record)
    status = (
        "complete"
        if (
            returncode == 0
            and normal
            and scf_converged
            and orca_version is not None
            and expected_engrad_complete
        )
        else "failed"
    )
    return {"task_id": task["task_id"], "status": status, "record": str(record_path)}


def run_manifest(
    manifest_path: Path,
    *,
    orca_path: Path | None = None,
    workers: int | None = None,
    nprocs: int | None = None,
    selected_tasks: set[str] | None = None,
    expected_runner_sha256: str | None = None,
) -> list[dict[str, Any]]:
    total_cpus = _parse_allocation()
    orca = _resolve_orca(orca_path)
    manifest_path = manifest_path.resolve()
    manifest, tasks = load_manifest_tasks(manifest_path)
    runner_identity = _verify_runner_policy(manifest, expected_runner_sha256)
    runtime_renderer_identity = _runtime_renderer_identity()
    _verify_runtime_renderer_policy(manifest, runtime_renderer_identity)
    if selected_tasks:
        known = {task["task_id"] for task in tasks}
        missing = selected_tasks - known
        if missing:
            raise TaskRunnerError(f"unknown selected tasks: {sorted(missing)}")
        tasks = [task for task in tasks if task["task_id"] in selected_tasks]

    if workers is None and nprocs is None:
        workers = min(3, total_cpus)
        nprocs = min(16, total_cpus // workers)
    elif workers is None:
        if nprocs is None or nprocs < 1:
            raise TaskRunnerError("nprocs must be positive")
        workers = min(3, total_cpus // nprocs)
    elif nprocs is None:
        if workers < 1:
            raise TaskRunnerError("workers must be positive")
        nprocs = min(16, total_cpus // workers)
    if workers is None or nprocs is None or workers < 1 or nprocs < 1:
        raise TaskRunnerError("allocation is too small for requested parallelism")
    if workers * nprocs > total_cpus:
        raise TaskRunnerError(
            f"workers*nprocs={workers*nprocs} exceeds allocated CPUs={total_cpus}"
        )

    manifest_hash = sha256_file(manifest_path)
    if manifest_hash is None:
        raise TaskRunnerError("manifest disappeared before execution")
    results: list[dict[str, Any]] = []
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                _run_one,
                task,
                manifest_path=manifest_path,
                manifest_sha256=manifest_hash,
                orca=orca,
                nprocs=task.get("nprocs_override") or nprocs,
                runner_identity=runner_identity,
                runtime_renderer_identity=runtime_renderer_identity,
            ): task["task_id"]
            for task in tasks
        }
        for future in as_completed(futures):
            task_id = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "task_id": task_id,
                    "status": "runner_error",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            results.append(result)
            print(f"{task_id}: {result['status']}", flush=True)
            if result["status"] not in {"complete", "already_complete"}:
                failures.append(task_id)
    if failures:
        raise TaskRunnerError(f"ORCA task failures: {sorted(failures)}")
    return sorted(results, key=lambda item: item["task_id"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--orca", type=Path)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--nprocs", type=int)
    parser.add_argument("--task", action="append", dest="tasks")
    parser.add_argument(
        "--expected-runner-sha256",
        help="fail unless this exact runner implementation is executing",
    )
    args = parser.parse_args()
    try:
        results = run_manifest(
            args.manifest,
            orca_path=args.orca,
            workers=args.workers,
            nprocs=args.nprocs,
            selected_tasks=set(args.tasks or []),
            expected_runner_sha256=args.expected_runner_sha256,
        )
    except TaskRunnerError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"completed/resumed {len(results)} manifested ORCA tasks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
