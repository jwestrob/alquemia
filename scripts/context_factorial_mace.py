"""Pinned, energy-only CPU replay of the two-boundary context factorial.

Login execution is explicitly authorized on 2026-10-12, within root's eight
threads of the aggregate sixteen-thread allowance. Historical GPU results are
reused only after coordinate, state, checkpoint and receipt checks. No optimizer.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import resource
import sys
import time
import traceback

CHECKPOINT_SHA = "9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469"
EV_TO_KCAL = 23.06054783061903
PROTOCOL = "context_factorial_native_omol_cpu_v1"


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": sha(path)}


def verify(ref):
    if sha(ref["path"]) != ref["sha256"]:
        raise ValueError(f"changed artifact: {ref['path']}")
    return Path(ref["path"])


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    with Path(path).open("x") as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write("\n")


def xyz(path):
    lines = Path(path).read_text().splitlines()
    rows = [(v[0], *(float(x) for x in v[1:]))
            for line in lines[2:] if (v := line.split())]
    if len(rows) != int(lines[0]) or any(len(r) != 4 for r in rows):
        raise ValueError("invalid XYZ")
    return rows


def key(row, metal):
    return "__".join((row["protein"], row["target_context"], row["geometry_id"], metal))


def reuse_cell(cell, task, model):
    receipt = read(verify(cell["MACE"]))
    archived_model = receipt.get("model")
    if archived_model is None:
        archived_model = read(verify(receipt["manifest"]))["model"]
    if archived_model != model:
        raise ValueError("archived model differs")
    state = receipt.get("state", receipt.get("input_state_check"))
    if (state["charge"] != task["charge"] or
            state["spin_multiplicity"] != task["multiplicity"] or
            state["head"] != "omol" or state["graph_count"] != 1):
        raise ValueError("archived physical state differs")
    a, b = xyz(verify(task["xyz"])), xyz(verify(cell["xyz"]))
    if len(a) != len(b) or [r[0] for r in a] != [r[0] for r in b]:
        raise ValueError("archived atom order differs")
    error = max(abs(x-y) for u, v in zip(a, b) for x, y in zip(u[1:], v[1:]))
    # Only decimal serialization noise; this is much tighter than preparation's
    # 1e-7-A motion mapping tolerance and is not structural approximation reuse.
    if error > 1e-12 or state["atoms"] != len(a):
        raise ValueError("archived numerical coordinates differ")
    energy = receipt["energy_eV"]
    if (receipt["status"] not in ("complete", "computed") or
            not math.isfinite(energy) or energy != cell["components"]["MACE_eV"]):
        raise ValueError("invalid archived energy receipt")
    return {"receipt": cell["MACE"], "xyz": cell["xyz"], "energy_eV": energy,
            "numeric_coordinate_max_error_A": error, "state": state}


def prepare(args):
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    prep = read(args.preparation)
    contexts = {}
    model = None
    for case in prep["cases"]:
        for context, data in case["contexts"].items():
            mp = verify(data["manifest"])
            m = read(mp)
            if model is None:
                model = m["model"]
            if m["model"] != model:
                raise ValueError("different archived model definitions")
            candidates = [mp.parent / n for n in ("collection_final.json", "COLLECTION_FINAL.json")]
            paths = [p for p in candidates if p.exists()]
            if len(paths) != 1:
                raise ValueError("ambiguous collection")
            collection = read(paths[0])
            matches = [c for c in collection["cases"] if c["case_id"].startswith(case["protein"])]
            if len(matches) != 1:
                raise ValueError("ambiguous source case")
            contexts[case["protein"], context] = (pin(paths[0]), matches[0])
    if (model["checkpoint"]["sha256"] != CHECKPOINT_SHA or
            model["head"] != "omol" or model["dtype"] != "float64" or
            model.get("execution_adapter") is not None or model.get("extra_dispersion") or
            model.get("charge_feature_adapter") or model["periodic"]):
        raise ValueError("only unchanged native vacuum OMOL is supported")
    verify(model["checkpoint"])
    tasks = []
    for row in prep["rows"]:
        for metal, endpoint in row["endpoints"].items():
            verify(endpoint["xyz"])
            t = {k: row[k] for k in ("protein", "target_context", "geometry_id")}
            t.update(endpoint, metal=metal, task_id=key(row, metal))
            t["reuse"] = None
            if row["previously_scored_same_context"]:
                cp, c = contexts[row["protein"], row["target_context"]]
                candidate = row["geometry_id"].removeprefix(row["target_context"] + "_")
                t["reuse"] = reuse_cell(c["matrix"][metal][candidate], t, model)
                t["reuse"]["collection"] = cp
            tasks.append(t)
    if len(tasks) != 40 or sum(t["reuse"] is not None for t in tasks) != 24:
        raise ValueError("unexpected factorial denominator")
    scouts = [t["task_id"] for t in tasks if t["protein"] == "q9z4j7-pqq-la_model" and t["geometry_id"] == "origin"]
    plan = {"protocol_id": PROTOCOL, "preparation": pin(args.preparation),
            "implementation": pin(__file__), "model": model, "tasks": tasks,
            "scouts": scouts, "threads": 8, "replay_tolerance_kcal_mol": 0.01,
            "coordinate_reuse_tolerance_A": 1e-12,
            "planned_new_calls": 20, "reused_main_cells": 24,
            "scope": "4 CPU/GPU replay checks; if all pass, 16 cross-context singlepoints; no optimization",
            "user_authorization": "2026-10-12 explicit <=16 total login numerical threads; root share <=8",
            "force_status": "not_requested", "calibration": None}
    write(out / "manifest.json", plan)
    print(json.dumps({"manifest": str(out / "manifest.json"), "tasks": 40, "reuse": 24, "new": 20}))


def execute(args):
    m = read(args.manifest)
    verify(m["implementation"])
    verify(m["preparation"])
    verify(m["model"]["checkpoint"])
    out = Path(args.manifest).resolve().parent
    if (out / "RUN_STARTED.json").exists():
        raise ValueError("run already attempted; collect preserved receipts instead of automatic retry")
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        if int(os.environ.get(name, "0")) < 1 or int(os.environ[name]) > m["threads"]:
            raise ValueError("thread environment missing or above declared bound: " + name)
    started = time.monotonic()
    cpu_start = time.process_time()
    write(out / "RUN_STARTED.json", {"pid": os.getpid(), "threads": m["threads"],
          "node": platform.node(), "python": pin(sys.executable), "manifest": pin(args.manifest),
          "thread_environment": {k: os.environ[k] for k in os.environ if k.endswith("NUM_THREADS")}})
    results = []
    try:
        import torch
        from ase import Atoms
        from mace.calculators import mace_omol
        torch.set_num_threads(m["threads"])
        torch.set_num_interop_threads(1)
        torch.set_default_dtype(torch.float64)
        calc = mace_omol(model=m["model"]["checkpoint"]["path"], device="cpu", default_dtype="float64")
        model = calc.models[0]
        if (type(model).__name__ != "ScaleShiftMACE" or list(model.heads) != ["omol"] or
                dict(model.embedding_specs) != m["model"]["embedding_specs"] or calc.energy_units_to_eV != 1.):
            raise ValueError("loaded model definition differs")
        versions = {n: p._version for n, p in model.named_parameters()}
        load_seconds = time.monotonic() - started
        software = {name: importlib.metadata.version(name) for name in ("torch", "mace-torch", "ase", "e3nn", "numpy")}
        jobs = [t for t in m["tasks"] if t["task_id"] in m["scouts"]]
        jobs += [t for t in m["tasks"] if not t["reuse"]]
        for t in jobs:
            rows = xyz(verify(t["xyz"]))
            atoms = Atoms([r[0] for r in rows], positions=[r[1:] for r in rows], pbc=False)
            atoms.info.update(charge=t["charge"], spin=t["multiplicity"])
            batch = calc._atoms_to_batch(atoms)
            if (batch["total_charge"].item() != t["charge"] or batch["total_spin"].item() != t["multiplicity"] or
                    len(batch["ptr"]) != 2 or calc.head != "omol"):
                raise ValueError("model batch state mismatch")
            tick = time.monotonic()
            with torch.no_grad():
                raw = model(batch.to_dict(), training=False, compute_force=False, compute_virials=False,
                            compute_stress=False, compute_displacement=False, compute_hessian=False,
                            compute_edge_forces=False, compute_atomic_stresses=False)
                energy = float(raw["energy"].item())
            if raw["forces"] is not None or not math.isfinite(energy):
                raise ValueError("invalid native energy-only result")
            if versions != {n: p._version for n, p in model.named_parameters()}:
                raise ValueError("weights mutated")
            r = {"task_id": t["task_id"], "status": "computed", "energy_eV": energy,
                 "model": m["model"], "xyz": t["xyz"], "charge": t["charge"], "multiplicity": t["multiplicity"],
                 "evaluation_seconds": time.monotonic()-tick, "versions": software,
                 "threads": torch.get_num_threads(), "forces": None, "parameter_versions_unchanged": True}
            if t["task_id"] in m["scouts"]:
                r["replay_error_kcal_mol"] = (energy-t["reuse"]["energy_eV"])*EV_TO_KCAL
                r["replay_pass"] = abs(r["replay_error_kcal_mol"]) <= m["replay_tolerance_kcal_mol"]
            write(out / (t["task_id"] + ".json"), r)
            results.append(r)
            print(json.dumps({k: r[k] for k in ("task_id", "energy_eV", "evaluation_seconds")}), flush=True)
            if r.get("replay_pass") is False:
                raise ValueError("CPU/GPU replay failed; no cross-context expansion")
        status, failure = "complete", None
    except Exception as exc:
        status, failure = "failed", f"{type(exc).__name__}: {exc}"
        traceback.print_exc()
        load_seconds = locals().get("load_seconds")
    final = {"protocol_id": PROTOCOL, "status": status, "failure": failure,
             "manifest": pin(args.manifest), "actual_new_calls": len(results), "results": results,
             "model_load_seconds": load_seconds, "wall_seconds": time.monotonic()-started,
             "process_cpu_seconds": time.process_time()-cpu_start,
             "peak_host_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
             "allocated_gpu_seconds": 0, "scheduler_allocation": None}
    write(out / "FINAL.json", final)
    if status != "complete":
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--preparation", required=True)
    p.add_argument("--output", required=True)
    p = sub.add_parser("execute")
    p.add_argument("--manifest", required=True)
    args = parser.parse_args()
    (prepare if args.command == "prepare" else execute)(args)
