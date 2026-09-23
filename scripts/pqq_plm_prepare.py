#!/usr/bin/env python3
"""Package explicit PLM triples for the existing motion-envelope executor.

Request validates inputs without protonation. Prepare uses the released source
preparer, then the unchanged envelope and scoring planners. No energy is run here.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import re
import time

from affordable_common import InvalidArtifact, read_json, record, verify, write_new
import pqq_fast_prepare as source
import pqq_three_source_envelope as envelope
import pqq_three_source_envelope_execution as execution

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "params/pqq_plm_envelope_v1.json"
SCHEMA = "Nikasha_PLM_three_source_operations_v1"


def package(path):
    p = read_json(path)
    if p["package_id"] != "Nikasha_PLM_three_source_envelope_operations_v1" or p["scientific_protocol_changed"]:
        raise InvalidArtifact("unsupported operations package")
    for key in ("release", "reference", "strict_qualification", "rank_qualification", "maxiter_qualification"):
        verify(p[key])
    release = envelope.original.standard.release(verify(p["release"]))
    envelope.reference_status(p["reference"], release["source_configuration"])
    return p, release


def workspace(path):
    out = Path(path).resolve()
    if "workspaces" not in out.parts:
        raise InvalidArtifact("compute products belong under workspaces")
    return out


def request(sources, protein_ids, agreement, output, package_path=PACKAGE, source_preparations=None):
    if len(sources) != len(protein_ids) or not sources or len(set(protein_ids)) != len(protein_ids):
        raise InvalidArtifact("one distinct protein ID per source request required")
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", p) for p in protein_ids):
        raise InvalidArtifact("safe protein IDs required")
    if source_preparations is not None and len(source_preparations) != len(sources):
        raise InvalidArtifact("one exact archive per source request required, or omit all for fresh preparation")
    p, release = package(package_path)
    archives = [record(x) for x in source_preparations] if source_preparations is not None else [None] * len(sources)
    out = workspace(output)
    out.mkdir(parents=True, exist_ok=False)
    write_new(out / "CONFIG.json", release["source_configuration"])
    groups = []
    for src, pid, archive in zip(sources, protein_ids, archives):
        target = out / "requests" / (pid + ".json")
        envelope.request(src, out / "CONFIG.json", verify(p["reference"]), pid, agreement, target,
                         release=verify(p["release"]))
        req = read_json(target)
        if archive:
            old = read_json(verify(archive))
            if old["config"] != req["config"]:
                raise InvalidArtifact("archived source configuration differs")
            for c in req["cases"]:
                if sum(envelope.original.matching_source(c, r.get("source", {})) for r in old["cases"]) != 1:
                    raise InvalidArtifact("exact archived source match unavailable: " + c["case_id"])
        groups.append({"protein_id": pid, "request": record(target), "source_archive": archive})
    ids = [c["case_id"] for g in groups for c in read_json(verify(g["request"]))["cases"]]
    if len(set(ids)) != len(ids):
        raise InvalidArtifact("source identifiers repeated across proteins")
    data = {"schema": SCHEMA, "package": record(package_path), "agreement": record(agreement),
            "groups": groups, "source_mode": "reuse-exact" if source_preparations is not None else "fresh",
            "source_denominator": len(ids), "new_energy_calls": 0,
            "implementation": {name: record(ROOT / "scripts" / name) for name in
                ("pqq_plm_prepare.py", "pqq_fast_prepare.py", "pqq_three_source_envelope.py",
                 "pqq_three_source_envelope_execution.py")}}
    write_new(out / "REQUEST.json", data)
    return dry_run(out / "REQUEST.json")


def checked(path):
    d = read_json(path)
    if d["schema"] != SCHEMA or d["new_energy_calls"] != 0:
        raise InvalidArtifact("unsupported preparation request")
    p, _ = package(verify(d["package"]))
    verify(d["agreement"])
    for pin in d["implementation"].values():
        verify(pin)
    for g in d["groups"]:
        req = envelope.check_request(read_json(verify(g["request"])))
        if req["protein_id"] != g["protein_id"] or req["reference"] != p["reference"]:
            raise InvalidArtifact("group/reference mismatch")
        if g["source_archive"]:
            verify(g["source_archive"])
    return d, p


def dry_run(request):
    d, p = checked(request)
    n = len(d["groups"])
    return {"request": record(request), "status": "input_preflight_passed", "proteins": n,
            "sources": d["source_denominator"], "source_mode": d["source_mode"],
            "new_energy_calls": 0, "new_protonation_calls": 0,
            "potential_execution_calls": execution.interface().declared_calls(3 * n),
            "prospective_source_protonations": 3 * n if d["source_mode"] == "fresh" else 0,
            "resources": p["resources"], "next_operation": "prepare"}


def prepare(request, output):
    d, p = checked(request)
    out = workspace(output)
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    rows, preparations = [], []
    for group in d["groups"]:
        pid = group["protein_id"]
        req = envelope.check_request(read_json(verify(group["request"])))
        target = out / pid
        target.mkdir()
        archive = verify(group["source_archive"]) if group["source_archive"] else target / "SOURCE_PREPARATION.json"
        if not group["source_archive"]:
            cases = []
            for c in req["cases"]:
                try:
                    row = source.prepare_case(c, req["config"], target / "source" / c["case_id"])
                except Exception as exc:
                    row = {"case_id": c["case_id"], "source": copy.deepcopy(c), "status": "unsupported",
                           "reason": str(exc), "exception_type": type(exc).__name__}
                cases.append(row)
            write_new(archive, {"config": req["config"], "cases": cases, "source_mode": "fresh",
                               "request": group["request"], "new_protonation": True, "new_energy_calls": 0})
        try:
            envelope.prepare(verify(group["request"]), archive, target / "envelope")
            prep = target / "envelope/PREPARATION.json"
            data = read_json(prep)
            row = {"protein_id": pid, "status": data["status"], "preparation": record(prep),
                   "missing_members": data["missing_members"], "source_preparation": record(archive),
                   "source_statuses": [{k: r.get(k) for k in ("case_id", "status", "reason")}
                                       for r in read_json(archive)["cases"]]}
            if data["status"] == "prepared":
                preparations.append(prep)
        except Exception as exc:
            row = {"protein_id": pid, "status": "unavailable", "reason": str(exc),
                   "exception_type": type(exc).__name__, "source_preparation": record(archive)}
        rows.append(row)
        print(json.dumps({k: row.get(k) for k in ("protein_id", "status", "reason")}), flush=True)
    summary = {"schema": SCHEMA, "request": record(request), "source_mode": d["source_mode"],
               "groups": rows, "denominator": len(rows), "available": len(preparations),
               "new_energy_calls": 0, "source_preparation_wall_seconds": time.monotonic() - started,
               "execution_plan": None, "status": "preparation_complete"}
    # Write failures even if downstream planning raises; incomplete triples never enter scoring.
    write_new(out / "PREPARATION.json", summary)
    if preparations:
        execution.prepare(preparations, verify(p["strict_qualification"]), verify(p["rank_qualification"]),
                          verify(p["maxiter_qualification"]), verify(d["agreement"]), out / "scoring")
        summary["execution_plan"] = record(out / "scoring/plan.json")
    summary["status"] = "ready" if preparations else "all_groups_unavailable"
    write_new(out / "READY.json", summary)
    return {"status": summary["status"], "proteins": len(rows), "prepared": len(preparations),
            "receipt": record(out / "READY.json"), "execution_plan": summary["execution_plan"], "new_energy_calls": 0}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="op", required=True)
    q = sub.add_parser("request")
    q.add_argument("--sources", nargs="+", required=True)
    q.add_argument("--protein-ids", nargs="+", required=True)
    q.add_argument("--source-preparations", nargs="+")
    q.add_argument("--package-path", default=str(PACKAGE))
    q.add_argument("--agreement", required=True)
    q.add_argument("--output", required=True)
    q = sub.add_parser("prepare")
    q.add_argument("--request", required=True)
    q.add_argument("--output", required=True)
    q = sub.add_parser("dry-run")
    q.add_argument("--request", required=True)
    args = vars(ap.parse_args())
    op = args.pop("op").replace("-", "_")
    print(json.dumps(globals()[op](**args), indent=2))


if __name__ == "__main__":
    main()
