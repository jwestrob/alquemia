#!/usr/bin/env python3
"""Describe recorded baseline ordering without fitting or changing decisions."""
import argparse
import collections
import hashlib
import json
import math
import os
from pathlib import Path
import time


def pin(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def read(p):
    if pin(p["path"]) != p:
        raise ValueError("Changed artifact: " + p["path"])
    return json.loads(Path(p["path"]).read_text())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--config-sha", required=True)
    args = ap.parse_args()
    if not os.environ.get("SLURM_JOB_ID", "").isdigit():
        raise RuntimeError("Saved scientific aggregation requires Slurm")
    start = time.monotonic()
    cp = {"path": str(Path(args.config).resolve()), "sha256": args.config_sha}
    c = read(cp)
    assert pin(__file__) == c["implementation"]
    read(c["ledger"])
    data = read(c["collection"])
    rows = {r["record_id"]: r for r in data["rows"]}
    assert len(rows) == len(data["rows"])
    result = {"protocol": "saved_baseline_ordering_v1", "config": cp,
              "job_id": os.environ["SLURM_JOB_ID"], "cohorts": {},
              "new_molecular_calls": 0, "threshold_changed": False,
              "interpretation": "Consumed-cohort ordering only; no new calibration or decisions.260 memberships represent257 contexts with correlated replicas."}
    for name, ids in c["cohorts"].items():
        assert len(ids) == len(set(ids)) and set(ids) <= set(rows)
        rr = [rows[i] for i in ids]
        assert all(r["expected_class"] in ("Ca", "La") for r in rr)
        cohort = {"declared_memberships": len(rr), "policies": {}}
        for policy in ("operational", "mathematical", "origin"):
            available = []
            missing = []
            inconclusive = []
            cal = data["calibrations"][policy]
            assert cal["status"] == "available"
            bands = cal["bands"]
            for r in rr:
                pool = r.get("pool")
                cell = pool.get(policy) if isinstance(pool, dict) else None
                if cell is None:
                    missing.append(r["record_id"])
                    continue
                score = cell["R_model_kcal_mol"]
                assert math.isfinite(score)
                entry = {k: r[k] for k in ("record_id", "source_case_id", "protein", "biological_group", "expected_class", "source_conditioning_metal")}
                entry.update(score_kcal_mol=score, recorded_decision=r["decisions"][policy])
                available.append(entry)
                if r["decisions"][policy]["outcome"] == "inconclusive":
                    inconclusive.append({**entry, "above_Ca_band_kcal_mol": score-bands["Ca_max"], "below_La_band_kcal_mol": bands["La_min"]-score})
            ca = [v for v in available if v["expected_class"] == "Ca"]
            la = [v for v in available if v["expected_class"] == "La"]
            ca_max = max(ca, key=lambda v: v["score_kcal_mol"]) if ca else None
            la_min = min(la, key=lambda v: v["score_kcal_mol"]) if la else None
            gap = la_min["score_kcal_mol"]-ca_max["score_kcal_mol"] if ca and la else None
            cohort["policies"][policy] = {
                "available": len(available), "unavailable_record_ids": missing,
                "recorded_outcomes": dict(collections.Counter(r["decisions"][policy]["outcome"] for r in rr)),
                "fixed_calibration": cal, "maximum_Ca": ca_max, "minimum_La": la_min,
                "observed_class_gap_kcal_mol": gap, "inconclusive": inconclusive,
                "Ca_above_minimum_La": [v for v in ca if la_min and v["score_kcal_mol"] >= la_min["score_kcal_mol"]],
                "La_below_maximum_Ca": [v for v in la if ca_max and v["score_kcal_mol"] <= ca_max["score_kcal_mol"]]}
        result["cohorts"][name] = cohort
    result["elapsed_seconds"] = time.monotonic()-start
    out = Path(c["output"])
    out.mkdir(exist_ok=False)
    (out/"RESULT.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
