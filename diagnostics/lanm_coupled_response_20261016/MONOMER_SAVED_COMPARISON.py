#!/usr/bin/env python3
"""Combine existing Mex/LBT matrices only. --config CONFIG --output FRESH."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time


def pin(p):
    p = Path(p).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def read(p):
    assert pin(p["path"])["sha256"] == p["sha256"], p["path"]
    return json.loads(Path(p["path"]).read_text())


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--config", type=Path, required=True)
    a.add_argument("--output", type=Path, required=True)
    args = a.parse_args()
    assert os.environ.get("SLURM_JOB_ID", "").isdigit()
    assert int(os.environ["SLURM_CPUS_PER_TASK"]) == 1
    started = time.monotonic()
    cfg = json.loads(args.config.read_text())
    assert pin(__file__)["sha256"] == cfg["implementation"]["sha256"]
    assert pin(cfg["plan"]["path"])["sha256"] == cfg["plan"]["sha256"]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    mex, lbt = read(cfg["mex_response"]), read(cfg["lbt_pool"])
    assert mex["status"] == "passed" and mex["slurm_job_id"] == "12717"
    assert lbt["status"] == "complete" and lbt["slurm_job_id"] == "12720"
    assert len(mex["rows"]) == 4 and all(r["status"] == "complete" for r in mex["rows"])
    assert len(lbt["rows"]) == 36 and all(r["status"] == "complete" for r in lbt["rows"])
    lc = read(lbt["config"])
    origin = read(lc["matched_origin_result"])
    oc = read(origin["config"])
    assert oc["mex_result"] == mex["static_result"]
    read(mex["static_result"])
    pool = mex["comparison"]["common3"]
    assert pool["status"] == "available" and mex["comparison"]["finite_pilot_passed"]
    assert all(pool["selections"][m]["geometry"] == m + "_proposal" for m in ("La", "Dy"))
    geo = {r["proposer"]: r for r in mex["comparison"]["geometries"]}
    names = ("MACE_vacuum", "GFN2_field_difference", "MM_ledger", "fixed_FF_charge_ALPB")
    delta = {k: geo["Dy"]["component_work_kcal_mol"]["Dy"][k] -
                geo["La"]["component_work_kcal_mol"]["La"][k] for k in names}
    menus = lbt["comparison"]["menus"]
    expected = {(s, m) for s in ("7CCO", "7CCN", "global")
                for m in ("origin", "prior", "control", "candidate", "union")}
    assert {(r["source_scope"], r["menu"]) for r in menus} == expected and len(menus) == 15
    assert all(r["status"] == "available" for r in menus)
    rows = []
    for r in menus:
        assert math.fsum(r["R_components_kcal_mol"].values()) == r["R_kcal_mol"]
        changed = {k: math.fsum((r["R_components_kcal_mol"][k], delta[k])) for k in names}
        rows.append({"source_scope": r["source_scope"], "menu": r["menu"],
                     "LBT_selected": r["selected"], "Mex_selected": pool["selections"],
                     "R_static_Mex_kcal_mol": r["R_kcal_mol"],
                     "R_accommodated_Mex_kcal_mol": math.fsum(changed.values()),
                     "R_accommodated_components_kcal_mol": changed})
    by = {(r["source_scope"], r["menu"]): r for r in rows}
    spreads = {m: abs(by[("7CCN", m)]["R_static_Mex_kcal_mol"] -
                      by[("7CCO", m)]["R_static_Mex_kcal_mol"])
               for m in ("origin", "prior", "control", "candidate", "union")}
    result = {"status": "complete", "slurm_job_id": os.environ["SLURM_JOB_ID"],
              "implementation": pin(__file__), "config": pin(args.config), "plan": cfg["plan"],
              "parents": {"Mex": cfg["mex_response"], "LBT": cfg["lbt_pool"]},
              "new_molecular_calls": 0, "new_geometries": 0, "rows": rows,
              "Mex_accommodation_component_change_kcal_mol": delta,
              "Mex_accommodation_change_kcal_mol": math.fsum(delta.values()),
              "source_local_exchange_spread_kcal_mol": spreads,
              "elapsed_seconds": time.monotonic() - started}
    (out / "RESULT.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    lines = ["# Saved combined monomer response", "",
             "Positive R favors Mex-La + LBT-Dy; kcal/mol. No new molecular calls.", "",
             "| LBT scope | Menu | Static Mex | Accommodated Mex |",
             "|---|---|---:|---:|"]
    lines += [f"| {r['source_scope']} | {r['menu']} | {r['R_static_Mex_kcal_mol']:.6f} | {r['R_accommodated_Mex_kcal_mol']:.6f} |" for r in rows]
    lines += ["", "| LBT menu | Difference between source-local contrasts |", "|---|---:|"]
    lines += [f"| {m} | {v:.6f} |" for m, v in spreads.items()]
    lines += ["", "All sources and menus retained. Conditional fixed occupancy and chemistry; no affinity, population or independent-protein accuracy claim."]
    (out / "REPORT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
