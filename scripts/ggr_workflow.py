"""GGR task packaging and comparison using the established ORCA runner/collector."""
from pathlib import Path
import argparse
import json
import re
import shutil

from affordable_common import (InvalidArtifact, HA_TO_KCAL, cache_key, energy,
                               paired, read_json, record, verify, write_new, xyz)
from affordable_benchmark import audit_reference, collect
from affordable_benchmark_set import preserve_implementation


ALLOWED_PROTOCOLS = {
    "generic_peptide_amide_vertical_native_r2scan3c_v3",
    "generic_peptide_alpha_caps_native_r2scan3c_dev_v1",
    "ggr_connected_segment_native_r2scan3c_dev_v1",
}


def prepare(config_path, root, output, agreement):
    root, output = root.resolve(), output.resolve()
    if output.exists() or not output.is_relative_to(root / "workspaces"):
        raise InvalidArtifact("new workspace directory required")
    cfg = read_json(config_path)
    if cfg["stage"] not in ("A", "B") or not cfg["cases"]:
        raise InvalidArtifact("named stage A/B cases required")
    verify(cfg["approved_plan"])
    reference = audit_reference(root)
    snapshot = preserve_implementation(output / "implementation")
    extras = {}
    for path in sorted((root / "scripts").glob("ggr_*.py")):
        dest = output / "implementation" / path.name
        shutil.copyfile(path, dest)
        extras[path.name] = {"live_source": record(path), "preserved_copy": record(dest)}
    write_new(output / "implementation/ggr_inventory.json", extras)
    tasks, cases = [], []
    seen = set()
    for item in cfg["cases"]:
        name = item["case"]
        if not re.fullmatch(r"[A-Za-z0-9_]+", name) or name in seen:
            raise InvalidArtifact("invalid/duplicate case identifier")
        seen.add(name)
        mp = verify(item["preparation_manifest"])
        m = read_json(mp)
        if m["protocol_id"] not in ALLOWED_PROTOCOLS:
            raise InvalidArtifact("unapproved representation protocol")
        ids, copied, charges = [], {}, {}
        for metal in ("La", "Ca"):
            ep = m["outputs"][metal]
            ip = verify(ep.get("input", ep.get("orca_input")))
            xp = verify(ep["xyz"])
            body = ip.read_text()
            lines = [s.strip().lower().split()[1:] for s in body.splitlines() if s.strip().startswith("!")]
            if len(lines) != 1 or set(lines[0]) != {"r2scan-3c", "noautostart", "cpcm(water)", "defgrid3"}:
                raise InvalidArtifact("stage A/B requires the frozen normal-SCF native recipe")
            if re.search(r"%basis|%pointcharges|%pal|\b(?:Opt|NumGrad|EnGrad|Freq|TightSCF)\b", body, re.I):
                raise InvalidArtifact("unexpected method/resource override")
            coord = re.search(r"^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+(\S+)\s*$", body, re.M | re.I)
            if not coord or int(coord[2]) != 1 or Path(coord[3]).name != xp.name:
                raise InvalidArtifact("coordinate binding or multiplicity")
            charges[metal] = int(coord[1])
            if charges[metal] != ep["charge"] or ep["multiplicity"] != 1:
                raise InvalidArtifact("declared endpoint state mismatch")
            d = output / name / metal
            d.mkdir(parents=True)
            i, x = d / ip.name, d / xp.name
            shutil.copyfile(ip, i)
            shutil.copyfile(xp, x)
            copied[metal] = x
            tid = f"{name}_{metal}"
            ids.append(tid)
            tasks.append({"task_id": tid, "case": name, "lane": item["representation"],
                          "metal": metal, "charge": charges[metal], "multiplicity": 1,
                          "input": record(i), "xyz": record(x), "output_path": str(i.with_suffix(".out")),
                          "source_input": record(ip), "source_xyz": record(xp),
                          "source_manifest": record(mp), "protocol_id": m["protocol_id"],
                          "task_type": "single_point", "stage": cfg["stage"]})
        inv = paired(copied["La"], copied["Ca"], charges["La"], charges["Ca"])
        if charges != {"La": 0, "Ca": -1}:
            raise InvalidArtifact("approved comparisons require unchanged ligand charge -3")
        row = {k: item[k] for k in ("case", "target_id", "source_structure_id", "biological_group", "evidence_stratum")}
        row.update(lane=item["representation"], protocol_id=m["protocol_id"],
                   source_manifest=record(mp), tasks=ids, paired_invariants=inv,
                   observation_group=item["target_id"], evaluation_role="method_development_already_consumed",
                   descriptors={"coordination": m["coordination"], "core_charge_La": 0, "core_charge_Ca": -1},
                   explicit_water_inventory=m["explicit_water_inventory"], stage=cfg["stage"],
                   historical_case=item.get("historical_case"),
                   atom_count=len(xyz(copied["La"])))
        row["cache_key"] = cache_key(row)
        cases.append(row)
    result = {"schema_version": "alquemia.ggr_tasks.v1", "status": "prepared_not_executed",
              "stage": cfg["stage"], "tasks": tasks, "cases": cases,
              "configuration": record(config_path), "agreement": record(agreement),
              "approved_plan": cfg["approved_plan"], "reference": reference,
              "release": record(root / "diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json"),
              "implementation": extras[Path(__file__).name]["preserved_copy"], "implementation_snapshot": snapshot,
              "ggr_implementation_snapshot": record(output / "implementation/ggr_inventory.json"),
              "execution_policy": {key: record(root / "scripts" / name) for key, name in (
                  ("task_runner", "run_orca_task_manifest.py"),
                  ("runtime_renderer", "render_orca_runtime_input.py"))},
              "orca": cfg["orca"], "compute_budget": None, "wall_time_limit": None,
              "decision_policy": "uncalibrated raw contrasts; existing aquo gauge only; no threshold fitting",
              "unsupported_preparations": cfg.get("unsupported_preparations", [])}
    verify(result["orca"])
    write_new(output / "manifest.json", result)
    return result


def checked_rows(path):
    result = read_json(path)
    verify(result["manifest"])
    for row in result["rows"]:
        if row["status"] != "complete":
            continue
        for ep in row["endpoints"].values():
            verify(ep["receipt"])
            if energy(verify(ep["output"])) != ep["energy_hartree"]:
                raise InvalidArtifact("archived endpoint energy mismatch")
        raw = row["endpoints"]["Ca"]["energy_hartree"] - row["endpoints"]["La"]["energy_hartree"]
        if abs(raw - row["score"]["R_hartree"]) > 1e-10:
            raise InvalidArtifact("archived contrast algebra mismatch")
    return result


def compare(collections, historical, output):
    prior = []
    for path in historical:
        prior.extend(checked_rows(path)["rows"])
    results = [checked_rows(path) for path in collections]
    rows = [row for result in results for row in result["rows"]]
    differences = []
    for row in rows:
        candidates = [p for p in prior if p["case"] == row.get("historical_case") and
                      p["protocol_id"] == "generic_peptide_amide_vertical_native_r2scan3c_v3"]
        if not candidates:
            continue
        if len(candidates) != 1:
            raise InvalidArtifact("ambiguous historical comparison")
        p = candidates[0]
        valid = row["status"] == p["status"] == "complete"
        differences.append({"case": row["case"], "historical_case": p["case"],
                            "status": "complete" if valid else "unavailable",
                            "delta_R_kcal_mol": (row["score"]["R_hartree"] - p["score"]["R_hartree"]) * HA_TO_KCAL if valid else None})
    nma_ggr = [r for r in rows if r.get("source_structure_id") == "1GLG" and r["lane"] == "nma"]
    ordering = []
    if len(nma_ggr) > 1:
        raise InvalidArtifact("duplicate primary NMA GGR")
    if nma_ggr:
        g = nma_ggr[0]
        for row in rows:
            if row["target_id"] != "ALACTA_BOVINE_STRONG_SITE" or row["lane"] != "nma":
                continue
            valid = row["status"] == g["status"] == "complete"
            ordering.append({"case": row["case"], "GGR_case": g["case"],
                             "status": "complete" if valid else "unavailable",
                             "alpha_minus_GGR_R_kcal_mol": (row["score"]["R_hartree"] - g["score"]["R_hartree"]) * HA_TO_KCAL if valid else None})
    result = {"schema_version": "alquemia.ggr_comparison.v1", "implementation": record(__file__),
              "collections": [record(p) for p in collections], "historical_collections": [record(p) for p in historical],
              "rows": rows, "representation_differences": differences, "alpha_minus_GGR": ordering,
              "completed_endpoints": sum(r["completed_endpoints"] for r in results),
              "total_endpoints": sum(r["total_endpoints"] for r in results),
              "default_changed": False, "decision": "uncalibrated_protocol",
              "response_status": "response_model_not_validated",
              "interpretation": "All cases are development; structural replicates are not independent observations. No sign-based selection or calibration."}
    write_new(output, result)
    lines = ["# GGR representation comparisons", "", "Existing aquo reporting gauge; no calibrated zero or inherited PQQ bands.", "",
             "| Case | Representation | Atoms | S, kcal/mol | Status |", "|---|---|---:|---:|---|"]
    for r in rows:
        value = "unavailable" if r["score"] is None else f'{r["score"]["S_kcal_mol"]:.6f}'
        lines.append(f'| {r["case"]} | {r["lane"]} | {r["atom_count"]} | {value} | {r["status"]} |')
    lines += ["", "## Change from the corresponding archived formamide model", ""]
    for d in differences:
        lines.append(f'- {d["case"]}: {d["delta_R_kcal_mol"]} kcal/mol; {d["status"]}.')
    lines += ["", "## Matched extended-amide alpha minus GGR ordering", ""]
    for d in ordering:
        lines.append(f'- {d["case"]}: {d["alpha_minus_GGR_R_kcal_mol"]} kcal/mol; {d["status"]}.')
    lines += ["", result["interpretation"], "", "Underlying endpoint receipts, units, raw contrasts and source mappings are retained in the JSON."]
    with output.with_suffix(".md").open("x") as f:
        f.write("\n".join(lines) + "\n")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="operation", required=True)
    q = sub.add_parser("prepare")
    for name in ("config", "root", "output", "agreement"):
        q.add_argument("--" + name, type=Path, required=True)
    q = sub.add_parser("collect")
    q.add_argument("--manifest", type=Path, required=True)
    q.add_argument("--output", type=Path, required=True)
    q = sub.add_parser("compare")
    q.add_argument("--collection", type=Path, action="append", required=True)
    q.add_argument("--historical-collection", type=Path, action="append", required=True)
    q.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.operation == "prepare":
        result = prepare(a.config, a.root, a.output, a.agreement)
    elif a.operation == "collect":
        result = collect(a.manifest.resolve())
        write_new(a.output, result)
    else:
        result = compare(a.collection, a.historical_collection, a.output)
    print(json.dumps({k: result[k] for k in ("status", "completed_endpoints", "total_endpoints") if k in result}))


if __name__ == "__main__":
    main()
