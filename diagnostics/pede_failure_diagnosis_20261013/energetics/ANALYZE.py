#!/usr/bin/env python3
"""Finite saved-artifact decomposition; no calculators, fitting or geometry search."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import time


ORDER = ("origin", "adaptive_Ca", "adaptive_La")
METALS = ("Ca", "La")


def read(path):
    return json.loads(Path(path).read_text())


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def verified(record):
    if pin(record["path"])["sha256"] != record["sha256"]:
        raise ValueError("changed input: " + record["path"])
    return Path(record["path"])


def write(path, value):
    with Path(path).open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def table(path, rows):
    with Path(path).open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def close(a, b, tolerance):
    if not (math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tolerance):
        raise ValueError(f"saved algebra differs: {a!r} versus {b!r}")


def coordinates(record):
    lines = verified(record).read_text().splitlines()
    rows = [line.split() for line in lines[2:] if line.strip()]
    if len(rows) != int(lines[0]) or any(len(row) != 4 for row in rows):
        raise ValueError("malformed pinned XYZ")
    return [(row[0], tuple(float(v) for v in row[1:])) for row in rows]


def paired(matrix, tolerance):
    errors = {}
    for candidate in ORDER:
        ca, la = [coordinates(matrix[metal][candidate]["xyz"]) for metal in METALS]
        if len(ca) != len(la):
            raise ValueError("paired atom inventory differs")
        swaps = [(a[0], b[0]) for a, b in zip(ca, la) if a[0] != b[0]]
        if swaps != [("Ca", "La")]:
            raise ValueError("paired metal substitution differs")
        error = max(abs(x - y) for a, b in zip(ca, la) for x, y in zip(a[1], b[1]))
        if error > tolerance:
            raise ValueError("paired coordinates differ")
        errors[candidate] = error
    return errors


def analyze(config):
    paths = {name: verified(record) for name, record in config["inputs"].items()}
    recovered = read(paths["recovered_collection"])
    native = read(paths["recovered_native"])
    rm = read(paths["recovered_manifest"])
    controls = read(paths["control_collection"])
    cm = read(paths["control_manifest"])
    canonical = read(paths["canonical_collection"])
    findings = read(paths["prior_findings"])
    if rm["model"] != cm["model"] or rm["settings"] != cm["settings"] or rm["runtime"] != cm["runtime"]:
        raise ValueError("controls have a different electronic/solvent protocol")
    if recovered["reference"]["source"] != config["inputs"]["canonical_collection"]:
        raise ValueError("frozen canonical reference changed")
    if controls["calibration"] != canonical["calibration"]:
        raise ValueError("transfer recalibrated independently")
    settings = rm["settings"]
    ev, ha = settings["EV_TO_KCAL"], settings["HA_TO_KCAL"]
    tolerance = config["saved_algebra_tolerance_kcal_mol"]
    actual = {r["source_case_id"]: r for r in recovered["rows"]}
    native_rows = {r["source_case_id"]: r for r in native["cases"]}
    historic = {r["source_case_id"]: r for r in controls["rows"]}
    historic_cases = {r["source_case_id"]: r for r in cm["cases"]}
    flags = {r["source_case_id"]: r for r in findings["rows"]}
    rows, cell_rows, work_rows = [], [], []
    common_signature = None
    for source in config["sources"]:
        sid, cohort = source["source_case_id"], source["cohort"]
        is_recovered = cohort == "recovered_two_Ca"
        score = (actual if is_recovered else historic)[sid]
        original = (native_rows if is_recovered else historic_cases)[sid]
        matrix = original["matrix"]
        if score["biological_group"] != "q88jh5" or score["expected_class"] != "Ca":
            raise ValueError("unexpected reference identity")
        if score["pool"]["status"] != "available":
            raise ValueError("required saved pool unavailable: " + sid)
        signature = read(verified(original["state_signature"]))
        if common_signature is None:
            common_signature = signature
        if signature != common_signature:
            raise ValueError("chemical state differs across comparisons: " + sid)
        pair_errors = paired(matrix, config["paired_coordinate_tolerance_A"])
        cells = {metal: {} for metal in METALS}
        for metal in METALS:
            for candidate in ORDER:
                cell = matrix[metal][candidate]
                component = score["components"][metal + "__" + candidate]
                if cell["charge"] != {"Ca": -3, "La": -2}[metal] or cell["multiplicity"] != 1:
                    raise ValueError("paired state differs")
                if not component["accepted"]:
                    raise ValueError("required component unaccepted")
                mace = read(verified(cell["MACE"]))
                raw_mace = cell["energy_eV"] if is_recovered else cell["MACE_eV"]
                if mace["energy_eV"] != raw_mace:
                    raise ValueError("archived OMOL value differs")
                solvent = component["primary_water_energy_hartree"] - component["primary_vacuum_energy_hartree"]
                close(solvent, component["primary_solvation_hartree"], 1e-12)
                value = {"OMOL_eV": raw_mace,
                         "GFN2_vacuum_Eh": component["primary_vacuum_energy_hartree"],
                         "GFN2_ALPB_Eh": component["primary_water_energy_hartree"],
                         "solvent_Eh": solvent, "OMOL_kcal_mol": raw_mace * ev,
                         "solvent_kcal_mol": solvent * ha,
                         "component_residual_Eh": component["component_error_hartree"],
                         "population_residual_e": component["maximum_population_error_e"]}
                value["composite_kcal_mol"] = value["OMOL_kcal_mol"] + value["solvent_kcal_mol"]
                cells[metal][candidate] = value
                cell_rows.append({"source_case_id": sid, "cohort": cohort, "metal": metal,
                                  "candidate": candidate, **value})
        works = {metal: {} for metal in METALS}
        for metal in METALS:
            origin = cells[metal]["origin"]
            for candidate in ORDER:
                cell = cells[metal][candidate]
                v = (cell["OMOL_eV"] - origin["OMOL_eV"]) * ev
                u = (cell["solvent_Eh"] - origin["solvent_Eh"]) * ha
                works[metal][candidate] = {"OMOL_kcal_mol": v, "solvent_kcal_mol": u,
                                          "composite_kcal_mol": v + u}
                work_rows.append({"source_case_id": sid, "cohort": cohort, "metal": metal,
                                  "candidate": candidate, "is_own_proposal": candidate == "adaptive_" + metal,
                                  "selected_operational": score["pool"]["operational"][metal]["candidate"] == candidate,
                                  **works[metal][candidate]})
        contrasts = {}
        for policy in ("origin", "operational", "mathematical"):
            old = score["pool"][policy]
            qa, qb = [old[metal]["candidate"] for metal in METALS]
            a, b = cells["Ca"][qa], cells["La"][qb]
            v, u = (a["OMOL_eV"] - b["OMOL_eV"]) * ev, (a["solvent_Eh"] - b["solvent_Eh"]) * ha
            close(v + u, old["R_model_kcal_mol"], tolerance)
            for metal, q in zip(METALS, (qa, qb)):
                close(works[metal][q]["composite_kcal_mol"], old[metal]["work_kcal_mol"], tolerance)
            contrasts[policy] = {"Ca_candidate": qa, "La_candidate": qb, "OMOL_kcal_mol": v,
                                 "solvent_kcal_mol": u, "composite_kcal_mol": v + u,
                                 "saved_decision": score["decisions"][policy]}
        minima = {}
        for metal in METALS:
            vacuum_q = min(ORDER, key=lambda q: works[metal][q]["OMOL_kcal_mol"])
            composite_q = min(ORDER, key=lambda q: works[metal][q]["composite_kcal_mol"])
            saved_q = score["pool"]["mathematical"][metal]["candidate"]
            close(works[metal][composite_q]["composite_kcal_mol"], works[metal][saved_q]["composite_kcal_mol"], tolerance)
            expected_operational = "origin" if works[metal][saved_q]["composite_kcal_mol"] >= -settings["origin_selection_tolerance_kcal_mol"] else saved_q
            if score["pool"]["operational"][metal]["candidate"] != expected_operational:
                raise ValueError("operational selection differs from frozen policy")
            minima[metal] = {"finite_pool_vacuum_minimum": vacuum_q, "finite_pool_composite_minimum": composite_q,
                             "saved_operational": score["pool"]["operational"][metal]["candidate"]}
        selected = score["pool"]["operational"]
        wa, wb = [works[metal][selected[metal]["candidate"]] for metal in METALS]
        change = {key: wa[key] - wb[key] for key in ("OMOL_kcal_mol", "solvent_kcal_mol", "composite_kcal_mol")}
        close(change["composite_kcal_mol"], contrasts["operational"]["composite_kcal_mol"] - contrasts["origin"]["composite_kcal_mol"], tolerance)
        rows.append({"source_case_id": sid, "cohort": cohort, "role": source["role"],
                     "state_signature": original["state_signature"], "paired_coordinate_max_error_A": pair_errors,
                     "contrasts": contrasts, "cells": cells, "works": works, "selection": minima,
                     "accommodation_contrast_change": change,
                     "component_input_pins": {k: v["inputs"] for k, v in score["components"].items()},
                     "proposal_diagnostics": flags[sid]["proposals"] if is_recovered else None})
    anchor = next(r for r in rows if r["source_case_id"] == config["canonical_anchor_source_id"])
    for row in rows:
        row["origin_minus_canonical_anchor"] = {
            key: row["contrasts"]["origin"][key] - anchor["contrasts"]["origin"][key]
            for key in ("OMOL_kcal_mol", "solvent_kcal_mol", "composite_kcal_mol")}
    spans = {}
    for cohort in ("recovered_two_Ca", "historical_La"):
        subset = [r for r in rows if r["cohort"] == cohort]
        lo = min(subset, key=lambda r: r["contrasts"]["origin"]["composite_kcal_mol"])
        hi = max(subset, key=lambda r: r["contrasts"]["origin"]["composite_kcal_mol"])
        spans[cohort] = {"minimum_source": lo["source_case_id"], "maximum_source": hi["source_case_id"],
                         "same_extrema_component_difference": {key: hi["contrasts"]["origin"][key] - lo["contrasts"]["origin"][key]
                            for key in ("OMOL_kcal_mol", "solvent_kcal_mol", "composite_kcal_mol")}}
    return {"status": "complete", "sources": config["inputs"], "rows": rows, "origin_spans": spans,
            "canonical_anchor_source_id": config["canonical_anchor_source_id"], "reference": recovered["reference"],
            "method_compatibility": {"model_equal": True, "settings_equal": True, "runtime_equal": True, "state_signature_equal": True},
            "settings": settings, "runtime": rm["runtime"],
            "numerical_residuals": {"maximum_abs_component_Eh": max(abs(r["component_residual_Eh"]) for r in cell_rows),
                                    "maximum_population_e": max(r["population_residual_e"] for r in cell_rows),
                                    "component_gate_Eh": settings["component_parity_Eh"],
                                    "population_gate_e": settings["population_parity_e"],
                                    "saved_algebra_tolerance_kcal_mol": tolerance},
            "controls_unavailable": config["controls_unavailable"], "new_molecular_calls": 0, "new_fit": False,
            "classification_policy_changed": False, "production_changed": False,
            "interpretation": "Components locate the source of the model contrast, not proof of bad folds or a unique physical mechanism. Policy margin analysis is separate."}, cell_rows, work_rows


def render(result):
    lines = ["# PedE: saved OMOL and solvent contributions", "",
             "All five recovered two-Ca sources and four exact-method La-conditioned controls are retained.",
             "No molecular calculation, fit, new vote, or source selection was performed.", "",
             "## Origin contrasts relative to the frozen canonical PedE source", "",
             "All entries below are kcal/mol. Components add to the raw contrast; these are not affinities.",
             "Origin and accommodated decisions use different frozen bands. Classification-margin attribution belongs to the separate parent analysis.", "",
             "| Source | OMOL difference | Solvent difference | Composite difference | Origin / accommodated call |",
             "|---|---:|---:|---:|---|"]
    for row in result["rows"]:
        d = row["origin_minus_canonical_anchor"]
        calls = [row["contrasts"][p]["saved_decision"]["decision"] for p in ("origin", "operational")]
        label = row["cohort"] + "/sample-" + row["source_case_id"].rsplit("sample-", 1)[-1]
        lines.append(f"| {label} | {d['OMOL_kcal_mol']:.6f} | {d['solvent_kcal_mol']:.6f} | {d['composite_kcal_mol']:.6f} | {' / '.join(calls)} |")
    lines += ["", "## Origin spread using the same composite extrema", "",
              "Component differences use the sources defining the composite minimum/maximum; independent component ranges are not added.", ""]
    for cohort, span in result["origin_spans"].items():
        d = span["same_extrema_component_difference"]
        lines.append(f"- {cohort}: OMOL {d['OMOL_kcal_mol']:.6f}, solvent {d['solvent_kcal_mol']:.6f}, total {d['composite_kcal_mol']:.6f} kcal/mol.")
        lines.append(f"  Sources: `{span['minimum_source']}` → `{span['maximum_source']}`.")
    lines += ["", "## Each selected accommodation work", "",
              "Negative work lowers an endpoint energy. Ca-minus-La work changes R. The vacuum search proposes candidates; composite energy selects them.", "",
              "| Source | Metal / chosen geometry | OMOL work | Solvent work | Composite work | Vacuum pool minimum |",
              "|---|---|---:|---:|---:|---|"]
    for row in result["rows"]:
        for metal in METALS:
            q = row["selection"][metal]["saved_operational"]
            w = row["works"][metal][q]
            label = row["cohort"] + "/sample-" + row["source_case_id"].rsplit("sample-", 1)[-1]
            lines.append(f"| {label} | {metal} / {q} | {w['OMOL_kcal_mol']:.6f} | {w['solvent_kcal_mol']:.6f} | {w['composite_kcal_mol']:.6f} | {row['selection'][metal]['finite_pool_vacuum_minimum']} |")
    numerical = result["numerical_residuals"]
    lines += ["", "## Numerical checks", "",
              f"Maximum absolute saved component residual: {numerical['maximum_abs_component_Eh']:.12g} Eh (frozen gate {numerical['component_gate_Eh']:.12g} Eh).",
              f"Maximum saved population residual: {numerical['maximum_population_e']:.12g} e (frozen gate {numerical['population_gate_e']:.12g} e).",
              "All required component records were accepted; paired coordinates, local state, model and runtime/settings compatibility passed.",
              "These checks do not establish biological correctness or independent electronic-method accuracy.",
              "", "## Limits and artifacts", "",
              "`CELLS.tsv` retains all 54 metal/geometry cells, including both cross evaluations; `WORKS.tsv` retains every origin-relative work.",
              "`RESULT.json` retains component residuals, paired/state checks, exact input pins and saved proposal/boundary diagnostics.",
              "The original absolute-v1 backend parity failure remains failed. Component-v2 acceptance is separate.",
              "Historical La sample 1 is absent from this exact corrected comparison manifest and remains explicitly unavailable here.",
              "These are correlated consumed source hypotheses. Different folds can change several geometric features together; an energy decomposition cannot establish whether a fold or Hamiltonian is wrong.",
              "No source was optimized or rescored, no calibration changed, and no new five-fold aggregation is defined."]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    job = os.environ.get("SLURM_JOB_ID", "")
    if not job.isdigit():
        raise SystemExit("Actual Slurm allocation required before saved-artifact scientific analysis")
    control = subprocess.run(["scontrol", "show", "job", job, "-o"], text=True, capture_output=True, check=True)
    if f"JobId={job} " not in control.stdout:
        raise SystemExit("Slurm allocation could not be verified")
    config = read(args.config)
    if pin(__file__) != config["driver"]:
        raise SystemExit("driver differs from frozen config")
    out = Path(config["output"])
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    receipt = {"config": pin(args.config), "driver": pin(__file__), "job_id": job,
               "allocation": control.stdout, "new_molecular_calls": 0, "status": "running"}
    try:
        result, cells, works = analyze(config)
        write(out / "RESULT.json", result)
        table(out / "CELLS.tsv", cells)
        table(out / "WORKS.tsv", works)
        (out / "REPORT.md").write_text(render(result))
        receipt.update(status="complete", outputs={n: pin(out / n) for n in ("RESULT.json", "CELLS.tsv", "WORKS.tsv", "REPORT.md")})
    except BaseException as exc:
        receipt.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        receipt.update(elapsed_seconds=time.monotonic() - start, process_CPU_seconds=usage.ru_utime + usage.ru_stime,
                       peak_RSS_KiB=usage.ru_maxrss)
        write(out / f"EXECUTION_{job}.json", receipt)


if __name__ == "__main__":
    main()
