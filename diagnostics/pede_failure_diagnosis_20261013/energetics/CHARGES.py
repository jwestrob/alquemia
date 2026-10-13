#!/usr/bin/env python3
"""Read archived origin populations and solvent contributions; no calculators."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import subprocess
import time


def read(path):
    return json.loads(Path(path).read_text())


def pin(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    if pin(record["path"]) != record:
        raise ValueError("changed artifact: " + record["path"])
    return Path(record["path"])


def dump(path, value):
    with Path(path).open("x") as out:
        json.dump(value, out, indent=2, allow_nan=False)
        out.write("\n")


def tsv(path, rows):
    with Path(path).open("x", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def groups(preparation, xyz, tolerance):
    """Full preparation QM indices include fixed PQQ H absent from motion maps."""
    atoms = [line.split() for line in xyz.splitlines()[2:] if line.strip()]
    if len(atoms) != int(xyz.splitlines()[0]) or atoms[0][0] not in ("Ca", "La"):
        raise ValueError("unexpected origin XYZ")
    assignment = {0: "metal"}
    max_error = 0.0
    for row in preparation["mapping"]["source_to_qm"]:
        i = row["qm_index"]
        if i in assignment:
            raise ValueError("duplicate mapped QM index")
        source = row.get("source")
        if row["kind"] == "sigma_link_H":
            group = "synthetic_caps"
        elif row["kind"] == "opaque_cofactor" and ":PQQ" in row["fragment"]:
            group = "PQQ_including_fixed_H"
        elif source:
            group = "residue:{chain}:{resname}:{resnum}{insertion_code}".format(**source)
        else:
            raise ValueError("unassigned actual atom: " + str(row))
        expected_element = row.get("element", source.get("element") if source else "H")
        if atoms[i][0] != expected_element:
            raise ValueError("map/XYZ element mismatch")
        error = max(abs(float(a) - b) for a, b in zip(atoms[i][1:], row["xyz_A"]))
        max_error = max(max_error, error)
        if error > tolerance:
            raise ValueError("preparation/XYZ coordinate mismatch")
        assignment[i] = group
    if set(assignment) != set(range(len(atoms))):
        raise ValueError("preparation does not map every QM atom")
    grouped = {}
    for i, group in assignment.items():
        grouped.setdefault(group, []).append(i)
    return grouped, max_error, len(atoms)


def iterations(text, backend):
    # Require the actual SCC row shapes; do not count orbital-population tables.
    rows = []
    for line in text.splitlines():
        fields = line.split()
        valid = (len(fields) == 7 and fields[-1] in ("T", "F")) if backend == "xtb" else len(fields) == 4
        if not valid or not fields[0].isdigit():
            continue
        try:
            values = [float(v) for v in fields[1:-1] if backend == "xtb"] if backend == "xtb" else [float(v) for v in fields[1:]]
        except ValueError:
            continue
        if len(values) != (5 if backend == "xtb" else 3):
            continue
        rows.append({"iteration": int(fields[0]), "line": line, "energy_printed_Eh": values[0],
                     "delta_energy_printed_Eh": values[1], "printed_residual": values[2]})
    if not rows or [r["iteration"] for r in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("could not identify one complete SCC table: " + backend)
    return {"count": len(rows), "last": rows[-1], "raw_rows": rows}


def analyze(config):
    inputs = {k: verify(v) for k, v in config["inputs"].items()}
    energetic = read(inputs["energetics"])
    inventory = read(inputs["inventory"])
    if inventory["source"] != config["inputs"]["energetics"]:
        raise ValueError("inventory belongs to a different energy analysis")
    manifests = [read(inputs[k]) for k in ("recovered_manifest", "historical_manifest")]
    cases = {c["source_case_id"]: c for m in manifests for c in m["cases"]}
    expected = {(r["source_case_id"], metal, medium) for r in energetic["rows"]
                for metal in ("Ca", "La") for medium in ("vacuum", "water")}
    actual = {(r["source_case_id"], r["metal"], r["medium"]) for r in inventory["rows"]}
    if expected != actual or len(expected) != len(inventory["rows"]):
        raise ValueError("origin inventory differs from declared nine sources")
    out, flat, atoms = [], [], []
    for row in inventory["rows"]:
        sid, metal, medium = [row[k] for k in ("source_case_id", "metal", "medium")]
        point = read(verify(row["result"]))
        xyz = verify(row["xyz"]).read_text()
        prep_pin = cases[sid]["source_preparation"]
        grouped, coordinate_error, atom_count = groups(read(verify(prep_pin)), xyz, config["coordinate_tolerance_A"])
        if point["task"]["xyz"] != row["xyz"] or point["task"]["charge"] != {"Ca": -3, "La": -2}[metal] or point["task"]["multiplicity"] != 1:
            raise ValueError("origin coordinate/state mismatch")
        result = {"source_case_id": sid, "cohort": row["cohort"], "metal": metal, "medium": medium,
                  "point": row["result"], "preparation": prep_pin, "groups_QM_indices": grouped,
                  "mapping_coordinate_error_A": coordinate_error, "backends": {}}
        for backend, records in row["backends"].items():
            saved = point["backends"][backend]
            execution = read(verify(records["execution"]))
            stdout = verify(records["stdout"]).read_text()
            population = verify(records["population_file"])
            scientific = read(population) if backend == "tblite" else None
            q = scientific["charges_e"] if scientific else [float(v) for v in population.read_text().split()]
            if q != saved["charges_e"] or len(q) != atom_count or not all(math.isfinite(v) for v in q):
                raise ValueError("raw population differs from admitted record")
            if execution["returncode"] != 0 or saved["status"] != "available" or execution["stdout"] != records["stdout"]:
                raise ValueError("unavailable execution")
            declared = saved["solver_evidence"].get("convergence_declared") if backend == "xtb" else saved["solver_evidence"].get("library_returned_without_error")
            if declared is not True:
                raise ValueError("saved convergence unavailable")
            grouped_q = {name: sum(q[i] for i in indices) for name, indices in grouped.items()}
            solver = iterations(stdout, backend)
            warnings = [line for line in stdout.splitlines() if re.search(r"warn|HOMO.LUMO GAP|SCC converged|convergence criteria satisfied", line, re.I)]
            occupations = scientific.get("orbital_occupations_e") if scientific else None
            partial = ([{"orbital_0based": i, "occupation_e": n} for i, n in enumerate(occupations) if 0.001 < n < 1.999]
                       if occupations is not None else None)
            result["backends"][backend] = {"charge_sum_e": sum(q), "charge_sum_residual_e": sum(q) - point["task"]["charge"],
                "groups_e": grouped_q, "convergence_declared": declared, "iterations": solver,
                "raw_warning_gap_convergence_lines": warnings, "fractional_occupations_0p001_to_1p999": partial,
                "partial_occupation_interpretation": "Descriptive finite-300K-electronic-temperature occupations; no state/branch assignment.",
                "artifact_pins": records}
            flat.append({"source_case_id": sid, "cohort": row["cohort"], "metal": metal, "medium": medium,
                         "backend": backend, "iterations": solver["count"], "last_delta_E_Eh": solver["last"]["delta_energy_printed_Eh"],
                         "last_printed_residual": solver["last"]["printed_residual"], "charge_sum_e": sum(q), **grouped_q})
            for group, indices in grouped.items():
                for i in indices:
                    atoms.append({"source_case_id": sid, "metal": metal, "medium": medium, "backend": backend,
                                  "qm_index_0based": i, "group": group, "charge_e": q[i]})
        out.append(result)
    solvent = []
    for row in energetic["rows"]:
        solvent.append({"source_case_id": row["source_case_id"], "cohort": row["cohort"],
                        **{metal + "_solvent_kcal_mol": row["cells"][metal]["origin"]["solvent_kcal_mol"] for metal in ("Ca", "La")}})
    recovered = {r["source_case_id"].rsplit("sample-", 1)[-1]: r for r in solvent if r["cohort"] == "recovered_two_Ca"}
    differences = {metal: recovered["1"][metal + "_solvent_kcal_mol"] - recovered["3"][metal + "_solvent_kcal_mol"] for metal in ("Ca", "La")}
    differences["Ca_minus_La"] = differences["Ca"] - differences["La"]
    return {"status": "complete", "inputs": config["inputs"], "origin_solvent": solvent,
            "sample1_minus_sample3_solvent_kcal_mol": differences, "rows": out,
            "new_molecular_calls": 0, "calibration_changed": False,
            "interpretation": "Population differences are diagnostics; convergence and cross-backend agreement do not establish the correct electronic branch or biological accuracy."}, flat, atoms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    job = os.environ.get("SLURM_JOB_ID", "")
    if not job.isdigit():
        raise SystemExit("Actual Slurm allocation required")
    allocation = subprocess.check_output(["scontrol", "show", "job", job, "-o"], text=True)
    if f"JobId={job} " not in allocation:
        raise SystemExit("Slurm allocation unverified")
    config = read(args.config)
    verify(config["driver"])
    if pin(__file__) != config["driver"]:
        raise SystemExit("wrong driver")
    output = Path(config["output"])
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    receipt = {"job_id": job, "config": pin(args.config), "allocation": allocation, "new_molecular_calls": 0}
    try:
        result, flat, atoms = analyze(config)
        dump(output / "RESULT.json", result)
        tsv(output / "ORIGIN_CHARGES.tsv", flat)
        tsv(output / "ATOM_CHARGES.tsv", atoms)
        tsv(output / "ORIGIN_SOLVENT.tsv", result["origin_solvent"])
        receipt.update(status="complete", outputs={n: pin(output / n) for n in ("RESULT.json", "ORIGIN_CHARGES.tsv", "ATOM_CHARGES.tsv", "ORIGIN_SOLVENT.tsv")})
    except BaseException as error:
        receipt.update(status="failed", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        receipt.update(elapsed_seconds=time.monotonic() - start, process_CPU_seconds=usage.ru_utime + usage.ru_stime,
                       peak_RSS_KiB=usage.ru_maxrss)
        dump(output / f"EXECUTION_{job}.json", receipt)


if __name__ == "__main__":
    main()
