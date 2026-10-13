#!/usr/bin/env python3
"""Exact saved xTB solvent decomposition; no electronic evaluation or new model."""
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

TERMS = ("Gelec", "Gsasa", "Ghb", "Gshift", "intrinsic_remainder")


def read(p):
    return json.loads(Path(p).read_text())


def pin(p):
    p = Path(p).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    if pin(record["path"]) != record:
        raise ValueError("changed pinned input: " + record["path"])
    return Path(record["path"])


def dump(p, value):
    with Path(p).open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def table(p, rows):
    with Path(p).open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def scalar(text, label, required=True):
    values = re.findall(r"^\s*::\s*(?:->\s*)?" + re.escape(label) + r"\s+([-+0-9.Ee]+)\s+Eh\s*::", text, re.M)
    if len(values) != 1:
        if not required and not values:
            return None
        raise ValueError("expected one summary scalar: " + label)
    value = float(values[0])
    if not math.isfinite(value):
        raise ValueError("nonfinite summary scalar")
    return value


def stdout_record(row, tolerance):
    point = read(verify(row["result"]))
    record = row["backends"]["xtb"]
    text = verify(record["stdout"]).read_text()
    receipt = read(verify(record["execution"]))
    if receipt["returncode"] != 0 or point["backends"]["xtb"]["status"] != "available":
        raise ValueError("unavailable saved endpoint")
    if receipt["stdout"] != record["stdout"]:
        raise ValueError("execution/stdout mismatch")
    total = scalar(text, "total energy")
    if abs(total - point["backends"]["xtb"]["energy_hartree"]) > tolerance:
        raise ValueError("summary/accepted total mismatch")
    components = {key: scalar(text, key, key in ("SCC energy", "repulsion energy", "add. restraining"))
                  for key in ("SCC energy", "repulsion energy", "add. restraining", "isotropic ES", "anisotropic ES", "anisotropic XC", "dispersion", "Gsolv", "Gelec", "Gsasa", "Ghb", "Gshift")}
    closure = total - sum(components[k] for k in ("SCC energy", "repulsion energy", "add. restraining"))
    if abs(closure) > tolerance:
        raise ValueError("printed total/SCC/repulsion closure failed")
    setup_keys = ("Parameter file", "Dielectric constant", "Reference state", "Free energy shift", "Temperature", "Interaction kernel", "Born radius scaling", "Born radii integrator", "Born offset", "H-bond correction", "Ion screening", "Surface tension", "Grid points", "electronic temp.", "# electrons", "net charge", "unpaired electrons")
    setup = [line for line in text.splitlines() if any(k in line for k in setup_keys)]
    thermo = [line for line in text.splitlines() if re.search(r"entropy|fermi|electronic temp|free energy|Gshift", line, re.I)]
    warnings = [line for line in text.splitlines() if "WARNING" in line or "Unknown reference state" in line]
    return {"point": row["result"], "stdout": record["stdout"], "task": point["task"], "total_Eh": total,
            "components_Eh": components, "printed_total_closure_Eh": closure, "setup_lines": setup,
            "temperature_free_energy_entropy_lines": thermo, "separate_electronic_entropy_Eh": None,
            "separate_electronic_entropy_status": "not separately printed; included in SCC according to upstream source; not assumed zero",
            "warnings": warnings}


def analyze(config):
    inputs = {k: verify(v) for k, v in config["inputs"].items()}
    inventory, energy = read(inputs["inventory"]), read(inputs["energetics"])
    if inventory["source"] != config["inputs"]["energetics"]:
        raise ValueError("mismatched energy/inventory")
    for record in config["source_pins"]:
        verify(record)
    records = {(r["source_case_id"], r["metal"], r["medium"]): r for r in inventory["rows"]}
    expected = {(r["source_case_id"], metal, medium) for r in energy["rows"] for metal in ("Ca", "La") for medium in ("vacuum", "water")}
    if set(records) != expected or len(records) != 36:
        raise ValueError("declared origin population changed")
    manifests = [read(inputs[k]) for k in ("recovered_manifest", "historical_manifest")]
    source_cases = {r["source_case_id"]: r for m in manifests for r in m["cases"]}
    conversion = energy["settings"]["HA_TO_KCAL"]
    rows, endpoints, flat = [], [], []
    for original in energy["rows"]:
        sid = original["source_case_id"]
        case = source_cases[sid]
        prep_pin = case["source_preparation"]
        prep = read(verify(prep_pin))
        state = read(verify(case["state_signature"]))
        charges = prep.get("source_formal_charge_ledger")
        cell = {"source_case_id": sid, "cohort": original["cohort"], "metals": {}, "source_preparation": prep_pin,
                "chemical_state": state, "source_formal_charge_ledger": charges,
                "caps": prep["mapping"]["cut_bonds_and_caps"], "protons_or_waters_added": prep.get("new_protons_or_waters"),
                "source_charge_ledger_status": "present" if charges is not None else "not present in archived preparation; state signature retained"}
        for metal in ("Ca", "La"):
            vacuum = stdout_record(records[(sid, metal, "vacuum")], config["closure_tolerance_Eh"])
            water = stdout_record(records[(sid, metal, "water")], config["closure_tolerance_Eh"])
            for key in ("xyz", "charge", "multiplicity", "uhf", "valence_electrons", "atom_count"):
                if vacuum["task"][key] != water["task"][key]:
                    raise ValueError("medium changed atom/state inventory")
            if water["task"]["charge"] != {"Ca": -3, "La": -2}[metal] or water["task"]["valence_electrons"] != 588:
                raise ValueError("unexpected PedE electron inventory")
            c = water["components_Eh"]
            if any(c[k] is None for k in ("Gsolv", "Gelec", "Gsasa", "Ghb", "Gshift")):
                raise ValueError("required actual ALPB decomposition unavailable")
            u = water["total_Eh"] - vacuum["total_Eh"]
            if abs(u - original["cells"][metal]["origin"]["solvent_Eh"]) > config["closure_tolerance_Eh"]:
                raise ValueError("saved solvent value changed")
            gclosure = c["Gsolv"] - sum(c[k] for k in TERMS[:-1])
            if abs(gclosure) > config["closure_tolerance_Eh"]:
                raise ValueError("Gsolv component closure failed")
            terms = {k: c[k] for k in TERMS[:-1]}
            terms["intrinsic_remainder"] = u - c["Gsolv"]
            residual = u - sum(terms.values())
            details = {"U_Eh": u, "Gsolv_Eh": c["Gsolv"], "terms_Eh": terms,
                       "terms_kcal_mol": {k: v * conversion for k, v in terms.items()},
                       "Gsolv_closure_Eh": gclosure, "U_decomposition_closure_Eh": residual,
                       "printed_component_medium_changes_Eh": {k: (c[k] - vacuum["components_Eh"][k]) if c[k] is not None and vacuum["components_Eh"][k] is not None else None
                           for k in ("SCC energy", "repulsion energy", "add. restraining", "isotropic ES", "anisotropic ES", "anisotropic XC", "dispersion")},
                       "vacuum": vacuum, "water": water}
            cell["metals"][metal] = details
            flat.append({"source_case_id": sid, "cohort": original["cohort"], "metal": metal,
                         "U_kcal_mol": u * conversion, "Gsolv_kcal_mol": c["Gsolv"] * conversion,
                         **details["terms_kcal_mol"]})
            endpoints.extend([vacuum, water])
        cell["contrast_terms_kcal_mol"] = {k: cell["metals"]["Ca"]["terms_kcal_mol"][k] - cell["metals"]["La"]["terms_kcal_mol"][k] for k in TERMS}
        rows.append(cell)
    extrema = energy["origin_spans"]["recovered_two_Ca"]
    by_id = {r["source_case_id"]: r for r in rows}
    lo, hi = [by_id[extrema[k]] for k in ("minimum_source", "maximum_source")]
    delta = {"from_source": lo["source_case_id"], "to_source": hi["source_case_id"],
             "Ca_minus_La_terms_kcal_mol": {k: hi["contrast_terms_kcal_mol"][k] - lo["contrast_terms_kcal_mol"][k] for k in TERMS},
             "by_metal_kcal_mol": {m: {k: hi["metals"][m]["terms_kcal_mol"][k] - lo["metals"][m]["terms_kcal_mol"][k] for k in TERMS} for m in ("Ca", "La")}}
    expected_delta = extrema["same_extrema_component_difference"]["solvent_kcal_mol"]
    delta["saved_solvent_double_difference_kcal_mol"] = expected_delta
    delta["sum_residual_kcal_mol"] = sum(delta["Ca_minus_La_terms_kcal_mol"].values()) - expected_delta
    if abs(delta["sum_residual_kcal_mol"]) > config["closure_tolerance_Eh"] * conversion:
        raise ValueError("extrema decomposition does not reproduce saved contrast")
    return {"status": "complete", "inputs": config["inputs"], "source_pins": config["source_pins"],
            "rows": rows, "sample3_to_sample1": delta, "new_molecular_calls": 0,
            "new_population_analysis": False, "new_correction": False,
            "interpretation": "Exact identity on saved self-consistent endpoints. Intrinsic remainder includes all terms outside printed Gsolv, including any electronic entropy and density-dependent dispersion. It is not a separately validated polarization/relaxation model."}, flat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    job = os.environ.get("SLURM_JOB_ID", "")
    if not job.isdigit():
        raise SystemExit("Actual Slurm allocation required")
    allocation = subprocess.check_output(["scontrol", "show", "job", job, "-o"], text=True)
    if f"JobId={job} " not in allocation:
        raise SystemExit("unverified Slurm allocation")
    config = read(args.config)
    if pin(__file__) != config["driver"]:
        raise SystemExit("changed driver")
    output = Path(config["output"])
    output.mkdir(parents=True, exist_ok=False)
    receipt = {"job_id": job, "config": pin(args.config), "allocation": allocation, "new_molecular_calls": 0}
    start = time.monotonic()
    try:
        result, flat = analyze(config)
        dump(output / "RESULT.json", result)
        table(output / "COMPONENTS.tsv", flat)
        receipt.update(status="complete", outputs={k: pin(output / k) for k in ("RESULT.json", "COMPONENTS.tsv")})
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
