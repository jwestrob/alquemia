#!/usr/bin/env python3
"""Independent read-only recount of the finite October 12 native controls.

No chemistry tools, runtime collector imports or new scientific outputs.
Raw stdout/property tables are parsed independently; solver comparisons use
the pinned archived macOS xTB/tblite records, not new upstream runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

NUM = r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def section(text, name):
    return re.findall(r"^\$" + re.escape(name) + r"\s*\n(.*?)^\$End", text, re.M | re.S)[-1]


def scalar(text, name):
    return float(re.findall(r"&" + name + r"\s+\[[^\]]+\]\s+(" + NUM + r")", text)[-1])


def parse_raw(directory):
    out = (directory / "endpoint.out").read_text()
    prop = (directory / "endpoint.property.txt").read_text()
    energy = scalar(section(prop, "Single_Point_Data"), "FinalEnergy")
    stdout_energy = float(re.findall(r"^FINAL SINGLE POINT ENERGY\s+(" + NUM + r")", out, re.M)[-1])
    population = section(prop, "SCF_Mulliken_Population_Analysis")
    atom_count = int(scalar(population, "NAtoms"))
    array = population.split("&AtomicCharges", 1)[1].split("&Method", 1)[0]
    charge_rows = re.findall(r"^\s*(\d+)\s+(" + NUM + r")\s*$", array, re.M)
    assert [int(i) for i, _ in charge_rows] == list(range(atom_count))
    charges = [float(q) for _, q in charge_rows]
    info = section(prop, "Calculation_Info")
    state = {key: int(scalar(info, field)) for key, field in
             (("charge", "Charge"), ("multiplicity", "Mult"), ("atoms", "NumOfAtoms"),
              ("valence_electrons", "NumOfElectrons"), ("basis_functions", "NumOfBasisFuncts"))}
    gates = re.findall(r"Last\s+(Energy|MAX-Density|RMS-Density)\s+change\s+\.\.\.\s+(" + NUM +
                       r")\s+Tolerance\s*:\s*(" + NUM + r")", out)
    assert len(gates) == 3
    failures = [{"criterion": name, "residual": float(value), "printed_tolerance": float(tol)}
                for name, value, tol in gates if abs(float(value)) > float(tol)]
    return {"energy_hartree": energy, "charges_e": charges, "state": state,
            "printed_gate_failures": failures, "printed_gates": gates,
            "stdout_energy_difference_hartree": abs(stdout_energy-energy),
            "normal": "ORCA TERMINATED NORMALLY" in out and '"NORMAL TERMINATION"' in prop,
            "converged": "SCF CONVERGED AFTER" in out and
                         bool(re.search(r'&Converged\s+\[[^\]]+\]\s+true', section(prop, "Single_Point_Data"))),
            "runtime_read_files": re.findall(r"^\s*Reading XTB parameter file\s+([^\r\n]+?)\s*$", out, re.M)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workspace", type=Path, required=True)
    ap.add_argument("--reported", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    assert not args.output.exists(), "Use a fresh report path"
    manifests = [read(args.workspace / f"manifest_{name}.json") for name in ("La8", "aquo50", "Ca2")]
    expected = {j["id"]: j for m in manifests for j in m["jobs"]}
    assert len(expected) == 60
    upstream_pin = manifests[1]["upstream_data"]
    upstream_path = Path(upstream_pin["path"])
    assert sha(upstream_path) == upstream_pin["sha256"]
    upstream = {r["id"]: r for r in read(upstream_path)["runs"]}
    files_checked = 0
    rows, raw = {}, {}
    for job_id, job in sorted(expected.items()):
        d = args.workspace / "results" / job_id
        receipt = read(d / "receipt.json")
        assert receipt["job"] == job
        for name, p in {**job["files"], **receipt["artifacts"]}.items():
            assert Path(name).name == name
            assert sha(d / name) == p["sha256"], f"changed artifact: {job_id}/{name}"
            files_checked += 1
        for name, p in job["files"].items():
            assert sha(d / name) == p["sha256"]
        x = parse_raw(d)
        assert x["charges_e"] == receipt["output"]["charges_e"]
        assert abs(x["energy_hartree"] - receipt["output"]["energy_hartree"]) <= 1e-9
        state_ok = all(x["state"][key] == job["expected"][key] for key in x["state"])
        charge_residual = abs(sum(x["charges_e"]) - job["expected"]["charge"])
        process = receipt["process"]
        admitted = (process["exit_code"] == 0 and process["error"] is None and process["stop_reason"] is None
                    and x["normal"] and x["converged"] and not x["printed_gate_failures"]
                    and state_ok and charge_residual <= 1e-5 and x["stdout_energy_difference_hartree"] <= 1e-9)
        assert admitted == receipt["passed"] == receipt["output"]["accepted"], f"admission differs: {job_id}"
        if job["variant"] not in ("default",):
            pnames = [n for n in job["files"] if n.endswith(".json")]
            assert len(pnames) == 1 and x["runtime_read_files"] and set(x["runtime_read_files"]) == set(pnames)
        else:
            assert not x["runtime_read_files"]
        rows[job_id] = {"accepted": admitted, "printed_gate_failures": x["printed_gate_failures"],
                        "state_matches": state_ok, "charge_sum_residual_e": charge_residual,
                        "receipt_sha256": sha(d / "receipt.json")}
        raw[job_id] = x
    pairs, parity, closed, replication, calcium_imports = [], [], [], [], []
    for job in manifests[0]["jobs"]:
        for key in ("archived_receipt", "archived_property"):
            assert sha(job[key]["path"]) == job[key]["sha256"]
        archive_directory = Path(job["archived_property"]["path"]).parent
        archived_receipt = read(job["archived_receipt"]["path"])
        assert archived_receipt["passed"]
        previous = parse_raw(archive_directory)
        current = raw[job["id"]]
        assert previous["state"] == current["state"]
        de = current["energy_hartree"] - previous["energy_hartree"]
        dq = max(abs(a-b) for a,b in zip(current["charges_e"],previous["charges_e"]))
        assert abs(de) <= 1e-6 and dq <= 1e-5
        replication.append({"id": job["id"], "energy_delta_hartree": de,
                            "max_charge_delta_e": dq, "archived_receipt": job["archived_receipt"],
                            "archived_property": job["archived_property"]})
    for job in manifests[2]["jobs"]:
        default_id = job["id"].removesuffix("corrected_file") + "default"
        current, default = raw[job["id"]], raw[default_id]
        assert current["energy_hartree"] == default["energy_hartree"]
        assert current["charges_e"] == default["charges_e"]
        calcium_imports.append({"id": job["id"], "default_id": default_id,
                                "energy_and_all_charges_identical": True})
    for job_id, job in expected.items():
        x = raw[job_id]
        if job["variant"] == "roundtrip":
            default_id = job_id.removesuffix("roundtrip") + "default"
            y = raw[default_id]
            pairs.append({"id": job_id, "energy_delta_hartree": x["energy_hartree"]-y["energy_hartree"],
                          "max_charge_delta_e": max(abs(a-b) for a,b in zip(x["charges_e"],y["charges_e"])),
                          "both_accepted": rows[job_id]["accepted"] and rows[default_id]["accepted"]})
        if job.get("stage") == "aquo50" and (job["variant"] == "corrected_tight" or job["metal"] in ("Ca", "Mg", "Sr", "Zn")):
            for engine, reference_id in job["upstream_accepted_ids"].items():
                ref = upstream[reference_id]
                assert ref["accepted"] and ref["engine"] == engine
                assert ref["input_sha256"] == job["files"]["context.xyz"]["sha256"]
                assert ref["charge"] == job["expected"]["charge"] and ref["multiplicity"] == job["expected"]["multiplicity"]
                assert len(x["charges_e"]) == len(ref["charges_e"])
                de = x["energy_hartree"]-ref["energy_hartree"]
                dq = max(abs(a-b) for a,b in zip(x["charges_e"],ref["charges_e"]))
                comparison = {"id": job_id, "engine": engine, "reference_id": reference_id,
                              "energy_delta_hartree": de, "max_charge_delta_e": dq,
                              "accepted_parity": rows[job_id]["accepted"] and abs(de) <= 1e-6 and dq <= 1e-5}
                (parity if job["variant"] == "corrected_tight" else closed).append(comparison)
    claimed = read(args.reported)
    accepted = sum(row["accepted"] for row in rows.values())
    assert accepted == claimed["accepted_native_endpoints"] == 52
    failures = [k for k,v in rows.items() if not v["accepted"]]
    assert set(failures) == set(claimed["failed_or_unaccepted"])
    assert len(parity) == 28 and all(r["accepted_parity"] for r in parity)
    result = {"schema": "independent_native_runtime_review.v2", "new_scientific_calls": 0,
              "workspace": str(args.workspace.resolve()), "reported_sha256": sha(args.reported),
              "upstream_data": upstream_pin, "independent_raw_parse_count": len(rows),
              "artifact_hash_checks": files_checked, "accepted": accepted, "failures": failures,
              "row_details": rows, "roundtrip_pairs": pairs, "corrected_vs_archived_upstream": parity,
              "La8_against_archived_native": replication, "Ca2_corrected_import_against_defaults": calcium_imports,
              "corrected_max_abs_energy_error_hartree": max(abs(x["energy_delta_hartree"]) for x in parity),
              "corrected_max_charge_error_e": max(x["max_charge_delta_e"] for x in parity),
              "closed_shell_vs_archived_upstream": closed,
              "scope": "Independent runtime/parameter-control recount; archived macOS xTB/tblite, no new upstream calls or protein validation."}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("row_details", "roundtrip_pairs", "corrected_vs_archived_upstream", "closed_shell_vs_archived_upstream", "La8_against_archived_native")}, indent=2))


if __name__ == "__main__":
    main()
