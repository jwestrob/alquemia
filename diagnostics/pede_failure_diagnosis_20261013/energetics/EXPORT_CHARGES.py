#!/usr/bin/env python3
"""Export immutable charge JSON without the v1 metal-name/charge collision."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time


def pin(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


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
    config = json.loads(Path(args.config).read_text())
    if pin(__file__) != config["driver"] or pin(config["input"]["path"]) != config["input"]:
        raise SystemExit("changed input or export driver")
    output = Path(config["output"])
    output.mkdir(parents=True, exist_ok=False)
    receipt = {"job_id": job, "config": pin(args.config), "allocation": allocation,
               "new_molecular_calls": 0, "new_population_analysis": False}
    start = time.monotonic()
    try:
        saved = json.loads(Path(config["input"]["path"]).read_text())
        if saved["status"] != "complete" or len(saved["rows"]) != 36:
            raise ValueError("expected immutable completed 36-point analysis")
        rows = []
        for r in saved["rows"]:
            if r["metal"] not in ("Ca", "La") or r["medium"] not in ("vacuum", "water"):
                raise ValueError("origin metal/medium identity unavailable")
            for backend, b in r["backends"].items():
                rows.append({"source_case_id": r["source_case_id"], "cohort": r["cohort"],
                    "metal": r["metal"], "medium": r["medium"], "backend": backend,
                    "iterations": b["iterations"]["count"],
                    "last_delta_E_Eh": b["iterations"]["last"]["delta_energy_printed_Eh"],
                    "last_printed_residual": b["iterations"]["last"]["printed_residual"],
                    "charge_sum_e": b["charge_sum_e"],
                    **{"charge_e:" + k: v for k, v in b["groups_e"].items()}})
        if len(rows) != 72 or any(r["metal"] not in ("Ca", "La") for r in rows):
            raise ValueError("unexpected exported identity/denominator")
        target = output / "ORIGIN_CHARGES.tsv"
        with target.open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        receipt.update(status="complete", exported_rows=len(rows), output=pin(target),
                       superseded_export=config["superseded_export"], authoritative_JSON=config["input"])
    except BaseException as error:
        receipt.update(status="failed", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        receipt.update(elapsed_seconds=time.monotonic() - start,
                       process_CPU_seconds=usage.ru_utime + usage.ru_stime, peak_RSS_KiB=usage.ru_maxrss)
        with (output / f"EXECUTION_{job}.json").open("x") as handle:
            json.dump(receipt, handle, indent=2)
            handle.write("\n")


if __name__ == "__main__":
    main()
