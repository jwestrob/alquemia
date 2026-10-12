"""Run the unchanged six-cell comparator once, inside an actual Slurm job."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import time
import traceback


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def verify(records):
    for item in records:
        if pin(item["path"]) != item:
            raise RuntimeError("Changed pinned input: " + item["path"])


def write_new(path, data):
    with Path(path).open("x") as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write("\n")


def interrupted(signum, frame):
    raise RuntimeError("Interrupted by signal " + str(signum))


def main(plan_path):
    # No scientific imports or subprocesses are permitted before this guard.
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Comparison execution requires Slurm")
    plan = json.loads(plan_path.read_text())
    output = Path(plan["output_directory"])
    receipt_path = output / "RECEIPT.json"
    result_path = output / "RESULT.json"
    if receipt_path.exists() or result_path.exists():
        raise RuntimeError("Existing comparison result/receipt; do not repeat")
    receipt = {"status": "running", "plan": pin(plan_path),
               "slurm_job_id": os.environ["SLURM_JOB_ID"], "host": socket.gethostname(),
               "allocated_cpus": os.environ.get("SLURM_CPUS_ON_NODE"),
               "started_epoch": time.time(), "new_molecular_calls": 0,
               "result_path": str(result_path), "returncode": None}
    write_new(output / "START.json", receipt)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    code = 1
    try:
        if int(os.environ.get("SLURM_CPUS_ON_NODE", "0")) != 1:
            raise RuntimeError("Expected the declared one-CPU comparison allocation")
        verify(plan["pins"])
        if receipt["plan"] != pin(plan_path):
            raise RuntimeError("Plan changed during startup")
        observed = [{"path": item["path"], "upstream_job": item["upstream_job"],
                     "available": Path(item["path"]).is_file(),
                     "pin": pin(item["path"]) if Path(item["path"]).is_file() else None}
                    for item in plan["future_collections"]]
        receipt["future_collections_at_start"] = observed
        query = subprocess.run(["sacct", "-n", "-P", "-j", "2665,2666",
                                "--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS"],
                               capture_output=True, text=True)
        receipt["upstream_accounting"] = {"returncode": query.returncode,
                                          "stdout": query.stdout, "stderr": query.stderr}
        command = [plan["python"], plan["comparator"], "--config", plan["config"],
                   "--output", str(result_path)]
        receipt["command"] = command
        env = os.environ.copy()
        env.update({name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS")})
        with (output / "comparison.log").open("x") as log:
            process = subprocess.run(command, cwd=plan["repository"], env=env,
                                     stdout=log, stderr=subprocess.STDOUT)
        receipt["returncode"] = process.returncode
        if process.returncode:
            raise RuntimeError("Existing comparator failed; inspect comparison.log")
        verify(plan["pins"])
        for item in observed:
            if Path(item["path"]).is_file() != item["available"]:
                raise RuntimeError("Collection availability changed during comparison")
            if item["available"]:
                verify([item["pin"]])
        if receipt["plan"] != pin(plan_path):
            raise RuntimeError("Plan changed during comparison")
        result = json.loads(result_path.read_text())
        receipt.update(status="comparison_complete", result=pin(result_path),
                       complete_cells=result["complete_cells"],
                       declared_cells=result["declared_cells"],
                       unavailable_cells=[{"source": row["source"], "metal": row["metal"],
                                           "status": row["status"], "reason": row.get("reason")}
                                          for row in result["rows"] if row["status"] != "complete"])
        code = 0
    except Exception as error:
        receipt.update(status="failed", error=str(error), traceback=traceback.format_exc(),
                       result_accepted=False)
    finally:
        receipt["finished_epoch"] = time.time()
        receipt["elapsed_seconds"] = receipt["finished_epoch"] - receipt["started_epoch"]
        receipt["allocated_cpu_seconds_through_collection"] = receipt["elapsed_seconds"]
        write_new(receipt_path, receipt)
    return code


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(main(args.plan.resolve()))
