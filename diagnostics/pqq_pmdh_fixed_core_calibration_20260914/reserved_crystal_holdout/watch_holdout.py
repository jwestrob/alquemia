#!/usr/bin/env python3
"""Persist terminal SLURM state, then wake the exact approved Codex thread."""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
import subprocess
import time


HERE = Path(__file__).resolve().parent
TERMINAL = {
    "BOOT_FAIL",
    "CANCELLED",
    "COMPLETED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "SPECIAL_EXIT",
    "TIMEOUT",
}


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def job_state(job_id: str) -> tuple[str, str]:
    completed = subprocess.run(
        ["sacct", "-j", job_id, "--format=JobIDRaw,State,ExitCode", "-n", "-P"],
        capture_output=True,
        text=True,
        check=False,
    )
    rows = [line.split("|") for line in completed.stdout.splitlines() if line]
    exact = next((row for row in rows if row[0] == job_id), None)
    if completed.returncode != 0 or exact is None:
        return "UNKNOWN", ""
    return exact[1].split()[0].rstrip("+"), exact[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("job_id")
    parser.add_argument("thread_id")
    args = parser.parse_args()

    watch_receipt = HERE / f"watch_{args.job_id}.json"
    execution_receipt = HERE / "prepared" / f"holdout_execution_{args.job_id}.json"
    result = HERE / "result" / "holdout_result.json"
    while True:
        state, exit_code = job_state(args.job_id)
        atomic_json(
            watch_receipt,
            {
                "job_id": args.job_id,
                "state": state,
                "exit_code": exit_code,
                "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            },
        )
        if state in TERMINAL:
            break
        time.sleep(120)

    message = (
        f"AUTOMATIC PQQ CRYSTAL HOLDOUT STATUS: SLURM job {args.job_id} ended "
        f"{state} (exit {exit_code}). Read {watch_receipt} and {execution_receipt}. "
        f"If present, read {result}; report the frozen-band primary PASS/FAIL and "
        "both 1H4I/4MAE S scores, then commit result receipts and update the vault."
    )
    command = [
        "/home/jwestrob/.local/bin/codex",
        "queue",
        "--thread",
        args.thread_id,
        "--message",
        message,
    ]
    for attempt in range(1, 13):
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        payload = json.loads(watch_receipt.read_text())
        payload["notification"] = {
            "attempt": attempt,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "command": command,
        }
        atomic_json(watch_receipt, payload)
        if completed.returncode == 0:
            return
        time.sleep(300)
    raise SystemExit("unable to queue Codex notification after 12 attempts")


if __name__ == "__main__":
    main()
