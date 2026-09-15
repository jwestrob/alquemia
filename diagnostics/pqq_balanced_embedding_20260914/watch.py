#!/usr/bin/env python3
"""Persist terminal SLURM state and wake the exact Codex thread once."""

import argparse
import datetime
import json
from pathlib import Path
import subprocess
import time


HERE = Path(__file__).resolve().parent
TERMINAL = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
    "NODE_FAIL", "BOOT_FAIL", "PREEMPTED",
}


def atomic(path: Path, value: object) -> None:
    temporary = path.with_name("." + path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def state(job: str) -> tuple[str, str]:
    result = subprocess.run(
        ["sacct", "-j", job, "--format=JobIDRaw,State,ExitCode", "-n", "-P"],
        capture_output=True,
        text=True,
        check=False,
    )
    rows = [line.split("|") for line in result.stdout.splitlines() if line]
    exact = next((row for row in rows if row[0] == job), None)
    if result.returncode or exact is None:
        return "UNKNOWN", ""
    return exact[1].split()[0].rstrip("+"), exact[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("job")
    parser.add_argument("thread")
    args = parser.parse_args()
    receipt = HERE / "watch_receipt.json"
    while True:
        status, exit_code = state(args.job)
        atomic(
            receipt,
            {
                "job_id": args.job,
                "state": status,
                "exit_code": exit_code,
                "checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            },
        )
        if status in TERMINAL:
            break
        time.sleep(120)
    message = (
        f"AUTOMATIC BALANCED-EMBEDDING STATUS — SLURM job {args.job} ended {status}. "
        f"Read {receipt}; if complete, read {HERE / 'result.json'} and report the frozen gates."
    )
    command = [
        "/home/jwestrob/.local/bin/codex", "queue", "--thread", args.thread,
        "--message", message,
    ]
    for attempt in range(1, 13):
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        payload = json.loads(receipt.read_text())
        payload["notification"] = {
            "attempt": attempt,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "command": command,
        }
        atomic(receipt, payload)
        if result.returncode == 0:
            return
        time.sleep(300)
    raise SystemExit("unable to queue Codex notification after 12 attempts")


if __name__ == "__main__":
    main()
