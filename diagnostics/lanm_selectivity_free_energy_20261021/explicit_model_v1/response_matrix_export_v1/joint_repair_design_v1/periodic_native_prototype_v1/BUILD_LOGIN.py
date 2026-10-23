#!/usr/bin/env python3
"""Build the isolated, pinned CPU prototype; run no molecular calculation."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import socket
import sys
import time


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


manifest_path, expected_hash, out_arg = sys.argv[1:]
manifest_path = Path(manifest_path).resolve()
out = Path(out_arg).resolve()
out.mkdir(exist_ok=False, parents=True)
start = time.monotonic()
start_cpu = time.process_time()
result = {"status": "failed", "manifest": pin(manifest_path), "commands": [],
          "new_molecular_calls": 0, "job_id": None, "execution_mode": "login",
          "hostname": socket.gethostname(), "controller_pid": os.getpid(),
          "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
try:
    if result["manifest"]["sha256"] != expected_hash:
        raise ValueError("Source manifest changed")
    manifest = json.loads(manifest_path.read_text())
    for item in manifest["source_outputs"]:
        if pin(item["path"]) != item:
            raise ValueError("Source pin changed: " + item["path"])
    threads = 1
    if os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_STEP_ID"):
        raise ValueError("Truthful login build requires no Slurm identity")
    if socket.gethostname() != "biotite":
        raise ValueError("Declared login host is biotite")
    if any(os.environ.get(k) != "1" for k in
           ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS")):
        raise ValueError("One numerical thread per compiler process")
    source = Path(manifest["source_root"])
    if pin(source / "cmake" / "CMakeLists.txt") not in manifest["source_outputs"]:
        raise ValueError("Build-system source not pinned")
    commands = [
        ["cmake", "-S", str(source / "cmake"), "-B", str(out / "build"),
         "-DTINKER_DIR=" + str(source), "-DCMAKE_BUILD_TYPE=RELWITHDEBINFO",
         "-DCMAKE_Fortran_COMPILER=/usr/bin/f95",
         "-DCMAKE_Fortran_FLAGS_RELWITHDEBINFO=-O2 -g -fbacktrace",
         "-DFFTW_LIB=/usr/lib/x86_64-linux-gnu/libfftw3.a",
         "-DFFTW_THREADED_LIB=/usr/lib/x86_64-linux-gnu/libfftw3_threads.a",
         "-DTINKER_OPENMM=OFF"],
        ["cmake", "--build", str(out / "build"), "--target", "analyze", "testgrad",
         "--parallel", str(threads)],
    ]
    for index, command in enumerate(commands):
        log = out / f"command_{index}.log"
        with log.open("w") as stream:
            process = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
        result["commands"].append({"argv": command, "returncode": process.returncode,
                                   "log": pin(log)})
        if process.returncode:
            raise RuntimeError("Build command failed; see " + str(log))
    result.update(status="complete", threads=threads,
                  executables={name: pin(out / "build" / name) for name in ("analyze", "testgrad")},
                  library=pin(out / "build" / "libtinker.a"),
                  cache=pin(out / "build" / "CMakeCache.txt"))
except Exception as error:
    result["error"] = repr(error)
finally:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    result.update(wall_seconds=time.monotonic() - start,
                  controller_CPU_seconds=time.process_time() - start_cpu,
                  child_cpu_seconds=usage.ru_utime + usage.ru_stime,
                  peak_child_rss_KiB=usage.ru_maxrss,
                  scheduler_allocated_CPU_seconds=0, scheduler_allocation=None,
                  maximum_parallel_build_workers=1)
    (out / "RESULT.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({key: result[key] for key in ("status", "wall_seconds", "new_molecular_calls")}))
sys.exit(0 if result["status"] == "complete" else 1)
