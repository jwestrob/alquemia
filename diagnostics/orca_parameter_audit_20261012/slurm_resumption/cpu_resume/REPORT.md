# Resume the existing Ca-source MACE batch on Slurm

This is a placement wrapper around the exact existing CPU batch snapshot and
manifest. It adds no optimizer, model, geometry or score change. Root owns the
graceful drain, submission and completion wake; this preparation made no
molecular calls, process signals or submissions.

## Execution choice

Use `standard-shared`, one node and 14 task slots. The unchanged two warm model
workers each use seven numerical threads; their optimizer parents compute while
their own worker is idle. Direct batch-shell execution avoids an `srun` task
binding that would leave the entire application one CPU. The runtime guard
requires exactly 14 allocated CPUs/tasks and at least 14 CPUs in shell affinity.
A larger whole-node CPU allocation or narrower affinity fails before MACE loads.

The script has no memory, CPUs-per-task, exclusive-allocation or time-limit flag.
The latest root instruction explicitly selects a shared CPU allocation for this
fixed two-lane workload. It does not change the general policy of fully using an
exclusive node. No 344-CPU node is reserved for fourteen active workers.

**Memory caveat:** actual `scontrol` inspection on 12 October shows cluster
`DefMemPerCPU=4000` but `standard-shared` has `DefMemPerNode=UNLIMITED`. Therefore
omitting memory flags does not prove a 56,000 MB reservation or actual RAM
sharing. The wrapper retains the complete job/node allocation and makes no such
claim. Root must inspect that actual record. The observed existing model peak
was 7,825,028 KiB per worker; that measurement is not a new resource cap.

The pinned historical warm GPU engine uses a different manifest/receipt contract,
and the current CPU worker hardcodes its device and child executable. A GPU
transport change is therefore separate work; root chose no GPU port in this turn.

## Drain and reuse

Required receipt:
`workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/DRAIN.json`.
It must name the exact original manifest and contain `status: "drained"`, with
empty `live_owned_processes`, `stopped_owned_processes` and
`incomplete_optimizer_task_ids` lists. Root determines these from real owned
processes and receipts. The wrapper cannot infer remote process liveness from
the compute node's local `/proc`.

Original input manifest SHA256:
`3b3c9c1fcd9c094790125a8354521954fbcb82973492c0b998b2386569fc44f8`.
Original CPU implementation SHA256:
`e3070a408777d8754506076440414f974ed01fc1207122bad812b88bad6718d0`.

The original shared `RUN.lock` and restart rules remain active. Every completed
source and compatible completed call is reused. Started incomplete optimizers
are not restarted. Each resumed attempt must pass the original four archived
Ca/La energy/force replays before any new-source calculation. These are actual
runtime checks, not assumed from this wrapper's syntax check.

## Commands and records

From the repository root, the completed no-compute dry-run is:

```bash
python diagnostics/orca_parameter_audit_20261012/slurm_resumption/cpu_resume/LAUNCH.py --manifest workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1/manifest.json --drain workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/DRAIN.json --receipts workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/cpu_resume --dry-run
```

Before the drain exists it explicitly reports that execution is not ready.
Python compilation and `bash -n` pass. No actual 14-slot Slurm placement or
compute-node MACE replay has been performed by this preparation agent.

Root's submission command after the verified drain is:

```bash
sbatch diagnostics/orca_parameter_audit_20261012/slurm_resumption/cpu_resume/resume.sbatch
```

Root must arm the existing completion-wake mechanism for the returned job ID.
Actual job/node/cpuset/thread environment, drain pin, frozen executable command
and wrapper costs are written under the explicit `cpu_resume` receipt directory.
That allocation evidence supersedes the frozen worker's historical
`scheduler_allocation=None` placement field only; scientific receipts remain
unaltered. The same existing `FINAL_COLLECTION.json` remains the solvent-stage
input after all 106 source records become terminal.
