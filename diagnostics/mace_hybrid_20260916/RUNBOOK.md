# MACE stage A: execution checkpoint

Approved 2026-09-16; active job **1200207**, queued on `gpu_h200` at this checkpoint.
No measured inference result yet. Eight real-artifact integrity/algebra tests
passed; these are not successful MACE integration tests.

**Execution unblocked:** Jacob explicitly accepted200000 MiB (~195 GiB):
“that's fine. let it run with the195GB.” [Resource amendment](RESOURCE_ACCEPTANCE.md).
The exact-share guard was replaced;1200197 was cancelled while still pending,
and1200207 requests the accepted allocation. Its task-owned continuationPID2628149
is active. The twelve scientific tasks and pinned implementation are unchanged.
The server-side rewrite mechanism remains unidentified; it no longer blocks this
approved initial allocation.
The [scheduler audit](SCHEDULER_MEMORY_AUDIT.md) also records that `CR_CPU`
does not track host-RAM reservation, and that controller logs are inaccessible.
That is a separate finding; it does not identify the cause of the200000MiB rewrite.

## Resources and pins

- Current allocation request: one H200,28 CPUs,**200000 MiB** host RAM (195.3125 GiB).
  Earlier proposed allocation: one H200, 28 CPUs, 257962 MiB (~251.916 GiB).
  Requested fraction: 2063701/8 MiB, rounded down by the scheduler's integral-MiB
  representation (0.625 MiB). This is one eighth to scheduler precision.
- Initial submission1200196 unexpectedly recorded 200000 MiB and seven days.
  Its pending memory edits did not persist. It was cancelled before allocation,
  with zero inference, and replaced by1200197 using explicit command-line resource
  arguments. The replacement recorded257962 MiB. A batch guard refuses inference
  if allocated memory differs from the prescribed share.
  `standard` QOS has a seven-day maximum: an unlimited-time request was pending
  for `QOSMaxWallDurationPerJobLimit`, so its external seven-day limit was restored.
  There is no project CPU/time stopping budget. No other jobs or priorities changed.
- Isolated Python3.11.15 environment, torch2.8.0, mace-torch0.3.16,
  graph-longrange0.4.4, e3nn0.4.4; complete dependency lock in the workspace.
- Backend source commit `6a86de5e3ed35fd86a55fc046aa085fe48a72764`.
- Official medium checkpoint SHA256
  `fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a`.
- Manifest SHA256
  `2fa2c928a3ecf5b7a8bcda76bf2c42014f81626983b118dcf1b72e70c8dab01c`.

Source-map audit: the archived full-system states differ only by binary floating
point roundoff <=7.11e-15 A on five Asp303 coordinates. Full MACE uses one common
copy; the frozen cores preserve their existing coordinates. No physical geometry
change. Core caps never enter the physical protein. All full pairs have 9141
atoms, La charge -8, Ca -9, closed-shell singlets. Installed source explicitly
interprets the `spin=1` input as multiplicity (`total_spin - 1`).

## Read-only status and collection

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MACE_WORK="$PWD/workspaces/mace_hybrid_20260916"
MACE_PY="$MACE_WORK/software_v1/venv/bin/python"
MACE_MANIFEST="$MACE_WORK/pilot_v1/manifest.json"
squeue -j 1200207
"$MACE_PY" scripts/mace_hybrid.py dry-run --manifest "$MACE_MANIFEST"
"$MACE_PY" scripts/mace_hybrid.py collect --manifest "$MACE_MANIFEST"
"$MACE_PY" -m unittest discover -s tests -p test_mace_hybrid.py -v
```

Raw task products, inputs, snapshots, timing, failed attempts and software are
under `workspaces/mace_hybrid_20260916/`; no production scorer is modified.
Batch exits collect `pilot_v1/collection_job_JOBID.json`. The separate task-owned
continuation uses the existing `scripts/affordable_watch.py` accounting watcher
and writes `pilot_v1/continuation/`. Its `completion.json` is the terminal record.
It never changes inputs, model, precision, electrostatic boundary or scientific
task list. Verified successful tasks are reused; the declared repeat tasks have
distinct identities and execute independently.

## Authorized technical recovery

The continuation first tries host offload if native autograd runs out of GPU
memory. Host OOM may request 2, 4, then 8 GPU shares with proportional CPUs/RAM.
Each still executes one model task at a time on one GPU. GPU memory does not
pool; GPU OOM after offload requires a kernel/implementation review, not more
unused GPU reservations. Unknown failures also stop for technical inspection.
All new launches and terminal resource receipts are retained. Do not start a
second continuation or remove its flock while it is active.

For explicit operator recovery **only after the owned job and continuation are
terminal**, the same manifest can be resumed without rerunning valid tasks:

```bash
sbatch --cpus-per-task=28 --gres=gpu:1 --mem=200000M --time=7-00:00:00 \
  --output="$MACE_WORK/pilot_%j.out" --error="$MACE_WORK/pilot_%j.err" \
  diagnostics/mace_hybrid_20260916/run_pilot.sbatch \
  "$MACE_PY" "$MACE_MANIFEST" native
```

Inspect the actual scheduler memory request. The first share accepts200000 MiB.
A host-OOM recovery must actually receive more RAM than the failed attempt;
reserving extra GPUs with the same RAM cap is not a successful resource increase.
A different memory mode or allocation is an explicit launch argument and recorded
in every task receipt.

## Interpretation

Energy is in eV; forces are negative energy gradients in eV/A. Reuse DFT in
Hartree, subtract before conversion, and report components separately. The
subtractive expression is in PLAN.md. There is no compatible reference or
solution-phase calibration, so S and class stay unavailable. Direct MACE,
hybrid partition behavior, runtime and physical checks answer different
questions. Baseline remains unchanged regardless of this pilot's result.
