# MACE stage A: completed memory and capability pilot

**Updated 2026-09-16.** All 12 approved medium-model calls completed on A5000
in job 1200381. Full 9,141-atom evaluations take about 58 seconds each and
9.55 GiB GPU allocation, with at most 1.65 GiB worker host RSS. Four original-core
comparisons and all 17 tests pass. Repeats, translations and charge closure pass;
rotation and the approximately -5 kcal/mol hybrid partition shift fail frozen
physical tolerances. See [final results](MEMORY_RESULTS.md).

The baseline/default is unchanged. No calibrated MACE class, solution-phase
reference or predictive improvement is established. Large-checkpoint execution
has not run. There is no active memory-work job or H200 continuation.

## Current campaign and pins

- Campaign: `workspaces/mace_hybrid_20260916/blocked_v6`.
- Manifest SHA256: `8d5fc630c11b50457cd30927f9bcb827f99bd1c67191b307ccb01bc0712c7189`.
- Medium checkpoint SHA256: `fab8b8713c832f31a2a853aaa22fd638be8a369cbf5095e6b3e982a18d10e93a`.
- Isolated Python 3.11.15, torch 2.8.0, mace-torch 0.3.16, graph-longrange 0.4.4,
  e3nn 0.4.4; dependency lock and backend source inventory in `software_v1`.
- GPU allocation: one RTX A5000, 16 CPUs, 64,474 MiB requested host RAM,
  partition `gpu`, node `node-128-512g-8gpu-1`. Execution mode `native`.
- Scientific protocol: `mace_polar_1m_vacuum_r2scan3c_subtractive_pilot_v1`.
- Implementation: isolated interface adapter plus pair/edge/node memory blocking;
  [design](BLOCKED_KERNEL.md), [agreement](MEMORY_AGREEMENT.md).

Both full endpoints contain 9,141 atoms, PQQ 3-, no waters, singlet multiplicity,
La total charge -8 and Ca -9. They use a common source geometry. The archived
cores retain their exact source coordinates; caps do not enter the full protein.
The source audit found only <=7.11e-15 Angstrom floating-point differences on
five Asp303 coordinates between the original archived states. No geometry search
or protonation change was introduced. Four archived vacuum DFT endpoints are
reused; new DFT endpoints: zero.

## Read-only verification and collection

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MACE_WORK="$PWD/workspaces/mace_hybrid_20260916"
MACE_PY="$MACE_WORK/software_v1/venv/bin/python"
MACE_MANIFEST="$MACE_WORK/blocked_v6/manifest.json"
MACE_RUNNER="$MACE_WORK/blocked_v6/implementation/mace_hybrid.py"
"$MACE_PY" "$MACE_RUNNER" dry-run --manifest "$MACE_MANIFEST"
"$MACE_PY" "$MACE_RUNNER" compare-cores --manifest "$MACE_MANIFEST"
"$MACE_PY" "$MACE_RUNNER" collect --manifest "$MACE_MANIFEST"
```

The terminal collection is `blocked_v6/collection_job_1200381.json`; its
human-readable report is `blocked_v6/REPORT.md`. Exact resource accounting,
including failures, is in `memory_cost_final.json` and `memory_sacct_final.tsv`.
Per-task receipts retain forces, density coefficients, model/adapter metadata,
GPU allocation/reservation, worker RSS and elapsed time. Accounting watcher
receipt: `blocked_v6/job_1200381_accounting.json`.

## Recovery of this exact manifest

No recovery is presently needed. The runner reuses valid tasks and preserves
failed attempts; its four-core gate prevents unverified full-system execution.
If recovering an interrupted copy of this approved campaign, first verify that
its executor is terminal and its pinned inputs and implementation match. The
original launch is recorded in `blocked_v6/full_launch.json`. With the variables
above, the resource-equivalent command is:

```bash
sbatch --partition=gpu --nodelist=node-128-512g-8gpu-1 \
  --cpus-per-task=16 --gres=gpu:1 --mem=64474M --time=7-00:00:00 \
  --export=ALL,MACE_MIN_MEMORY_MIB=64474 \
  --output="$MACE_WORK/blocked_recovery_%j.out" \
  --error="$MACE_WORK/blocked_recovery_%j.err" \
  diagnostics/mace_hybrid_20260916/run_pilot.sbatch \
  "$MACE_PY" "$MACE_MANIFEST" native
```

Seven days is the cluster QOS limit, not a project compute stopping budget.
Do not automatically offload or reserve extra shares for a GPU OOM: inspect the
specific tensor allocation. The successful implementation needs neither. The
requested RAM is not claimed to be an enforced per-process RSS cap.

## Historical campaigns and interpretation

`pilot_v3` contains the original four successful core calls (job 1200308) and
is the pinned numerical reference; [core result](CORE_RESULTS.md). Its queued
H200 job 1200309 was cancelled without allocation after the first blocked-core
validation. Its continuation was stopped. Earlier H200 resource negotiation is
preserved in [acceptance](RESOURCE_ACCEPTANCE.md) and
[scheduler audit](SCHEDULER_MEMORY_AUDIT.md); do not relaunch that obsolete queue
entry. `blocked_v1` through `blocked_v5` retain every memory-debug failure.

`E_hybrid = E_MACE(full) + E_DFT(core) - E_MACE(core)` is a vacuum diagnostic.
Energy is in eV; forces are negative gradients in eV/Angstrom. Archived DFT
energies are Hartree; subtract before conversion. The full term cancels in the
partition contrast. No compatible aquo reference exists, so S and class remain
null. Direct scoring, hybrid partition behavior, numerical consistency, runtime
and biological accuracy are separate questions.

Proposed follow-up: locate the rotation sensitivity and verify large-checkpoint
compatibility on the same cores. Neither is an executed result or authorization
for additional scientific variants. Agree consequential method changes with Jacob.
