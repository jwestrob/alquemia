# MACE core checks: completed 2026-09-16

**Four real MACE-POLAR-1 medium core evaluations completed on RTX A5000.**
The vacuum subtractive hybrid's partition shift is **-4.998674456 kcal/mol**:
substantial cancellation of the archived61.738603357 kcal/mol DFT-core shift,
but **FAIL** against the frozen2 kcal/mol diagnostic. This is not evidence of
improved La/Ca classification. Baseline/defaults are unchanged.

## Exact experiment and outputs

Consumed1H4I chain A, archived qm33/qm36 cores (47/54 atoms), original charges,
singlets, PQQ state and coordinates; no waters, solvent or external field.
Float64, isolated realspace, unchanged medium weights. Zero new DFT endpoints.
Interface compatibility defects were repaired in a versioned adapter; see
[repair](REALSPACE_INTERFACE_REPAIR.md). No numerical electrostatic kernel changed.

| Core | Metal | MACE E (eV) | Evaluation seconds | Peak CUDA allocated bytes |
|---|---|---:|---:|---:|
| qm33 | La | -47764.98110144369 | 1.9296248331665993 | 1145917952 |
| qm33 | Ca | -65345.66162812736 | 0.8742979764938354 | 1145917952 |
| qm36 | La | -53984.178300211475 | 0.9167513959109783 | 1283105280 |
| qm36 | Ca | -71561.96482462324 | 0.8138684332370758 | 1283105280 |

Charge sums pass: maximum absolute error1.7763568394002505e-15 e.
Peak GPU allocation1.067–1.195 GiB; maximum worker host RSS1626912 KiB
(1.552 GiB). Evaluation times include energies and analytic forces; separate
worker processes take7.14–10.93 seconds each including imports/model setup.
The35-second allocation covers all four tasks and batch overhead.

## Partition result without the full-system calculation

`R = E_Ca - E_La`, energies converted exactly once.

| Component | qm33 (kcal/mol) | qm36 (kcal/mol) | qm36 minus qm33 |
|---|---:|---:|---:|
| Archived DFT core R | -405349.0793037981 | -405287.34070044087 | 61.73860335722566 |
| MACE core R | -405420.12418042123 | -405353.3869026081 | 66.73727781313937 |
| DFT minus MACE R | 71.04487662314204 | 66.04620216722833 | -4.998674455913715 |

For `E_hyb(M)=E_MACE(full,M)+E_DFT(core,M)-E_MACE(core,M)`, the identical
full-system term cancels exactly between partitions. Thus the hybrid partition
shift is already determined, although full-system energies, runtime, direct R
and hybrid absolute R remain unavailable. Do not interpret the large elemental
offsets or correction signs as calibrated affinity. No aquo reference/S/class.
This does not identify which method best describes the added residue or caps.

## Actual cost and verification

| Job | Result | Allocated wall seconds | CPUs | Allocated CPU-seconds |
|---|---|---:|---:|---:|
|1200302|Unnecessary reciprocal-grid OOM; no energy|23|16|368|
|1200306|Missing checkpoint metadata; no energy|11|16|176|
|1200308|All four core tasks computed|35|16|560|

Total69 GPU allocation seconds and1104 allocated CPU-seconds, including both
failed attempts. All jobs requested one A5000 and64474 MiB host RAM. Slurm reports
75.820 summed actual CPU-seconds across the three parent job records; accounting
and process-level memory receipts are retained. These development failures are
not part of an ordinary successful four-core inference cost. Queued/cancelled
1200207 and metadata probe1200301 used zero allocation time.

Twelve tests pass: archived energy/sign/unit replay; exact input/state and cap
mapping; transform inventory; cache rejection/invalidation; versioned retry
inputs; bitwise checkpoint preservation; adapter vs original realspace feature,
energy and gradient equality on all four actual computed densities. No fabricated
MACE outputs. The four GPU calls are the actual scientific integration runs.
Full-protein repeat/rotation/translation checks have **not run** at this checkpoint.

## Provenance and continuation

- Scientific protocol remains `mace_polar_1m_vacuum_r2scan3c_subtractive_pilot_v1`.
- Execution adapter `polar0316_graph044_isolated_interface_v2`.
- Current manifest `workspaces/mace_hybrid_20260916/pilot_v3/manifest.json`, SHA256
  `eff0b6bbe28d045d7e4fcf904dfd302d72fda67d50148d45c0cf3b7ba8414dee`.
- Collection: `pilot_v3/collection_job_1200308.json`; task receipts under
  `pilot_v3/execution/`; every failed attempt remains under v1/v2.
- Test log: `workspaces/mace_hybrid_20260916/interface_and_core_tests_v3.log`.
- Previously approved remaining eight full-system calls queued as1200309:
  one H200,28 CPUs,200000 MiB RAM. It reuses the four successful core receipts.
  Task-owned continuationPID3463060; memory-only recovery remains authorized.

Recommendation: continue the already approved full-system feasibility checks;
retain the baseline. Core execution is affordable, hybrid partition consistency
misses its target, and predictive usefulness remains untested. Neither global
affordability nor direct MACE discrimination follows from these core results.
Commands and current status: [runbook](RUNBOOK.md).
