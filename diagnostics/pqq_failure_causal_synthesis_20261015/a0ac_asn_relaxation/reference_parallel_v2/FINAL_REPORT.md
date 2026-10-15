# DFT supports the large La benefit, with a Ca discrepancy

All four fixed endpoints completed and pass the unchanged energy-diagnostic
checks. Native r2SCAN-3c independently supports the large La stabilization and
the direction of the MACE-predicted selected-state response. It does **not**
support the small Ca stabilization: the new Ca geometry is higher in DFT energy.
Retain this discrepancy rather than describing the two models as equivalent.

This closes the declared four-cell local-energy test. The full uniform benchmark
transfer remains pending; no additional DFT or production promotion follows.

## Actual response

Values below are kcal/mol, from the executed collector, rounded only for display.
Unrounded energies, differences, source hashes and checks are in [RESULT.json](RESULT.json).

| New selected minus old selected | Native DFT | MACE | DFT minus MACE |
|---|---:|---:|---:|
| Ca | +1.906580 | −0.286036 | +2.192615 |
| La | −17.400253 | −18.274380 | +0.874127 |
| Ca work minus La work | +19.306832 | +17.988344 | +1.318488 |

The La discrepancy lies within the predeclared 1 kcal/mol interpretive resolution.
The Ca and paired discrepancies exceed it. This resolution is neither an admission
threshold nor an error bar for DFT. A positive paired response alone is not a
qualification test.

The main gain is therefore supported independently of MACE for this conditional
local representation. The small Ca gain and exact magnitude are not. Relaxation
moves several donors: this test cannot assign the response solely to Asn oxygen
orientation or removal of one H–H contact. The earlier two-H sensitivity used
origin-based paths, not these selected old-to-new paths, and does not bound this
increment's hydrogen sensitivity.

## Scope and numerical evidence

One consumed XoxF1 A0ACD6B9F2 source: Ca-conditioned seed1/sample4,
envelope `2eda778295ffed58`; 204 atoms, 14 mapped caps, dry oxidized PQQ3−,
original proton inventory. Ca−4 and La−3 are singlets. Each metal compares its
own saved `adaptive_M` geometry against `relaxed_flip_M`. No geometry was
optimized with DFT and no new geometry was chosen using these DFT energies.

Protocol: `A0AC_Ca4_selected_old_new_native_r2SCAN3c_vacuum_v1`;
ORCA 6.1.1, r2SCAN-3c/DefGrid3/TightSCF/NoAutostart, native def2-mTZVPP and
auxiliary basis, La 46-electron ECP. No GFN2 parameter file or solvent term is
part of this calculation. The staging-v2 revision changed scheduling only.

All four cells pass all 20 fixed checks and terminate normally. Their last
TRAH/NR orbital-gradient residuals are below the actual printed 1e−5 threshold:
old Ca 5.627194e−7, old La 8.035475e−7, new Ca 9.834354e−7, new La 4.762669e−7.
Large retained density/DIIS summary entries have unresolved post-switch provenance;
near-final frontier and electronic-Hessian warnings remain. Independent electronic
stability, agreement between electronic branches and force accuracy were not tested.
See [the convergence interpretation](../../review/REPORT.md#how-to-interpret-the-native-reference-checks).

These are distinct metal-specific paths, not a common displacement or a complete
cross-scored DFT candidate pool. They do not establish DFT pool minima, solvent
accuracy, basin populations, an affinity, a class decision or biological transfer.
In particular, DFT's higher new Ca energy must not be silently replaced by the
old Ca energy and presented as the declared four-cell response.

## Execution and provenance

Workspace: `workspaces/pqq_failure_causal_synthesis_20261015/asn_selected_reference_parallel_v2`.
Manifest SHA256: `08a9b8965e0749e1ff9381d3e524fdf150c38835ce5a185403019ff58371878e`.
Final collection: `attempts/5868_collect_all_four/COLLECTION.json`, SHA256
`4f99815e0e94421dd9cfeb64c8cd8faa74b2b3338aadee74a6119ef739e2658f`.
The matching application EXECUTION receipt is complete. Exact launch commands,
decisions and watcher commands remain in workspace SUBMISSION files; no rerun is
needed to reproduce extraction from these retained results.

| Endpoint | Job | Scheduler wall seconds | Allocated CPUs |
|---|---:|---:|---:|
| Old Ca | 5762 | 6,736 | 12 |
| Old La | 5763 | 7,026 | 12 |
| New Ca | 5866 | 8,073 | 12 |
| New La | 5867 | 7,628 | 12 |

Four molecular evaluations cost **353,556 allocated CPU-seconds**, zero GPU.
Preflight5761 adds 22 seconds on 12 CPUs (264 allocated CPU-seconds), zero molecular
calls. These are development-validation costs, not routine classifier latency or
measured utilization. Collector5868 wrote its complete application receipt after
23.33 seconds. Its scheduler allocation lingered after all application work had
finished. A same-allocation process probe failed at Slurm task launch; the exact
shutdown/communication cause is unresolved. Root cancelled only this completed
zero-molecular collector allocation, preserving its complete results. Scheduler
state is **CANCELLED by 601**, 363 seconds on 12 CPUs (4,356 allocated CPU-seconds),
not normal scheduler completion. The whole four-cell stage including preflight
and collector consumed **358,176 allocated CPU-seconds**, zero GPU. See workspace
`COLLECTOR_ALLOCATION_RELEASE.json` and `COLLECTOR_SCANCEL_RECEIPT.json`.

Independent read-only review confirmed the four source/output joins, unchanged
admission, and the collector's same-metal new-minus-old then Ca-minus-La algebra.
No calculation, tolerance, band or experimental label was changed after inspection.

**Next decision:** interpret the already-running full transfer5661 and its queued
report5666/figure5748. The present result strengthens the amide-search mechanism;
it does not establish that uniform application improves the whole benchmark.
