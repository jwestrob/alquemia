# Remaining local LBT reference endpoints — complete

All three selected endpoints completed and passed the unchanged native
state/energy/analytic-gradient admission. Each executed exactly one new PBE0-D4
vacuum EnGrad call. The original 7CCO/La scout 12826 was reused, not rerun.
No source contrasts or method comparisons were calculated for this report;
the separate saved-data comparator owns that calculation.

| Endpoint | Job | SCF cycles | Effective / physical multiplicity | ECP core |
|---|---:|---:|---|---:|
| 7CCN/La | 12844 | 17 | 1 / 1 | 46 |
| 7CCO/Dy | 12845 | 18 | 1 / 6 | 55 |
| 7CCN/Dy | 12846 | 18 | 1 / 6 | 55 |

The native outputs confirm charge −1, 200 explicit electrons, molecular
multiplicity 1, the expected metal ECP and normal ORCA 6.1.1 termination. The
archived 48-atom inputs and prepared cap maps remain pinned. Each receipt confirms
the expected complete gradient file, and each selected result retains 48 gradient
vectors, explicitly `gradient_not_force`. Dy is the frozen 4f9 effective-singlet
model with separate physical-sextet metadata, not an explicit-f sextet or
spin–orbit calculation. The atomic guess's printed multiplicity does not change
the molecular state's effective multiplicity.

## Convergence evidence and limits

All three use native `ConvCheckMode=Total+1el-Energy`, `ConvForced=0`, with SOSCF
initialized. Energy and orbital-gradient residuals pass their printed tolerances.
MAX/RMS density and DIIS residuals do not; 7CCN/La's orbital rotation also exceeds
its printed tolerance, while both Dy orbital rotations pass. Consequently,
**all three `all_printed_residuals_pass` flags remain false**. These observations
do not replace the predeclared native admission rule.

Each output includes a small-gap warning while skipping pre-diagonalization
(reported 0.039 for 7CCN/La and 7CCO/Dy, 0.038 for 7CCN/Dy). These are the actual
SCF warning messages, not independently validated final-state gaps. No negative
final gap, electronic stability, numerical force consistency or physical-state
qualification is inferred. `force_consistency_qualified` and
`physical_state_qualified` remain false throughout.

The singleton results retain all four declared slots but only one newly selected
task per allocation. Their `full_four_available=false` is not a new failure of
the other completed endpoints; the four-cell collector must join their owning
results and the unchanged scout explicitly.

## Actual scheduler cost and preserved failures

| Job / role | State | Elapsed seconds | CPUs | Allocated CPU-seconds | Measured TotalCPU | Batch MaxRSS |
|---|---|---:|---:|---:|---|---|
| 12844, 7CCN/La | COMPLETED | 1068 | 12 | 12816 | 03:29:21 | 3397672K |
| 12845, 7CCO/Dy | COMPLETED | 866 | 12 | 10392 | 02:50:32 | 3108564K |
| 12846, 7CCN/Dy | COMPLETED | 995 | 12 | 11940 | 03:16:16 | 3242888K |
| 12837, original precall failure | FAILED | 2 | 12 | 24 | 00:00.854 | 63508K |
| 12838, original precall failure | FAILED | 3 | 12 | 36 | 00:01.145 | 112544K |
| 12839, original precall failure | FAILED | 1 | 12 | 12 | 00:00.983 | 76716K |
| 12840, saved replay diagnosis | COMPLETED | 2 | 1 | 2 | 00:00.922 | 110976K |

Each successful endpoint used 12 MPI slots and the actual 48000 MiB grant, with
grant-derived MaxCore 3100 MB per rank. Batch MaxRSS is scheduler reporting, not a
claim of combined MPI peak memory; the adapter's parent RSS is also not that peak.
No GPU, new MACE/GFN calls, optimization or changed scientific settings occurred.

The original 12837–12839 records remain `execution_failure`, each with zero
molecular calls. Diagnosis 12840 localized their replay failure to JSON tuple/list
serialization of `executed_ecp`. V2 normalized that representation while retaining
the original source, chemistry, numerical settings and admission semantics.
These failures and their costs are preserved separately, not relabeled as SCF
failures or erased by the successful replacement jobs.

## Exact result pins

All paths/pins, raw scheduler query, selected SCF/state evidence and original
failure identities are retained in `remaining_result_pins.json`.

| Endpoint | Owning `RESULT.json` SHA256 |
|---|---|
| 7CCN/La | `c7d9d92d6225775dd2ccb323852ec6bad6312fa2b6ba7a3c682bf974ac556558` |
| 7CCO/Dy | `00d810c5e1940ebf1b5ff02ec8e478eb19e434ff120007e277a97b8aa41f47b7` |
| 7CCN/Dy | `1db006eca709e0af38de8df50936a025bc05d0d1ec99d5d7919275b0ec722b9d` |

These are under
`workspaces/lanm_coupled_response_20261016/lbt_local_reference_remaining_v2/<task>/execution/`.
All three use source SHA `10cd06ec4d1d569adfba5c56b3b37b5f8a7e61627dff5d591b41883d2d989aa9`
and plan SHA `2d4bfe0dc8114a8577edae98dcc96a377b08ba711027c94b5d196ee5c617fa54`.
This report establishes endpoint availability with explicit limits, not a
preference result, force qualification or production-default decision.
