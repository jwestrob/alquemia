# Full polarization accounting passes the native operator check

Four real GGR Ca/La endpoints pass all **28 numerical checks** after a wrapper
initialization repair. The stationary bilinear functional reproduces the
archived induction energy within 4.60e-8 kcal/mol. Two independent evaluations
of the full GK reaction term agree within 2.73e-11 kcal/mol. Independently
evaluated response residuals are below 3.00e-10 Debye RMS.

This validates the implemented frozen-response energy accounting against the
existing native GK model. It does not establish a conductor correction,
self-consistent new response, improved discrimination or compatible calibration.

| Actual state | Full GK cross reaction term, kcal/mol |
|---|---:|
| 2FW0 Ca | -2441.7874393427032 |
| 2FW0 La | -2381.751218575202 |
| 2FVY Ca | -2231.8411100320427 |
| 2FVY La | -2160.185624092197 |

These include the permanent environment and both induced-dipole representations;
they are not site scores and cannot be compared directly as affinities.

## What ran

The corrected job 1201310 performed four operator replays: four Born-radius
evaluations, four permanent-field evaluations, four mutual-field evaluations,
twelve static GK electrostatic energy evaluations, and four nonpolar evaluations.
No induced-dipole iterations, DFT, MACE, continuum solves, forces or new labels.
The pinned native Tinker library remained unchanged. Actual inputs, driving
fields and archived dipoles were read at full available precision.

The first job 1201309 had a missing native cutoff initialization in the new
wrapper. Its failed checks and outputs remain intact. The repair restores
the same `switch('MPOLE')` call used by native `induce0c`; the corrected wrapper
also prints and asserts the cutoff. No physical model or tolerance changed.
See FROZEN_RESPONSE_RECOVERY.md.

Each job used 29 allocation-wall seconds on 64 CPUs. Including the failed
attempt: **58 wall seconds, 3,712 allocated core-seconds, 234.977 reported
CPU-seconds**, zero GPU. Local isolated builds took 1.164453 and 1.165101 seconds
(1.116202 and 1.125109 CPU seconds), separately from cluster cost. Slurm sampled
no useful main-process RSS for these short runs; process resource files remain
available. Three real-input/output tests pass in 22.090 seconds, no skip.

## Reproducible artifacts

Protocol: `native_AMOEBA_GK_frozen_response_functional_check_v1`.
Corrected manifest SHA256:
`ebd1cb961adb52bb99392a6006065b9f5e2e7db5d117fc7a2423e1347d0dcec2`.
Under `workspaces/mace_omol_20260917/`:

- `frozen_response_v2/collection_job_1201310.json`: passing native outputs.
- `frozen_response_v1/collection_job_1201309.json`: retained wrapper failure.
- `frozen_response_software_v{1,2}/`: separate source/build receipts.
- `frozen_response_cost_v1/`: both allocation receipts and summed cost.

The full derivation and unchanged acceptance rules are in
FROZEN_RESPONSE_TRANSFER_PLAN.md. Next qualify a matched conductor representation
of these same permanent and induced multipoles. No full correction is available
yet; baseline/default and references remain unchanged.
