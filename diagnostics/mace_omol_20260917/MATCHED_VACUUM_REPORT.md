# Solvent mismatch explains most of the masked hybrid partition failure

**Eight matched vacuum DFT endpoints completed.** Removing only CPCM reduces
the raw-zero masked hybrid's GGR partition discrepancy from20.1603 to
**1.8298 kcal-scale**, within the predeclared2 diagnostic tolerance. The solvent
contribution is18.3305kcal/mol. This is a concrete explanation for most of that
specific failure, not evidence of improved biological classification.

| Connected minus extended GGR contrast | kcal-scale |
|---|---:|
| Native DFT/CPCM | -7.343500873 |
| Native DFT/vacuum | -25.674002864 |
| Raw-zero masked MACE core | -27.503809106 |
| Shared learned-neutral MACE core | -18.255784798 |
| CPCM minus vacuum contribution | **+18.330501991** |
| Vacuum DFT minus raw-zero MACE | **+1.829806242** |
| Vacuum DFT minus shared-neutral MACE | **-7.418218065** |

Both component identities close exactly at reported precision. Unlike subtracting
the printed CPCM term from an already polarized SCF energy, these are separate
self-consistent vacuum calculations on identical nuclei. The difference includes
the density's response to solvent. It cannot be assigned uniquely to cavity
shape, reaction-field energy or density polarization without further decomposition.

The raw-zero model passes this vacuum partition diagnostic; the shared-neutral
model still fails. Both earlier CPCM hybrid candidates remain rejected, with
their conditional whole calls unexecuted. A separate vacuum hybrid requires a
new protocol and declared test. Vacuum consistency alone supplies no aqueous
affinity, calibrated zero, threshold, mechanical response or broad accuracy.
The r2SCAN-3c versus OMOL training-functional mismatch also remains.

## Gradients and scientific state

All eight analytic gradients pass actual output, element/order, coordinate,
energy, native D4/gCP and La46-electron ECP checks. All original source heavy/H
coordinates, cap mappings, charges(Ca−1/La0), singlet states and waters are
unchanged. The new parser explicitly labels `isolated_vacuum_endpoint`;
old CPCM extraction and the production baseline remain unchanged.

The result retains each endpoint's vacuum/CPCM gradient and projections onto
the exact physical metal translation and peptide rotation. For example,
GGR_extended La metal-direction gradients are−6.1478(CPCM),+5.6627(vacuum)
and+5.3820(masked MACE), in their declared kcal/Angstrom scales. This shows why
the previous cross-Hamiltonian gradient discrepancy cannot all be attributed
to MACE force error. Other projection errors remain; no curvature, relaxation,
uncertainty covariance or force-accuracy validation follows from eight centers.

## Execution and measured cost

Protocol `native_r2scan3c_matched_vacuum_response_diagnostic_v1`; ORCA6.1.1 native
r2SCAN-3c/DefGrid3/TightSCF/EnGrad. Eight new DFT calls, no failures, retries,
new MACE/solver/training calls, optimization or numerical DFT derivatives.
Bounded archive audit found no matching vacuum input among165 inputs in five
named project archives;25 shared the exact required coordinates.

Job1200905 completed on64CPUs with four concurrent16-rank endpoints:
**618 seconds wall,39,552 allocated core-seconds,35,369 reported CPU-seconds**.
No GPU allocation. This is a one-time diagnostic cost; no matched-hardware
production speedup is claimed. Per-endpoint and matched CPCM source receipts
are retained in `matched_vacuum_cost_v1.json`. Local preparation/report/test
resources are recorded separately.

Four distinct real-fixture tests pass: exact solvent-only input transformation,
rejection of real CPCM output as vacuum, corrupted actual numerical-gradient
input rejection, and actual executed-output/partition closure. The last was
explicitly skipped before outputs and then passed in1.812seconds.

## Artifacts and next step

Under `workspaces/mace_omol_20260917/`: `matched_vacuum_v1/manifest.json`
SHA256 `f43b66d7dd708082e643ef42f9ed939424787afabdcbd85f319cf7b237a98960`,
its eight actual receipts/outputs, `matched_vacuum_report_v1/result.json`,
`matched_vacuum_cost_v1.json` and the bounded archive audit. The complete result
has unrounded endpoint/gradient components; [compact result](MATCHED_VACUUM_RESULT.json)
and [commands](MATCHED_VACUUM_COMMANDS.md) support replay.

Recommendation: **pursue a consistently defined vacuum hybrid as a separate
development candidate**, retaining the production baseline. Its usefulness and
eventual solvent treatment still need testing; no default promotion is justified.
