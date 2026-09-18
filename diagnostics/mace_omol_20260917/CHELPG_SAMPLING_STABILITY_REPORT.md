# Charge sampling matters, but most of the GK discrepancy persists

All 40 native fixed-orbital property replays and 40 saved-density queries
completed. All state, energy, density and potential-fit checks pass. The
predeclared 0.5 kcal source-self sensitivity screen fails 5 of 16 comparisons.
No SCF optimization, geometry change, MACE inference or baseline change.

| Source-self Ca-minus-La, kcal/mol | Default 0.3/2.8 | Finer 0.2/2.8 | Finest 0.15/2.8 | Extent 0.3/3.5 |
|---|---:|---:|---:|---:|
| GGR 1GLG extended | -16.0354 | -15.9744 | -16.0201 | -16.2542 |
| GGR 2FW0 | -14.0271 | -14.6936 | -14.6334 | -15.0361 |
| GGR 2FVY | -32.1561 | -31.0789 | -30.7874 | -30.0033 |
| Alpha 1F6S | -44.5125 | -44.2101 | -44.4145 | -44.0782 |
| Alpha 6IP9 | -36.9860 | -37.7568 | -37.4854 | -36.5436 |
| **GGR 2FW0 minus 2FVY** | **18.1290** | **16.3853** | **16.1540** | **14.9673** |

Headers give grid spacing/fitting extent in Angstrom. These are source-self GK
components, not full hybrid scores or affinity predictions. All eight finer-to-
finest checks pass; three of eight extent checks pass. Changing the fitting
region changes which exterior potentials constrain the atomic charges, even
though the electronic density stays fixed. No setting was selected per protein.

## Interpretation

Sampling changes the problematic GGR difference by up to 3.16174 kcal, but it
remains 14.9673–18.1290 across all four declared settings. Sampling does not
explain most of the discrepancy. The remainder could reflect physical effects
or representation/model errors; this experiment does not distinguish those.
Direct, induction and full hybrid terms were not recomputed.

Passing the ordinary potential-fit screen does **not** certify the solvent
energy to 0.5 kcal precision. Retain both checks when developing this model.
Increasing grid density alone would not address its fitting-extent dependence.

Maximum energy identity error: 3.23739e-8 Hartree; potential identity error:
1.03001e-13 au. Default charges reproduce archived values exactly at printed
precision; source-self energies reproduce the native-verified decomposition
within 1.599e-13 kcal. Largest sampled charge change: 0.049527 e. All charge,
ECP, source-projection and paired-coordinate checks pass.

## Execution and limits

Jobs 1201173/74 qualify the route; 1201177 runs the remaining 38 high-level
NoIter property evaluations; 1201211 runs the remaining 38 density queries.
Initial results are reused exactly. No scientific retry. Raw missing-SCF-
convergence flags remain separate from the explicit successful NoIter contract.
See [qualification details](CHELPG_SAMPLING_QUALIFICATION_REPORT.md).

Total cluster cost: **1030 summed job-wall seconds, 65920 allocated core-seconds,
53718.534 reported CPU seconds, zero GPU**. The largest job reports CPU only
to whole seconds. Peak reported RSS: 11791988 KiB. Local report replay took
18.373 seconds. Other local preparation timings are incomplete, so this is
not a complete end-to-end runtime or an ordinary production-score cost.
All 11 real-fixture tests pass in 9.732 seconds, with no final skip.

Full output: `workspaces/mace_omol_20260917/chelpg_sampling_report_v1/result.json`.
Compact pins/values: `CHELPG_SAMPLING_STABILITY_RESULT.json`.
Cluster receipts: `chelpg_sampling_cost_v1/` under the same workspace root.
Five geometries remain two consumed biological groups. No blind validation,
new threshold or demonstrated predictive improvement follows.

Next, test the reaction-field approximation on the same physical cavities and
fixed default sources with a maintained PCM solver. The isolated ddX build
imports successfully; its scientific convergence and cost remain untested.
Exact-density coupling requires the additional integral documented in
[the interface assessment](DDX_CAPABILITY_NOTE.md).
