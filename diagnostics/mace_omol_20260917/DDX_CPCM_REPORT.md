# Affordable conductor solves reveal strong model dependence, but need refinement

The full20-role conductor-like inventory completed in285allocated wall seconds
on64CPUs, without GPU, DFT or MACE calls. All20native solves succeeded.
The numerical report passes420/430checks and **fails its frozen overall gate**.
This is a source-self reaction component, not a new discriminator score.

| Resolution | 2FW0 Ca-minus-La | 2FVY Ca-minus-La | 2FW0 minus2FVY |
|---|---:|---:|---:|
| 6/194 | -44.2737617941 | -46.1807409929 | 1.9069791988 |
| 9/302 | -44.0757818952 | -46.1739813174 | 2.0981994222 |
| 12/590 | -44.2598926681 | -46.1517241429 | 1.8918314748 |

All values are kcal/mol. The archived GK between-structure source-self contrast
is18.1290209498kcal. Coordinates, source charges and intrinsic atom radii are
held fixed. The GK neck/descreen approximation and conductor sphere boundary
are different physical models; this does not identify one isolated cause or
prove that the conductor answer is correct. It does show that the large GK
source-self discrepancy is strongly model dependent, beyond the earlier
charge-sampling effect. No classification was used to choose this model.

## What fails and what passes

The nine failures beyond one coarse reciprocity failure are refinement/rotation
screens. Primary-to-refined endpoint changes reach1.05777kcal; the2FW0contrast
changes0.18411kcal and the between-structure contrast0.20637kcal, exceeding0.1.
The2FVYCa rotation change0.05496kcal narrowly exceeds0.05. Coarse2FW0reciprocity
is0.08840kcal versus0.05. All remaining reciprocity, zero/repeat, source,
energy-accounting, parameter, charge and paired-state checks pass.

The first inventory remains failed. DDX_CPCM_REFINEMENT_PLAN.md declares a new
18/974primary versus24/2030refinement, with all original thresholds unchanged.
It reuses four actual12/590endpoints and executes16new roles. No charge/radius/
water/geometry change, fitted threshold or favorable-case selection.

## Actual cost and provenance

Job1201299:285allocated wall seconds,18240allocated core-seconds,
13388reported CPU-seconds(whole-second precision),721624KiB Slurm-sampled peak
RSS. Summed group execution253.9456wall/13357.2776CPU seconds. Preparation,
reporting and local tests are separate: read-only report45.5377seconds; three
actual-source/output tests pass30.145seconds. A fourth refinement-reuse test
subsequently passes with the full four-test set in45.049seconds, no skip.

Primary endpoint solves take8.1–9.3seconds each;12/590solves take24.0–26.6seconds.
Add actual model setup/source preparation (roughly2seconds at9/302 and3.7–3.8
at12/590), runner validation and collection. These timings do not establish
ordinary production cost or full-hybrid predictive value.

Protocol`fixed_source_full_protein_ddCPCM_component_v1`. Energy is the native
conductor energy multiplied by(78.3-1)/78.3, then converted once from Hartree.
Raw energies, coefficients and scale are retained. No old aquo reference or
baseline bands are inherited. Baseline/default unchanged.

Actual collection:
`workspaces/mace_omol_20260917/ddx_cpcm_source_v1/collection_job_1201299.json`.
Manifest SHA256`e5ac44778fa284028c1895903f334d8d355533569128c314a7fbfd275ebf968d`.
Read-only numerical report:`ddx_cpcm_source_report_v1/result.json`.
Higher-resolution manifest:`ddx_cpcm_source_v2/manifest.json`,job1201302.

Recommendation: pursue this component's numerical refinement and coherent
full-energy accounting. It is computationally promising; numerical acceptance
and improved biological discrimination remain unestablished.
