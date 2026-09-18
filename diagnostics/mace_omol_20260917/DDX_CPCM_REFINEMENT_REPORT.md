# Conductor contrasts stabilize; endpoint qualification still fails

The higher-resolution inventory completed all sixteen new native solves and
reused four matching prior endpoints. **444 of 450 checks pass.** All paired
contrast and between-structure refinement checks pass; four individual endpoint
refinement checks and two rigid-transform endpoint checks fail. The model is
not numerically qualified for a full score.

| Resolution | 2FW0 R | 2FVY R | 2FW0 minus 2FVY |
|---|---:|---:|---:|
| 12/590, reused | -44.2598926681 | -46.1517241429 | 1.8918314748 |
| 18/974 | -44.3826110003 | -46.3084560052 | 1.9258450049 |
| 24/2030 | -44.4608903299 | -46.4055241837 | 1.9446338538 |

R is the source-self Ca-minus-La component in kcal/mol, not a complete affinity
score. The between-structure contrast changes 0.01878885 kcal/mol between the
two new resolutions. Its large difference from archived GK's 18.12902095 kcal
persists. This shows model dependence; it does not prove which model predicts
binding better.

Endpoint refinement changes are -0.416525/-0.338245 kcal for 2FW0 Ca/La and
-0.362134/-0.265065 for 2FVY. Rotation changes fail for 2FW0 Ca (+0.055733) and
2FVY La (+0.057281), against the unchanged 0.05 kcal gate. All reciprocity,
zero/repeat, source, charge, scaling and receipt checks pass. No tolerance was
relaxed; both earlier frozen conductor inventories remain failed.

## Actual cost

Job 1201302: 1,446 allocation-wall seconds, **92,544 allocated core-seconds**,
81,415 reported CPU-seconds (whole-second precision), 2,444,996 KiB Slurm-sampled
peak RSS; no GPU, DFT, MACE or forces. Primary endpoint solves take roughly
59–65 seconds, refined solves 162–165 seconds, plus setup and reporting.
This is development validation, including rigid/zero/repeat controls; it is
not the cost of one ordinary score. Read-only numerical report: 120.139343
seconds, no solver calls.

Manifest SHA256:
`b6225d23440423c917009009c5e7d6581445fd5207057e3668a40df60f45d300`.
Actual artifacts under `workspaces/mace_omol_20260917/`:
`ddx_cpcm_source_v2/collection_job_1201302.json`,
`ddx_cpcm_source_report_v2/result.json`, and `ddx_cpcm_cost_v2/`.

Next: DDX_RESOLUTION_PLAN.md separates angular-basis and surface-integration
errors with eight solves on the same 2FW0 pair. It preserves the physical model
and thresholds and does not substitute for the full qualification inventory.
Baseline/default, labels, aquo reference and calibration remain unchanged.
