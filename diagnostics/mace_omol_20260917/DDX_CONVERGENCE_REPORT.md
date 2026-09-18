# The first La iteration failure is recoverable

The unchanged GGR 2FW0 La coarse calculation converged with the same 1e-10
iterate-change tolerance after raising only maxiter from300 to1200. The native
log records337 dielectric-system iterations and79 subsequent single-layer
iterations. Thus the original failure is slow convergence for this state,
not evidence of an unstable or impossible electrostatic solution.

| Actual result | Value |
|---|---:|
| Reaction component, Hartree | -0.06510015720125943 |
| Reaction component, kcal/mol | -40.85096540267961 |
| Direct source-potential check, au | 2.706168622523819e-16 |
| Monopole-integral normalization error | 1.7763568394002505e-15 |
| Independent energy contraction error, kcal/mol | 0 |
| Setup, seconds | 1.825818743556738 |
| Forward solve, seconds | 264.1825735978782 |
| Full driver, seconds | 267.609367005527 |
| Process CPU, seconds | 16831.655566545 |
| Process peak RSS, KiB | 376568 |

Job1201280 completed:280allocated wall seconds,64CPUs,17920allocated core-seconds,
16842reported CPU seconds (whole-second precision),330096KiB Slurm-sampled peak
RSS. The process peak and scheduler sample differ; both are retained. No GPU,
new DFT, MACE, geometry, charge-fit or force calculation. Preparation/local
inspection costs are separate and not fully measured.

The charge sum remains the printed source value,0.000003999999999840936e;
it was not renormalized. The result is a projected-charge source-self component
in a full physical protein cavity. It supplies no full hybrid score, affinity
classification, calibration or proof of numerical convergence with resolution.
The native stopping criterion measures iterate changes, not an independently
evaluated equation residual.

## Provenance and next execution

- Plan: DDX_CONVERGENCE_DIAGNOSTIC_PLAN.md, declared before execution.
- Manifest: `workspaces/mace_omol_20260917/ddx_convergence_v1/manifest.json`,
  SHA256`38a3368bed372eda419918ab6a610a8bcda1384809cb396a26fedaac3ea98020`.
- Actual receipt: same directory,`collection_job_1201280.json`.
- Native log and coefficients: `attempt_0001/solver.log` and`arrays.npz`.
- Original300-iteration failure: direct-source job1201279; preserved unchanged.

DDX_ITERATION_RECOVERY_PLAN.md declares independent fresh Models for missing
roles in the original numerical inventory and exact reuse of successful states.
The original job must finish before its final reusable inventory is prepared.
Fresh Models also prevent an existing native error flag from blocking later
endpoints before their solve starts. Baseline/default and scientific tolerances
remain unchanged. No production or predictive improvement is claimed yet.
