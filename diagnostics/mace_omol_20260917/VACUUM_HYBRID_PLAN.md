# Separate vacuum hybrid: six whole-protein calls

Declared after the matched-vacuum diagnostic and before these new inference
outputs. Active-goal authorization applies. Both earlier CPCM candidates remain
failed. This is a separate protocol, not execution of their conditional stages.

Protocol: `masked_omol_subtractive_context_DFT_vacuum_v1`.

## Question and expression

Does a source-matched whole-protein context improve the consumed alpha/GGR
ordering when both sides of the core replacement use vacuum energies?

    H_M = E_DFT,vacuum(core,M) + T_mask(full,M) - T_mask(core,M)
    R_H = H_Ca - H_La.

Use unmodified native r2SCAN-3c vacuum endpoints from matched_vacuum_v1 and the
same pinned raw-zero masked OMOL descriptor. Fixed unit coefficients, Hartree
and eV each converted once. Retain baseline/CPCM, vacuum core, learned core,
learned full, context and final contrast side by side. This is an empirical
vacuum hybrid descriptor: no continuum, aqueous reference, physical zero,
calibrated decision, entropy, mechanical correction or binding-free-energy claim.
No additional electrostatic term is silently set to zero or double counted.

The already declared diagnostic found a1.829806242kcal-scale GGR partition
residual, passing its fixed2 limit. That result motivates this experiment but
is not an independent validation case. No charge-feature/category changes.

## Exact input inventory and reuse

Use the three `physical_cases/.../grids/center` full systems in the exact
mechanics preparedV2 manifest, SHA256
`361d2ac85b93360c8560a761e0d76d622b8f314a7c238013f8fb7faa184ae99c`:
GGR_1GLG, ALPHA_1F6S and ALPHA_6IP9, each bound Ca/La. Six new energy-only MACE
calls, using the existing qualified float64 raw-zero feature and exact
1024-edge/product adapter. No new DFT, solvent solver, gradient, training or
optimization calls. Reuse all eight actual raw-zero core-center calls and eight
newly completed vacuum DFT centers. The two GGR cores share exactly one full
protein; their full contribution must cancel in the partition contrast.

Preserve original source H and heavy coordinates, charge/spin, assembly,
protonation and water inventory. The old context report's source/cap mapping
checks apply to these exact unchanged geometries. Verify their actual pins;
do not substitute the differently H-normalized direct-score benchmark. Known
original H-bond-length defects remain an explicit limitation, not a new silent
repair. Any later uniform H normalization needs new matched DFT and a new version.

## Fixed reporting gates

Require native-readout/component closure <=0.01model kcal and independent
endpoint-versus-component algebra <=1e-7kcal-scale. Report the same<=2kcal
partition diagnostic. Missing/incompatible endpoints block the paired result;
no baseline substitution or missing correction filled with zero.

The consumed development ordering gate requires all four alpha-minus-GGR
contrasts to exceed0.02kcal-scale: both alpha structures versus each GGR core.
Report every contrast, denominator and the learned full term. GGR is one direct
Ca-favoring observation; both alpha structures are one condition-qualified
La-favoring observation from another study. This is not a common-assay
cross-protein affinity scale or broad validation. No old PQQ band transfers.

A pass permits planning further real-structure and uniformly prepared tests;
it does not promote this model. A failure is retained without fitting a weight,
selecting a favorable GGR core, changing a label or adding solvent to rescue
these outputs. Any later solvent-inclusive expression is a separate model.

## Execution

Existing manifested MACE runner and oneA5000/16CPU/64474MiB allocation. Comparable
GGR forwards cost about13seconds; expect a few minutes including loading and
validation. Record all actual preparation, inference, failure and allocated
costs. No project CPU/time limit; scheduler QOS and concurrency policies remain.
Baseline, source experiments, references, PLM jobs and pending native H200 work
remain unchanged. Actual scientific integration tests must wait for real outputs.
