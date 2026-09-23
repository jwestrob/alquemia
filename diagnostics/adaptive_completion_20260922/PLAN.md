# Complete the adaptive family under one numerical policy

Jacob approved completing adaptive accommodation under its own canonical-only
calibration and then structural transfer. Parent delegated this numerical
completion on 2026-09-22. The original four-angle search produced 59/60 valid
candidates. Its MMOL1770 La origin is valid; an oversized first SLSQP trial
created an atom clash before any accepted iteration. Preserve that failure.

## Fixed numerical version

Protocol: `common_four_angular_native_OMOL_SLSQP_hartree_units_v2`.

Use the same physical four-angular-mode selector at the same original paired
q0, native MACE energy/analytic forces, charge/protonation/water inventory,
original source maps, 0.8 Å final source-heavy displacement constraint, ±0.8 rad
component bounds, one original start per endpoint, and maxiter200. Add no metal,
water, cofactor or scaffold degree of freedom. Keep every chemistry/overlap guard.

The **only optimizer numerical change** is a constant positive unit scaling:

`f_H(q) = [E_MACE(q) − E_MACE(0)]_eV × EV_TO_KCAL / HA_TO_KCAL`.

Multiply the analytic objective gradient and the former 1e−9 eV ftol by the
same fixed factor. Existing project conversion factors define this factor;
it is not selected using labels or score outputs. Store and report actual MACE
energies in their original eV units and raw gradients in their existing units.
No energy or gradient is evaluated at a substituted/clipped geometry.

Installed SciPy uses an identity initial BFGS matrix and requests objective
values before updating nonlinear constraints. Thus objective scaling changes
its initial proposed displacement and numerical path, while preserving the
physical objective, its stationary points and feasible set. SciPy's single ftol
also enters constraint/step tests; scaling it makes those internal checks
stricter, not exactly invariant. The final 0.8 Å acceptance/tolerances remain
unchanged. SLSQP may still request an invalid trial; retain an explicit failure
if so. This version adds no fallback, restart or alternate optimizer.

## Declared sequence and selection

1. Geometry-only first-step probe on all 30 original sources and both metals,
   using the exact archived q0 gradient. Stop before returning any nonzero-point
   objective value. This inspects numerical proposals, not scientific energies.
2. Eight-endpoint pilot: 1H4I, 4MAE, `mmol_1770-pqq-la_model`, and
   `q9z4j7-pqq-la_model`, both metals. Parent receives preflight/manifest before
   coordinated launch. The same fixed unit conversion applies to every point.
3. Apply this numerical version to all remaining 26 original sources (52
   endpoints) on the same declared30/60-task manifest. The pilot comprises four
   sources/eight endpoints, not eight sources. Existing pilot completions are
   reused when the full60 phase runs. Do not merge old-optimizer candidate
   labels into v2 completion. The parent's new five-candidate pool contains
   q0, the earlier terminal-only Ca/La proposals, and the two scaled-angular
   proposals. Previous unscaled four-angle scores remain separate; do not add
   their geometries only for the failed case. Retain all failures/denominators.
4. Prepare primary225 structural transfer with the same selector/numerical rule,
   retaining unavailable inputs and exact source states. Root owns pool scoring,
   canonical-only calibration and frozen-reference transfer timing. Crystals,
   PLM predictions and noncanonical folds do not set the calibration.

No scientific calls are used to choose another scale. If the single declared
policy remains inadequate, report that result before proposing another version.
No flexible fit or label-selected candidate/parameter is permitted.

## Execution and interpretation

Existing warm native-MACE worker, H200, 32 CPUs, 200000 MiB host memory, one GPU.
Finite manifested endpoint searches; no artificial aggregate budget or time
limit. Existing exact evaluations may be reused only with a traced matching
state/model/geometry receipt. No DFT or GFN2 call in this adapter, no new folds,
ensemble or chemical states. Root owns separate solvent manifests and scores.

Numerical completion is necessary for a fair adaptive calibration and transfer;
it is not evidence of improved classification. Preserve baseline, production,
original failed records and the separate joint-metal branch.
