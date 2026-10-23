# Bounded joint refinement of the same four water parameters

Start from the actual17575 selected `bounded_LS_full` point, scout RESULT `dfb777e38e152f7333960ac0989243c2c911af677903a16073e87ffb2546b13c`. Its successful coupled improvement establishes that the earlier coordinate stop did not exhaust this local response. La's maximum-component control still fails. This single finite refinement preserves the four parameters, original bounds/epsilons, four water geometries, PBE0 energy/force targets, symmetric weights, separate controls and method/state limitations. No donor fit, validation opening or production release.

## Exact iteration

Allow at most12 outer steps. At each current point:

1. Stop immediately if all original selected-point development controls already pass.
2. Evaluate the four parameter derivatives at ±0.005 of their original normalized interval widths. Clip to the original bounds, reuse completed identical physical inputs and divide by the actual normalized endpoint separation. Complete all eight endpoint slots before building the weighted-residual Jacobian. A clipped endpoint equal to the current point is reused.
3. Solve the same weighted linear least-squares problem with SVD relative cutoff1e−6. Form the uniformly scaled Gauss–Newton direction and direct box-constrained least-squares direction, each restricted to max|Δz|≤0.025 and the original parameter box. The latter uses the previously reviewed81-pattern active-set method. A zero direction is a reused point, not a fallback or expanded trust region.
4. Evaluate exactly the half/full points of scaledGN, then half/full points of boundedLS. Complete all four slots before selection; preserve fixed order independent of completion order. Select the lowest **unchanged** objective with the current point first and strict improvement required. Keep current on ties. Report every point's separate controls, including any passing but unselected trial.
5. Stop if the selected point passes all existing controls; otherwise stop on no actual improvement after the complete trial batch, or after12 completed outer steps. No automatic continuation, smaller derivative step, larger trust region, different optimizer, epsilon fit or threshold change.

Maximum144 distinct derivative/trial parameter requests beyond the initial point and1152 new native program invocations (12×12×4×2). Reuse can reduce both. The two inherited baseline admissions are excluded from this new-vector count and make zero native calls. Runtime independently enforces the new bounds in addition to the unchanged evaluator's241-vector/1920-invocation ceilings.

## Cache and implementation

Use a fresh workspace/cache/zero-count ledger. Reuse the frozen scout's saved-admission helper to admit the original pilot, then add read-only-use links to all16 completed scout native-cell directories. Re-admit the exact selected scout model/score before any new native call. Preserve the original pilot and scout results/ledgers and verify their hashes after execution. Do not continue either old ledger.

The evaluator, native parser, rendering, polarization settings, body response and objective sources remain byte-identical. Only the new runner selects the declared finite joint steps and adds per-program1152 and new-vector144 guards. One root-owned eight-CPU allocation runs up to eight independent one-thread derivative workers; four trial workers use the same allocation. First native/admission failure stops new scheduling, preserves all iteration/slot records and yields no retry. Costs use aggregate start/end child CPU plus parent CPU and scheduler accounting.

Saved-only preparation pins source/config, current four-cell baseline, cache provenance and limits before root release. It makes zero native calls. Numerical controls-pass stopping is distinct from a global minimum, independently validated parameters or physical accuracy. Final method sensitivity is reported against the same saved CC/CP derivatives; no alternative reference surface enters the loss.
