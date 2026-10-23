# Six-parameter paired water representability test

Implement the root declaration in this directory as a fresh finite development run. The motivation is the completed17593 local four-parameter stop and17596's independent, parameter-linear epsilon responses. Its projected improvement is not an evaluated fitted candidate. Preserve both results, all ledgers and earlier protocols.

## Data, family and initial point

Use the exact eight consumed water cells and unchanged training result `3aefeb74`, radial objective `5b83cf5f`, rigid-body helper `bf65a20a` and native Tinker programs. The four equal loss blocks remain energy/La/Dy/paired. Energy retains its seven weighted residuals and nine controls; every individual and paired fixture retains RMS≤2 and maximum body component≤5 kcal/mol/Å. Energy controls remain≤1 kcal/mol. Keep original four-center diagnostics and reference-method sensitivity separate.

Parameter order is `[La_R,La_damping,Dy_R,Dy_damping,La_epsilon,Dy_epsilon]`. Start at exact17593 selected radius/damping `[3.832755583060068,0.31164088100092946,3.525530001514476,0.32135326991694113]`, followed by original epsilons `[0.24429636467571708,0.37332330078584813]` kcal/mol. The original four bounds are unchanged. Each epsilon has the positive interval0.5–1.5 times its original value; these are exploration limits, not uncertainty estimates.

Only explicit metal–water-O `VDWPR` radius/epsilon and `POLPAIR` damping may change. All geometry, reference states/settings, charges, bonded terms, permanent multipoles, polarizabilities, water-H and non-water pairs remain fixed. No donor parameter, solvent correction or target changes. Selected17593 native E/F and scores must replay exactly with zero new calls before release.

## Unchanged algorithm in six dimensions

Normalize each parameter to its frozen interval. At each origin evaluate all12 clipped central±0.005 parameter probes, recording actual widths. Build the weighted-residual Jacobian in the unchanged241-entry order. Use the existing least-squares cutoff1e−6, uniform-scaled Gauss–Newton direction and bounded least-squares direction under the same componentwise±0.025 trust box and original bounds. Generalize only the dimension: the bounded solver enumerates3^6 active-set patterns instead of3^4, preserving pattern order, least-squares rule and strict minimum selection.

Evaluate all four proposals in fixed order: scaled-GN half/full, bounded-LS half/full. Finish the complete derivative batch before directions and the complete trial batch before selection. The current point wins exact ties; select only an actual decrease of the unchanged objective. Derivative probes never replace the candidate. Stop when the selected point passes every original control, a complete trial batch produces no improvement,24 outer steps finish, a finite ceiling is reached, or native/admission failure occurs. No fallback, threshold/bound change or automatic extension.

## Finite calls, cache and resources

At most24 steps,384 new full parameter vectors (385 including the initial vector),4096 new native program invocations. A complete step has at most12 derivative and4 trial vectors; active-metal reuse gives at most96 derivative plus64 trial native programs, hence3840 for24 complete steps before reuse.4096 remains the independently enforced conservative hard ceiling.

Full-vector request identity includes all six values. A native cell's effective identity includes only that cell's metal radius, damping and epsilon, exact rendered input and scientific/engine provenance. Reuse exact admitted17593/17589 and diagnostic17596 physical cells where rendered input identity matches, without altering their source identities or old ledgers. The opposite metal is never recalculated solely because another metal's parameters changed. Use a fresh cache and ledger; preserve full eight-cell denominators and raw authorities.

One eight-CPU standard_expand task, up to eight workers, one native process/thread per worker, default memory/time. Stop new starts immediately on native/admission failure, retain completed/running siblings and their receipts. Persist full result at batch boundaries and failure; per-cell receipts/cache remain durable. Record actual program counts, wall time and aggregate parent/child CPU. Root owns release/submission.

## Admission and interpretation

Preparation and bounded saved controls verify exact initial replay, epsilon bounds/pair-only rendering, opposite-metal reuse, six-dimensional solver/clipping and hard limits with zero native calls. Root release binds the exact source/config/preparation and independent review.

This fits one consistently uncorrected PBE0 E/F development surface; CC/CP sensitivity remains descriptive and is not mixed into force targets. Passing consumed fixtures is neither physical validation nor a global optimum. Report bound contacts, rank/conditioning, residuals and actual stop. No reserved24 outcomes, donor fit/transfer, protein dynamics or production promotion is authorized here.
