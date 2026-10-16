# Uniform displacement-envelope continuation

Declared16October2026 after complete scout12533, continuation12538, paired pool12541
and saved diagnostic12542, before any extended-range molecular evaluation.
Protocol: lbt3_whole_composite_range04_continuation_v1. Root owns execution.

## Question and reason

Does the original0.20Å maximum atom displacement artificially truncate useful
response in the existing seven-dimensional motion set? All four previous searches
stopped when the next four allowed trial steps exceeded this envelope at Asp7
OD1/OD2; those trials passed the separate covalent/contact checks. Active forces
remain substantial, particularly7CCN/La. The shared-pool result reduced source-local
contrast spread72.966739→57.895436kcal/mol but did not establish robust preference.

Change ONLY maximum displacement of any physical atom from its original source:
0.20→0.40Å. Retain target displacement≤0.20Å, each torsion±0.20rad relative to
its original source, trial maximum atom step0.05Å, all covalent/protected-group/
chirality/peptide-basin/fixed-atom/severe-contact checks, seven original common
coordinates, numerical tolerances and complete energy/analytic forces. Do not
reinterpret an inherited contact baseline at the last accepted geometry.

## Four-state continuation and reuse

Use exactly7CCO/La,7CCO/Dy,7CCN/La,7CCN/Dy in that declared order. Every search
starts at its admitted old terminal point with its exact energy/gradient and
original-source coordinate vector. Keep its ORIGINAL prepared source as the
geometric comparison baseline and work-reference; report incremental work from
the continuation start separately. The initial source-metric whitening remains
unchanged. Reset L-BFGS history; no changed geometry, smearing, state, proton
inventory, water inventory, model parameters or solvent recipe.

The underlying molecular evaluator/energy-cache identity stays physically
compatible with the prior calculation. The changed search protocol and envelope
receive this new identity; no old scores/results are relabeled. Reuse compatible
saved calls exactly, including a trial that duplicates a prior physical point.
No previously successful molecular evaluation is repeated just to initialize.

Each continuation uses at most12 additional accepted steps and24 new distinct
composite evaluations, the existing four Armijo multipliers and coefficient,
existing gradient stopping criterion and one warm model across states. This is
a finite experiment, not a total project compute cap. New limiting target/torsion
bounds, active forces and no-admissible-trial outcomes are results. No automatic
second bound expansion, internal bond flexibility or score-sign rescue.

Run7CCO/La as the first actual continuation scout. If execution and physical
admission pass, continue the other three unchanged irrespective of improvement
size or preference direction. Backend/state/gradient failures leave that new
search unavailable and preserve earlier valid points separately; no hidden
successful fallback. Distinct unaffected states may still be evaluated.

## Paired comparison

Preserve the complete prior six-candidate/twelve-cell pool. Add each newly
admitted terminal point with exact coordinate deduplication and cross-score
both metals on every unique physical geometry. At most four new opposite-metal
composite cells are needed after the own searches; together the four searches
and cross-evaluation allow at most100 new MACE and200 new GFN2 calls. Planned
maxima are not implied executions. All numerical work uses shared Slurm and
actual completion wakes. No quantum campaign follows from this declaration.

Report old versus expanded per-source common menus and the global common menu,
individual work, selected source/geometry, residual forces, all active bounds,
coverage and actual calls/resources. A global shared menu removes source choice
by construction; only source-local comparisons diagnose remaining source bias.
The experimental two-source LBT3 panel is consumed development evidence; there
is no affinity, population, classifier threshold, Mex-site label or equilibrium
claim. The natural-protein LanM goal remains active.

Omitted sidechain/backbone loads justify a later independent common-subspace
comparison, but those freedoms are not added here. Large carboxylate internal
forces are not proof of defective deposited geometry; bond-stretch freedom needs
separate physical qualification. Preserve PQQ, old failures and all source states.
