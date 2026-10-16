# Matched seven-versus-ten coordinate continuation

Declared16October2026 after range searches12544/12545 and before any new
basis/control molecular call. Range pool12546 is running; this rule does not
select motions by its result. Protocol: lbt3_common_basis_response_v1.
Root owns submission/integration. This remains consumed peptide development.

## Scientific question

Do coupled donor/backbone freedoms reduce source dependence beyond further
search in the existing seven coordinates? Saved diagnostic12542, predating the
range experiment, records omitted own-metal loads on both sources/metals.
Add exactly existing mode9 A3:sidechain:CB-CG, mode38 A11:sidechain:CB-CG and
mode32 A10:phi:N-CA to the old seven. They move complete mapped donor groups
or the complete downstream peptide through Ala17/OXT; no detached fragments.
Use the same ten names/mappings for both metals and both sources, retain all
256 atoms in the energy, and never optimize the La/Dy energy difference.

Compare two arms: old seven-coordinate control and ten-coordinate candidate.
Both start at the SAME exact range terminal for each source/metal. Reset
optimizer history. Preserve source displacement .40Å, target .20Å, all torsions
±.20rad relative to the original source, all-atom trial .05Å, all existing
covalent/contact/stereochemical checks, original atom/state/water/proton/model
and solvent definitions, numerical energy/force gates and same maximum12 new
accepted steps/24 new distinct composites per search. The small budget is a
finite experimental definition, not a total project resource cap.

## Coordinates and geometry

The absolute64-component q vector remains authoritative. New mode entries start
at zero. Derive each arm's transform from its ORIGINAL-source Jacobian using
physical-displacement normalization and the existing relative1e-10 rank cutoff.
Each metal shares its source/arm transform. Write updates as
q_active = q_start_active + T_arm @ v, with v initially zero. Never reinterpret
old seven-component scaled coordinates as ten-component coordinates or move
Cartesian atoms when changing metric. Reuse exact starting coordinates,
energies and full gradients; retain original source as every physical bound
and contact reference. Recompute moving/fixed membership from the actual active
union. Use the full existing64-mode kinematic order, not a manually appended
backbone rotation. Add zero-call actual-fixture centered finite-difference
Jacobian checks for the three added modes at original and range geometries:
steps1e-5 and5e-6radian, maxabsolute Cartesian derivative residual≤1e-6Å/rad.
This checks mapping/Jacobian, not the accuracy of the molecular Hamiltonian.

## Feasible motion without relaxing bounds

The current unconstrained direction may point outside a limit even when other
downhill torsions remain. Apply the SAME deterministic feasible fallback to
both arms. First attempt the existing four Armijo multipliers on the L-BFGS
direction. Only if all four trials fail geometry BEFORE molecular evaluation,
consider one-coordinate own-gradient descent directions in the physical q basis.
For each active coordinate use the downhill sign, maximum tangent atom step
.05Å and the same four multipliers. Replay exact geometry/guards first; keep the
largest admissible multiplier for that coordinate. Rank admitted candidates by
predicted own-energy decrease from the current analytic gradient; stable tie
order is ascending numeric active-mode index. Try at most the top two admitted directions, using
the same Armijo rule and remaining24-call allowance. Reset L-BFGS history on an
accepted fallback step. No entry into fallback after any energy-evaluated trial on the initial
L-BFGS branch, including a cached evaluation. Once fallback has been entered
after four geometry-only rejections, the second ranked coordinate may be
tried if the first fails Armijo. A backend failure stops the search rather
than triggering another direction. No recursive fallback, random starts,
extra displacement or unbounded geometric search. This interpretation is
clarified before the first basis/control evaluation.

A blocked fallback is a finite-search outcome, not constrained stationarity.
Record failed guards, proposed/accepted fallbacks, residual forces, steps/calls
and original/incremental work. Backend failures leave the new arm unavailable;
valid old states remain separately named. No score-sign rescue or hidden fallback.

## Execution and paired comparison

Internal zero-call preflight covers both arms/allfour states first. Scout is
7CCO/La control then candidate with one warm native model and exact cache reuse.
Execution/physical admission permits the unchanged remaining six searches;
no gain/sign threshold. A new arm's failure need not stop independent states,
but a complete arm comparison requires its whole declared denominator.

Keep prior ten-candidate/twenty-cell range pool (or actual deduplicated counts).
Add each of eight new admitted terminal points, deduplicate exact coordinates,
and score both metals on every unique candidate. At most eight new opposite-state
composites after eight searches: planned maximum200MACE/400GFN calls, not implied
execution. Reuse all compatible existing physical points and cells exactly.

Report prior/control/candidate per-source common menus and global menus separately,
with paired source spread, individual work, selected geometry, remaining bounds,
coverage, actual costs and optional full union. Do not select whichever arm or
subset gives a desired label or call a larger common menu equilibrium sampling.
This is not a calibrated La/Dy classifier, binding free energy or natural LanM
validation. No extra bond flexibility, larger displacement, quantum campaign,
dimer/library rescore or PQQ change follows from this declaration.
