# Two-seed recovery of one failed native scalar cell

Approved research, September 23, 2026. Root assigned this bounded experiment
under Jacob's current discretionary overnight authorization: “alright. keep on
goin buddy, launch subagents to do what experiments you like so they ping you on
completion and keep you awake. i'm goin to bed. email me when cool stuff happens.”
Root's exact instruction: prepare exactly two seeded attempts at the failed exact
coordinates; origin-seeded value is proposed recovery, adaptive_La-seeded value
is the agreement check. Both must converge and agree within 0.1 kcal/mol.

## Fixed question and scope

Can a different initial electronic state resolve the actual 500-cycle oscillation
in A0ACD6B9F2 La-conditioned sample0, La at adaptive_Ca, vacuum, under its unchanged
strict native scalar Hamiltonian? The original full100 result remains unavailable;
this is a separately named recovery sensitivity, never a primary replacement.

- Exact case: `a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-0__envelope_2eda778295ffed58`.
- Target is the exact archived adaptive_Ca coordinates, 204 atoms, charge −3,
  singlet, 636 native electrons, 573 basis functions/325 shells.
- Exactly two starts, ordered `origin` then `adaptive_La`; both seed files
  (matching GBW and native xtbw) must come from the actual successful vacuum
  calculation of this same source/context/metal/state and atom ordering.
- The documented matching-runtime-basename AutoStart route activates native
  XTBRESTART. Require that actual printed marker, native mixer and parameter/state
  audit, normal termination, SCF convergence, rank1 and printed TolE1e-10.
- Same ORCA6.1.1 nativeGFN2, TolE1e-10 Hartree, MaxIter500, 300 K, vacuum,
  maxcore2000 MiB, one rank/cell. Only NoAutostart is removed to permit the seeds.
- Two independent concurrent workers on 2 CPUs/8 GiB, CPU-only GPU-partition host;
  no GPU, new MACE, DFT, geometry, protonation or other molecular work.
- Both must converge and their energy difference must be ≤0.1 kcal/mol. Always
  retain both. Use origin-seeded energy only if the gate passes, regardless of
  which energy is lower. Otherwise recovery is unavailable; no third attempt.
- Preserve raw original failure, all seed/input/output/receipt pins, both attempt
  costs and failure reasons. No shared collector, candidate or reference edits.

Strict32 prior fresh/seed agreement supports scalar reference compatibility;
it does not itself establish this failed cell's recovery, force accuracy, or a
new classification. Root may separately join a qualified recovery sensitivity.
