# Live status — 2026-09-19

This is a progress record, not a completed accuracy claim.

## Context-supported geometry

- Ten native MACE searches completed in 45.05 worker seconds, job1202426.
- SLSQP converged in12–19 iterations. Every proposal meets covalent/source/water
  checks but touches the declared0.20Å boundary; these are constrained proposals.
- Ten original-core native CPCM-DFT checks are running as job1202428,
  four concurrent16-rank tasks on64 allocated CPUs.
- First job1202425 failed at import before model load/energy evaluation. The
  implementation-only fix defers preparation dependencies outside the GPU path.
- The first four native DFT results retain PQQ ordering but narrow the
  XoxF-minus-MxaF gap30.0719→26.2757kcal/mol from the matched starting geometries.
  All endpoint energies decrease; this does not demonstrate accuracy improvement.
  Six alpha/GGR endpoints remain pending at this checkpoint.

## Chemical states

The reusable contextual-water preparation replay has reproduced four archived
proposals and30zero-water identity operations without new molecular calculations.
Scoring replay/negative tests completed: nine tests pass; no new molecular calls.
The reusable component is documented in diagnostics/contextual_water_20260919.

Next approved finite test: all existing site-water/carboxylate-O contacts≤2.60Å
in the full-water contextual alpha preparations (one1F6S and two6IP9contacts),
bothCa/La, four new internal-proton-transfer path points percontact:24native
EnGrad calculations with archived neutral centers reused. Same atoms/charge and
all other coordinates fixed, so the proton reservoir cancels. Exact fractions,
H identity and acid geometry checks are frozen in the owning agent's plan before
energies. An uphill fixed-heavy path does not rule out a separately relaxed basin.
No pH population, entropy or full occupancy claim follows from this electronic
test. The failed wider-basin approximation remains unavailable.

Identity correction:1F6Sand6IP9arebothalpha-lactalbumin structures. Prior wider
water-basin prose mistakenly called6IP9GGR; numerical PDB IDs and energies are
unchanged. The owning agent is correcting interpretation and the vault note.

## Independent electronic diagnostic

Frozen first scope: six DLPNO-CCSD(T1)/TightPNO/CPCM-PTES single points on exact
real GGR1GLG and alpha1F6S original/context-prepared water-H Ca/La cores.
Submitted as job1202429 on64CPUs/256GiB, four16-rank tasks concurrently. The
agent is verifying the actual parsed Hamiltonian before interpreting energies.

Selected reference: diffuse def2-TZVPPD on ligands/La, cc-pwCVTZ onCa,
La46-electron def2 ECP, explicit Ca10/La46 total frozen-core policy,
AutoAux correlation fitting, TightSCF/DefGrid3. No composite-method D3/gCP
addition. This is a finite-basis diagnostic, not a claimed complete-basis truth.
Preserve SCF/correlation/triples terms and diagnostics. A changed prediction
cannot be attributed uniquely to correlation when basis/solvent treatment also
differs. Measure small-core throughput before executing larger PQQ controls.

## Existing noncatalytic guardrail

`diagnostics/benchmark_set_20260915/SCORING_RESULTS_1199508.md` and
`workspaces/benchmark_set_20260915/ready_tasks_v4/manifest.json` retain the
completed50-atom Hans-LanM EF1/EF2/EF3 Ca/La pairs. All rank aboveGGR under the
baseline. These are three ordered sites within one protein-level coupled
folding/dimerization measurement, not three independent site-affinity labels.
They are available for relevant follow-ons without new folds or invented labels.

## Completed scoped review

Commit `e085d55` contains completed response/second-shell code and reports plus
the new parallel-scope record. Parent reran26real-fixture tests successfully.
The detailed~1MBsecond-shell atom inventory remains preserved locally; it was
not added to the compact diagnostic commit. Other agents' dirty changes remain
untouched. No default, published reference or production score was changed.
