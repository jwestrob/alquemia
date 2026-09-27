# Frozen numerical diagnosis — 2026-09-27

Question: is the small joint rigid-transform failure numerical integration error,
and is the measured environmental response stable under refinement?
Declared before new evaluations. Reuse all completed six-cell scout and twenty-cell
force results. No new molecule, state, label, optimization, or environment selection.

Eight new endpoints: original-grid Ca/La A translated by (0.173,0.117,0.231) A;
refined-grid Ca/La at A, B, and A rotated0.37rad about z then translated as above.
Use identical real 1H4I source, QM membership, links, charges, singlets, ECP,
9087 point charges and Thr159 modeled B perturbation. r2SCAN-3c EnGrad,
NoAutostart DefGrid3 TightSCF and DoEQ false remain. Only refined cells additionally
set AngularGrid7, IntAcc7.0, GridPruning Unpruned, HGridReduced false in %method.
Official ORCA6.1 numericalintegration manual documents these grid controls;
archive: workspaces/metal_environment_response_20260926/reference_docs/.
The actual output headers must establish executed settings; unexpected headers
remain unqualified until inspected, not silently accepted.

Same rigid thresholds: energy1e-5Eh, maximum core/external gradient1e-4Eh/bohr.
Original-to-refined individual A/B response and double-difference changes must
be <=0.05kcal/mol, the original numerical uncertainty target. Report all residuals.
Original failed qualification remains failed; refined directional finite differences
are unrun and not inherited. No automatic further grid rounds.

Run eight workers x43 MPI ranks=344, exclusive full node with --mem=0,
verified full RealMemory and standard25% ORCA operational headroom. Normal priority,
no PQQ dependency/startup guard, no test queue. Existing Monday shutdown cutoff.
Collector and actual Codex-queue wake armed. Record actual costs and failures.
This diagnoses reference numerics, not full hybrid correctness or discrimination.
