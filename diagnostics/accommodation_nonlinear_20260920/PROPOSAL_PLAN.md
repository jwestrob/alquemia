# MACE proposes geometry; the composite selects and scores it

2026-09-20. Separate contained experiment under Jacob's overnight goal. Declared
before these optimizations/results. The running fully composite pilot is preserved.

## Why this comparison

Native DFT now confirms both compressed PLM displacement responses. Native MACE
alone describes their differential work more closely than the composite on these
particular paths; this does not justify deleting solvent from the classifier.
Solvent is useful in the reference classification, but evaluating two quantum
solvent endpoints at every optimizer step may not be necessary. Test whether
MACE can cheaply propose a physical donor arrangement, then let actual composite
energies choose between that proposal and the original. This is a restricted
two-geometry electronic descriptor, not a continuous composite minimum.

## Fixed sources and proposal rule

Use all28 exactly prepared canonical/crystal PQQ contexts from
`reference_preparation_v1/READY.json`, plus the two original sample0 PLM contexts
from the torsion design. Thirty unique contexts,60 metal-specific starts. Keep
the known labels only for reporting; the two PLM labels remain unknown.

Use precisely the same terminal Glu chi3 and actual extra-Asp chi2 mappings as
the nonlinear pilot/reference preparation. All other coordinates, chemistry,
waters, PQQ, charges and spin remain fixed. Start at q=0 and use native float64
OMOL100M vacuum energy/analytic forces with L-BFGS-B and the same angular bounds,
maxiter,maxfun,maxls,gtol and relative-energy ftol as PLAN.md. No new starting
points, tuned springs, geometry targets or score-directed moves.

Keep actual optimizer success, raw active derivatives, boundary distance, geometry
checks, energies and forces. A successful bounded optimizer and supported geometry
can provide a proposal even at the boundary, explicitly flagged; it is never
called an unconstrained or full-protein minimum. An unsuccessful/missing proposal
remains unavailable. No entropy, curvature-based correction or native-DFT call.

## Actual composite selection and score

For each metal evaluate exactly q0 and its one MACE proposal under the existing
primary native-GFN2 ALPBwater-minus-vacuum composite recipe. Reuse matching old
q0 outputs, with actual method/state/coordinates checked. The primary recipe is
the released scorer's original numerical policy; the running full-gradient pilot
uses Tight and must report that difference explicitly in any comparison.

At most120 new GFN2 singlepoints cover the sixty proposals, two media each. No
new q0 calls are needed: all thirty contexts have archived primary q0 energies.
Fresh MACE origin evaluations check replay and provide forces for the optimizer;
all actual model evaluations/failures are counted. If exact q0 compatibility is
missing, report that endpoint unavailable, without silently recomputing it.

Choose the proposal only if its full composite energy is lower than q0 by more
than0.10kcal/mol, the pre-existing solvent-transfer numerical scale. Otherwise
retain q0 as the explicitly evaluated member of this fixed candidate set. This
is energy selection, never label selection. Require both complete members;
solver failure does not silently fall back to a successful baseline score.

Report all original/proposed/selected endpoint energies and
R_selected=E_Ca,selected-E_La,selected, component work, coordinates, boundary flags
and old-band transfer decisions. Frozen bands are a developmental compatibility
check, not a new calibration; no fit on these outputs. Report canonical25,
consumed crystals3 and unknown PLM2 separately, with all denominators.

## Execution and interpretation

Reuse existing native GPU model and ORCA task runners. One warm1GPU/32CPU/
200000MiB allocation for MACE proposals; a finite nativeGFN2 manifest runs on
64CPUs with8concurrent8-rank tasks, using the same primary input recipe. No
project compute/time cap or uncontrolled optimization. Save failed attempts;
no outcome-selected extra poses or retries. No production/default change.

This experiment tests whether a cheap proposal earns physical/classification
utility while avoiding solvent calculations throughout the search. It is not
proof of a DFT minimum or of accuracy on unknown PLM proteins. Compare with the
fully composite pilot on the same four sources once available; differences in
search objective and primary/Tight numerical policy stay explicit. Additional
fold transfer or native candidate validation requires its own recorded manifest.
