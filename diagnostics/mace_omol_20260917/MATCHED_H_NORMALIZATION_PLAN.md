# Match the already established H preparation across the vacuum hybrid

Declared after the original-H vacuum hybrid failed1/4 ordering, before new
normalized-core coordinates or energies. Active-goal authorization applies.
The failed original-H candidate remains immutable. This experiment addresses
an independently documented geometry defect; no coordinates are chosen using
scores, and a favorable result is not a validity test for the H preparation.

## Fixed physical preparation

Use exactly the pre-existing normalized `physical_atoms` in the three original
global preparations referenced by mechanics preparedV2: GGR_1GLG, ALPHA_1F6S,
ALPHA_6IP9. Their pinned preparations predate both recent vacuum experiments.
Keep their existing ff19SB radial protein-H correction and recorded water-H
preparation; do not introduce a new H rule, optimization or rotation.

Transfer those existing source H coordinates into all four original core
representations(GGR_extended/GGR_connected/ALPHA_1F6S/ALPHA_6IP9) through the
recorded source-atom graph. Preserve every heavy coordinate, element, donor,
formal charge, multiplicity, explicit water, proton identity and atom order.
Synthetic sigma-link H coordinates remain fixed because their heavy anchors
remain fixed. Reject an unmapped/ambiguous H or any changed heavy atom.

The read-only compatibility audit is
`workspaces/mace_omol_20260917/normalized_H_reuse_audit_v1.json`.
It finds21/47/16/18 changed source H atoms in extended GGR/connected GGR/
alpha1F6S/alpha6IP9, respectively. No new coordinates or scores were produced
by that audit. All six exact normalized whole endpoints already have accepted
masked-MACE receipts. Alpha1F6S has duplicate numerical-control matches:
require their energies to agree within the existing0.01model-kcal tolerance,
then use the lexicographically first task ID, without looking at a class label.
Retain all duplicate/reuse provenance. A mismatch blocks reuse.

Protocol: `masked_omol_matched_normalized_H_vacuum_hybrid_v1`.

## Finite new calculations

Eight native ORCA6.1.1 vacuum r2SCAN-3c/NoAutostart/DefGrid3/TightSCF/EnGrad
endpoints, one Ca/La pair for each of the four normalized cores. Eight matching
raw-zero masked OMOL energy-only core calls with the unchanged checkpoint,
float64 and existing qualified adapter. Reuse six exact normalized whole
energies; no new whole-protein call when reuse passes. No new cofactor chemistry,
DFT functional, basis/ECP, solvent solver, training, numerical DFT gradient,
Hessian, relaxation or entropy. Gradients are retained but not converted to a
mechanical correction or projected through an incompatible old-H Jacobian.

Use the existing runners and policies:64CPU/four16-rank ORCA tasks; oneA5000/
16CPU/64474MiB for eight small MACE cores. Comparable recent work costs about
ten minutes of64CPU wall time and a short GPU allocation. These are expectations,
not limits. There is no application CPU/time budget. Record every actual attempt.

## Same expression and fixed gates

H_M = E_DFT,vacuum(normalized_core,M)
      + T_mask(normalized_full,M) - T_mask(normalized_core,M).

R_H=H_Ca-H_La. Hartree/eV converted once. Preserve all components and compare
with original-H values as a named preparation effect, not additional physics.
No aqueous reference, calibrated zero or inherited decision band exists.

Retain native-readout closure<=0.01model kcal, direct/component algebra<=1e-7
kcal-scale, GGR partition diagnostic<=2kcal-scale, and all four alpha-minus-GGR
differences>0.02kcal-scale. Both GGR representations and both alpha structures
remain in the denominator. Report missing/unscorable cases, all failures and
unrounded values. No threshold change, selected core, fitted coefficient or
rescue of the previous candidate. These remain two consumed biological groups,
with qualified cross-study direction comparisons rather than shared-condition
affinity data. A pass would still require broader real-structure tests.

Prepare, dry-run, execute, collect, compare and test using pinned real artifacts.
Scientific integration tests stay explicitly unrun until actual outputs exist.
Keep baseline/default, old references and concurrent work untouched. No push.
