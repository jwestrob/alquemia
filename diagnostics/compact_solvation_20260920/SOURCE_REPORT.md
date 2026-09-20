# Exact compact MACE sources and solvent comparison

The inventory contains **all 33 cases, 66 case/representation pairs and 132 real
native vacuum MACE endpoint receipts**. No geometry or scientific energy was
created during source assembly. Source coordinates, charge, spin, checkpoint,
prepared task, execution result and source protocol are pinned. Every MACE
result is native unmasked OMOL total vacuum energy; no whole-protein masked or
subtractive score is substituted.

## Cases, roles and preserved chemistry

- 25 canonical PQQ calibration cases, with their original accession groups and
  functional-class evidence. The source has 14 Ca and 11 La labels.
- 1H4I and 4MAE retain `retrospective_structural_transfer`; 1KB0 retains its
  distinct `retrospective_external_class_transfer` role. All three are consumed.
  The first two overlap calibration accessions; they are not new biological
  observations.
- Two bovine alpha-lactalbumin structures (1F6S, 6IP9) and three GGR structures
  (1GLG, 2FW0, 2FVY) retain the relative alpha-minus-GGR evidence direction.
  Their six comparisons represent one biological direction, not six independent
  tests. No canonical PQQ threshold is used to label those proteins.

Both original core and complete context are available for every case. The alpha
inputs are the actual successful contextual-water preparations: four water-H
coordinates differ between their 1F6S Ca/La endpoints, and six between 6IP9's.
Those metal-specific orientations are preserved exactly in each representation;
no attempt is made to force paired coordinates to match. Their heavy atoms,
atom inventory and protonation remain as archived. The other cases have matching
nonmetal coordinates. No coordination-proposal optimization was mixed into this
static-context inventory.

Original native PQQ bands remain saved under `native_PQQ_bands`; they are used
only on the native score. Their gaps are79.03060796106001 model-kcal (core) and
2.3154591089696623 (context). The native alpha/GGR replay is2/6 for core and6/6
for context. All source records, original preparations and software versions
are retained in INVENTORY.json, which the root executor pins unchanged.

## Comparison and missing-data behavior

`compact_solvation_compare.py compare` verifies low-level output/receipt hashes,
actual endpoint energies, the exact source MACE endpoint/XYZ and declared medium,
solver and electronic state. It reports native MACE, GFN2 vacuum, GFN2 ALPB,
endpoint transfer, Ca-minus-La transfer and composite contrasts separately:

`R_composite = (E_MACE,Ca − E_MACE,La) × eV_to_kcal
             + [(E_ALPB,Ca − E_vac,Ca) − (E_ALPB,La − E_vac,La)] × hartree_to_kcal`.

Each unit conversion occurs once. Larger R remains more La-like. Partial
collections leave the native score available and composite/correction fields
null. Missing, failed and unsupported rows remain in the denominator and failure
history. Duplicate conflicting successful calculations require explicit source
selection rather than silent replacement. Absolute aquo reference and S remain
null. This composite is not a binding free energy or a CPCM correction.

The predeclared class-extrema procedure can generate distinct composite bands
only after all25 calibration rows are available and the gap exceeds the inherited
procedure's0.02 model-kcal numerical criterion. Transfer labels do not set bands.
Core/context and native/ordinary-SCF assessments remain distinct. Raw gap, every
cross-class ordering, within-class spread, three transfers and the complete
alpha/GGR structural matrix remain reportable. Numerical qualification is a
separate root-owned audit; computed preliminary contrasts do not establish it.

## Actual initial pilot replay

The root's real job1203145 supplied eight GFN2 endpoints (two contexts, Ca/La,
ALPB/vacuum). This comparison makes no new scientific calls. All66 output rows
are retained: two composite values are available and64 are null. No calibration
or transfer classification is available from this two-Ca-case pilot.

| Context | Native R | Solvation contribution to R | Composite R |
|---|---:|---:|---:|
| q9z4j7-pqq-la_model | -405347.73205159 | -116.45569669 | -405464.18774828 |
| 1H4I | -405428.47335146 | -78.17463049 | -405506.64798195 |

The same-Ca-class difference contracts80.74129987→42.46023367 model-kcal. This is
preliminary observed dispersion, not evidence of improved La/Ca discrimination;
there is no La example in this pilot. Native xTB does not print a separately
identified ALPB component/settings block here. Inputs, converged totals and
actual exported properties support the matched calculation; component and
numerical qualification are updated below from the actual root audit.

## Tests and use

SOURCE_TESTS_FINAL.txt records eight real-source/parser/algebra tests. They replay
all132 receipts, preserve the exact native2/6→6/6 result, retain null corrections
without low-level outputs, check the actual eight-call GFN2 pilot's units/sign and
identity, reject explicitly corrupted copies of real receipts/collections, and retain the
actual failed numerical cross-check alongside available native results.
No invented successful scientific output is used. The first preflight exposed a
test assumption that all three transfers shared one role string; the assertion
was corrected to preserve the distinct1KB0 role. INVENTORY.json did not change.

Run from the repository root (all calculations remain the root executor's task):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/compact_solvation_compare.py compare --inventory diagnostics/compact_solvation_20260920/INVENTORY.json --collection workspaces/compact_solvation_20260920/pilot_v1/audited_collection_v2.json --numerical-collection workspaces/compact_solvation_20260920/numerical_v1/failed_collection_v1.json --solver native --output workspaces/compact_solvation_20260920/pilot_comparison_replay_v1.json
```

The collector accepts repeated `--collection` arguments, reuses identical pinned
receipts and refuses conflicting duplicates. A full-panel collection can be
passed through the same interface. Creation refuses existing output paths.
No baseline/default change, commit, push or new scientific submission was made
by the source/comparison subagent. Root owns orchestration and final promotion
judgments; see PLAN.md for authorization, numerical tolerances and execution.

## Final source/comparison checkpoint

PILOT_COMPARISON_v3.json is the latest source-agent pilot comparison. It pins the
root's audited native collection v2 and actual failed ordinary-SCF collection.
Native charge sums, valence states and parameter exports passed the root's
source audit; the largest absolute atomic charge was about0.85e. Separately
printed ALPB terms/radii are unavailable, not fabricated or filled with zero.
All eight ordinary-TightSCF endpoints failed convergence after124iterations;
no numerical-comparison energies are accepted. The report therefore records
`numerical_qualification = numerical_crosscheck_unavailable`, retains all eight
failure artifacts and leaves numerical energy differences null. Native pilot
scores remain available as preliminary descriptors; no solver was substituted.

Root alone owns the continuing256primary-call expansion. The source inventory
and its132archived native endpoints are unchanged. The same comparison interface
will collect full-panel raw contrasts, distinct protocol bands, transfers and
all paired core/context effects while carrying the failed numerical check.
The source/comparison subagent launched zero scientific calculations. Root's
execution receipts, not this inventory assembly, determine aggregate compute cost.
