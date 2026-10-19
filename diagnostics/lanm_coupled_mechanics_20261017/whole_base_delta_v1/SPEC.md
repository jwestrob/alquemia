# Fixed whole-GFN base plus regional vacuum correction: first execution version

The immutable [root declaration](../whole_base_delta_design_v1/SPEC.md), its
design review, and `FEASIBILITY_INPUTS.json` govern the question and eight-cell
denominator. This adapter implements **zero-call preflight and one CPU scout
only**. Development/challenge execution and saved comparison use separately
reviewed source. This version cannot launch those later stages.

## Bound inputs and native settings

Freeze all eight 256-atom cells from the two original reference inventories,
seven exact regional MACE/isolated-GFN pairs, the single missing 7CCO after-Dy
pair, and all eight existing fixed-charge ALPB receipts. Both source-local
before/after pairs are immutable. Every regional map has 172 atoms, two caps,
and the original full physical-coordinate derivative map. Full and regional
charges are −2; physical MACE multiplicity is La1/Dy6; GFN effective UHF is zero.
Full750 and regional538 valence counts are derived from the exact qualified
parameter profile at Slurm preflight, alongside physical spin-parity checks.

The exact `635e111b…` corrected La/Dy profile, `f6f05c6c…` complete parameter
text, xTB6.7.1 executable, libraries, environment and native admission are
pinned. New calls use vacuum GFN2, accuracy0.0001, electronic temperature300K,
500iterations, no restart, JSON/population/native analytic gradient outputs.
No damping modification or fallback is permitted. Native threads remain one.

The maximum candidate budget remains eight whole-GFN, one regional-GFN and one
regional-MACE call. This executable permits **one** whole-GFN call, exclusively
7CCO before-La. It permits no ALPB, MM, embedded-GFN, DFT, model loading,
optimization or coordinate changes. Subsequent GPU MACE requires separate
source and allocation; no GPU is requested here.

## Saved-only preflight

The Slurm guard executes before numerical/scientific imports. Replay the original
GFN raw admissions and MACE state/energy/force provenance; match their exact
regional XYZ and physical-source bytes. Verify all ALPB raw receipts using the
unchanged whole-reference `admit_alpb`, with its exact sequential unit conversion.
The original after-Dy inventory's geometry-parent alias and separately completed
ALPB supplement are admitted explicitly. Preserve missing whole-MACE force
arrays as missing; none are needed for this native preflight.

Replay the certified 7CCN after-La solver exception using the frozen
`RECOVER.certify_scc` saved-only function and all four real control outputs.
Retain the original failed default result, certificate and primary .05 receipt;
never repeat that failed solve. No HybridEvaluator or complete-pool constructor
is invoked.

Check archived mapping identities/coordinates/caps by the existing V2 replay
policy, which permits only its previously declared bounded derived-Jacobian
roundoff. Archived maps remain the derivative authority. No coordinate or
physical acceptance threshold changes. The candidate regional derivative is a
new assembly, not a reconstruction of the old hybrid total: convert raw MACE
force to gradient with `-f*EV`; convert raw isolated GFN with `(g*HA)/BOHR`;
scatter each through both cap anchors separately, negate the scattered GFN,
then add in that order. Store both full256 components and their sum. No old MM,
embedded-GFN, or full-hybrid gradient is used.

Preflight writes eight explicit states, prepared whole/regional XYZ pins, map
checks, native/solvent reuse and full component arrays. All call counters remain
zero. Failure preserves the eight-slot denominator and blocks the scout.

## Scout and stopping

`run-scout` requires an actual passed preflight with identical scientific binding
and prepared signatures. Fresh output only; root owns one invocation and its
actual scheduler receipt. Stage the complete parameter file in the native
working directory and use the pinned original native runner/admission.

Preserve every raw output and native lifecycle failure. Require accepted SCC,
finite full256 gradient, and native wall time ≤600seconds to publish
`scout_passed`. A failed admission or excess cost stops continuation; no retry,
solver sweep or automatic follow-up. The cost bound is not a scheduler timeout.
The remaining seven whole-GFN slots remain `not_requested` in this version.

## Later interface and interpretation

The mandatory base-only comparator is full-GFN+sameALPB. The candidate adds
regional MACE−isolatedGFN. Future completed rows will retain raw native energies
and four independently mapped gradient components in fixed order:
`GFN_full`, `MACE_region`, `negative_GFN_region`, `fixed_FF_ALPB`.
Per-metal deformation work uses raw native before/after differences before
conversion. The predeclared 7CCO development gate is ≤1kcal/mol absolute work
discrepancy for **both** metals; failure closes the four-call 7CCN challenge.

This is mechanically embedded subtraction, not mutual high/low polarization.
Fixed-charge solvent does not feed back into GFN. No derivative qualification,
physical spin/state accuracy, ground-state stability, model accuracy or optimizer
permission follows from native or scout admission. Existing whole-PBE0 vacuum
dianion, frozen-f, SCF residual and finite-deformation limitations remain.

## Execution

Root reviews and submits `RUN.sbatch` in the owned preflight workspace. It asks
for one CPU, no GPU and no requeue, retaining scheduler default time/memory.
Use the pinned original MACE Python environment for saved-array replay; no model
is imported or loaded. `CONFIG_SCOUT.template.json` has a null future preflight
pin and cannot execute until root binds the actual passed result.
