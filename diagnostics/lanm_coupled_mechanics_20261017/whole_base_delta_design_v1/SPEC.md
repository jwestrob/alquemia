# Whole electronic base with a regional correction: fixed LBT test

Declared 19 October 2026 before any new evaluation for this candidate. The
active mechanics goal authorizes this bounded test. Root owns submissions.
Independent design review is [REVIEW.md](REVIEW.md); exact source and cache
inventory must be pinned before preflight. Existing DFT job 14795 continues.

## Question and energy

Can a whole-system self-consistent electronic base retain scaffold and remote
electronic response while a matched local correction retains the stronger
metal-site description, without the challenged classical scaffold ledger?

For the unchanged physical coordinates R and existing capped regional map Q(R):

`E_candidate = GFN2_vac(full R) + MACE_vac(Q) - GFN2_vac(Q) + fixed_FF_ALPB(R)`

Also report the indispensable no-extra-call comparator:

`E_base = GFN2_vac(full R) + fixed_FF_ALPB(R)`

Whole-system GFN2 replaces the internal scaffold/cross-region base. Do not add
the old MM ledger, embedded-minus-vacuum GFN term, separate Coulomb or additional
dispersion. The regional terms must have identical atoms, caps, coordinates and
charge. Full forces require the original cap Jacobian for the regional difference.

This is a subtractive mechanically embedded model. The vacuum regional
correction does not see the environment field or feed back into the full GFN
density. Fixed-charge ALPB is additive. GFN-level self-consistency is therefore
not mutual high/low polarization or electronic solvation. GFN2 includes
self-consistent multipole electrostatics and density-dependent dispersion
([primary method](https://doi.org/10.1021/acs.jctc.8b01176)); that description is
not evidence of accuracy on these states. Cap cancellation can remain imperfect.

## Frozen systems and denominator

Use the exact four-cell 7CCO before/after × La/Dy reference inventory from
`../lbt_whole_reference/INVENTORY.json`, and the exact four-cell 7CCN inventory
from `../lbt_whole_7ccn_reference_v1/INVENTORY.json`. They are separate fixed
deformations, not interchangeable basins. No alignment, atom replacement beyond
the already-declared metal pairing, coordinate search, protonation change,
source selection or new optimization is allowed.

Each whole system has 256 physical atoms, target index 255 and charge −2. Reuse
the existing baseline region rule: 172 atoms, residues 3–14 with two caps,
regional charge −2. Do not substitute the old 48-atom/charge −1 reference. Bind the
actual physical identities, cap parents and Jacobians for each coordinate pair.
Use unchanged physical MACE multiplicities La 1 / Dy 6. Both GFN terms use the same
qualified published valence model and effective UHF 0/f-in-core convention;
this is explicitly distinct from MACE spin and the DFT ECP implementation.

The complete denominator is eight candidate cells. Source inspection identifies
seven available regional MACE/GFN pairs and all eight fixed-FF ALPB receipts.
The 7CCN after-La isolated GFN component must reuse its exact certified damping
recovery, preserving original failure and certificate. Never repeat its failed
default SCC solve. The only missing regional pair is 7CCO after-Dy. A preflight
may discover more exact reuse, but may not enlarge the declared new-call budget.

## New-call limits and stages

Maximum ten new molecular evaluations across this candidate:

- Eight whole 256 vacuum GFN2 gradient evaluations, one per declared cell.
- One regional 172 native MACE energy/gradient and one isolated GFN2 gradient,
  both for the missing 7CCO after-Dy pair.
- Zero ALPB, MM, embedded GFN, DFT, geometry search or optimizer calls.

First run a saved-data preflight on Slurm. It must bind all eight inventories,
parameter/runtime/model pins, exact compatible reuse and full-force shapes;
verify regional map/cap consistency and source-state bookkeeping; replay raw
admissions for reused components under their original criteria; and preserve
known unavailable arrays as unavailable. No model is loaded in this stage.

Then one whole 256 7CCO before-La GFN scout, on a CPU allocation. Continue only
after native admission, finite full gradient and actual wall time ≤600 s. This is
a measured continuation bound, not a scheduler time override. Preserve failure
or excess-cost outcome; no fallback, solver sweep or automatic retry.

If the scout passes, complete the other three whole 7CCO GFN cells and the two
missing regional components. Reuse the scout exactly. Root may place the single
MACE evaluation in a separate one-GPU allocation; all work uses at most 12 CPUs
per job and actual scheduler defaults. No waiting collector allocation.

Collect the four consumed-7CCO cells first. The remaining four prebound 7CCN
whole-GFN cells may run only if all components admit and the candidate's
absolute per-metal deformation-work discrepancy from whole-PBE0 is ≤1 kcal/mol
for **both** La and Dy. This is a newly frozen development triage rule using the
existing reporting resolution, not a reference uncertainty, force-accuracy
threshold or proof of physical validity. Do not tune it after seeing results.
Either work failure closes this candidate's challenge stage without changing
the running DFT panel. Report the base-only comparator regardless of outcome.

Stop this candidate's new scheduling on its first native/admission/component
closure/budget failure. Preserve all declared slots, partials and successes;
no survivor work comparison. Technical adapter defects may receive separately
declared source-only recovery, never silent changes to numerical protocols.

## Native lineage and arithmetic

Use the exact existing xTB 6.7.1 executable/library/parameter/runtime and La/Dy
qualification 12488, full corrected profile SHA 635e111b… and standalone
parameter text SHA f6f05c6c…. Stage the complete qualified parameter text for
every new GFN call; no native ORCA/default parameter fallback or global shell
reordering. Preserve accuracy 0.0001, 300 K, 500 iterations, UHF 0, no restart and
native gradient/population outputs. Reused certified solver variants keep
their exact separate provenance. New calls use the qualified default settings.

Use existing immutable MACE checkpoint/native float64 route and force sign.
Retain original energy/coordinate/charge/native admission tolerances. Derive
valence/electron/parity checks from the actual full and regional inventories
on Slurm. Whole atom order and regional target-first order have distinct maps;
neither may be assumed from the task name.

Use raw before/after energy differences before conversion and the established
units. Collect GFN whole, regional MACE, negative regional GFN and ALPB work
separately, plus base-only/candidate totals and Dy-minus-La work. Scatter the
regional gradient difference with the full existing cap chain rule, then add
the whole GFN and unchanged ALPB gradients. No component is counted twice.

## Interpretation and delivery

Report all eight slots with failure/stage status, both comparators, exact
reused/new counts, parameter and state lineage, work discrepancies, full 256
force differences and the original fixed chord projections wherever the
independent reference is available. The pending 7CCN reference outcomes are
not used to select this method, region, geometry pair or triage rule.

No optimizer follows this test. Independent derivative qualification and
physical-state/force-accuracy evidence remain necessary for any later response
search. Whole-PBE0 vacuum dianion, frozen-f, residual and state limitations
remain. Two peptide deformations cannot establish mixed-ion Mex physics or
natural-protein selectivity. Preserve the failed MM hybrid and all source
evidence; record actual allocation and application costs; commit scoped results.
