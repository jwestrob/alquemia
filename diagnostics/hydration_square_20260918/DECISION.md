# Hydration execution decision — 2026-09-18

## Authorization

Jacob approved the matched observed-water experiment with “Yallah!” and then
explicitly removed the stale approval gate: “remove that AGENTS.md thing bro.
I thought we got rid of that.” His standing full discretionary permissions for
contained discriminator improvements apply. The pending choice in AGREEMENT.md
and SOURCE_CONFIG.json is historical; both remain immutable. This decision and
EXECUTION_CONFIG.json resolve it without rewriting those source records.

## Question and fixed scope

Does each experimentally positioned water change the Ca−La electronic contrast,
and how much did preparation-generated water hydrogen geometry affect it?
Inputs are the two already-consumed bovine alpha-lactalbumin cores 1F6S (waters
A211/A212) and 6IP9 (A310/A322/A326), all five waters, one biological group.

Normalize every water to the pinned bulk_water.xyz internal geometry (about
0.9572 A O−H and 104.52 degrees), retaining its oxygen, plane and bisector.
Every protein/cap coordinate stays byte-identical. The same repaired coordinates
are used for Ca and La. This fixes observed 1.163–1.190 A water O−H distances;
it is deterministic preparation, not an optimized water orientation or basin.

Recompute the two full-water parents for both metals (4 endpoints), and delete
each individual water from its repaired parent for both metals (10 endpoints).
Total: 7 paired states, 14 new single points. Original full-parent energies are
reused only to isolate the normalization effect. No additional water states,
optimization, functional change, classifier fit, or replacement reference.

Native ORCA 6.1.1 r2SCAN-3c, CPCM(Water), DefGrid3, NoAutostart; singlets,
La charge 0 and Ca charge −1. Retain the archived SCF policy. Protocol:
`native_r2scan3c_cpcm_single_water_square_v1`, with the explicit
`reference_internal_geometry` preparation policy in each manifest/cache key.

## Outputs and interpretation

Report all unrounded endpoint energies, failures, paired invariants, components,
receipts and costs. For each water report:

`DeltaS_add = (E_Ca(full) − E_Ca(deleted)) − (E_La(full) − E_La(deleted))`.

Positive means retaining that water makes the contrast more La-like. Bulk-water
chemical potential and common aquo offset cancel here; neither is assigned a
value. This diagnostic does not determine equilibrium occupancy, binding free
energy or an absolute classification. Report normalization separately from
water deletion. Algebra closure tolerance: 1e-7 kcal/mol, set before execution.

## Execution

Existing runner, 64 allocated CPUs, four concurrent endpoints of 16 MPI ranks,
256 GB memory. Recent archived endpoints took 74–99 s each on 16 ranks; this
is a small CPU pilot, with actual receipts required and no GPU request.
No artificial compute budget or requested wall-time cap. Scheduler rules still
apply. Keep baseline inputs/results/defaults and all unrelated work untouched.
