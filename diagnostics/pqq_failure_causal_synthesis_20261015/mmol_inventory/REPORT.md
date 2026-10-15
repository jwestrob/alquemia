# MMOL1770 Ca sample0: a specific donor-response gap remains

**Six fixed-geometry vacuum DFT cells can directly test the response behind this
inconclusive case.** The exact current pool and all cheap-model receipts exist;
no compatible DFT energy or atomic gradient was located in the declared archive
search. Run the target's Ca/La origins first, then its four saved proposal cells
only after the actual reference execution/admission gate. The first available
same-conditioning control is sample1; it remains inventoried, not scheduled.

## What the existing calculation actually selected

The target is `mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-0`, in the
fixed `mmol_1770-pqq-la_model__envelope_5063ec4e8b07a456` context. Both sources
have193 atoms/13 caps, dry oxidized-PQQ3−, Ca charge−3/La charge−2, and singlet
multiplicity. Exact source, state, atom/cap maps, origins, both proposal receipts,
MACE energies/forces, xTB vacuum/water charge files and qualified frozen-charge
native ALPB results are pinned in [INVENTORY.json](INVENTORY.json). All12 cells
exist; metadata inspection verified208 referenced file hashes. No numerical
geometry analysis, scoring, scientific imports or molecular calls ran here.

The current operational chooser selects `adaptive_Ca` for Ca and `adaptive_La`
for La. Saved results, in kcal/mol:

| Target donor arrangement | Ca work | La work |
|---|---:|---:|
| origin | 0 | 0 |
| adaptive_Ca | −16.944880 | −8.529883 |
| adaptive_La | −10.260699 | −12.758787 |

The mixed-selection contrast changes **−4.186093**, comprising MACE vacuum
**−4.004941** and native ALPB at frozen vacuum charges **−0.181152**. This is not
a single common-path differential. The six-cell reference retains both metals
at all three geometries, so it can distinguish each common path from the actual
mixed selection. Values are copied from completed score job4037, not recomputed
on login: `workspaces/pqq_inconclusive_review_20261014/scores_v1/execution/RESULT.json`.

The target's four selected coordinates are Glu200 chi2/chi3, Asn318 chi2 and
Asp362 chi1. Metal/PQQ/scaffold coordinates remain fixed. The Ca-generated
proposal reaches the existing0.8Å heavy-atom bound and has a stored selected-mode
MACE residual3.400850kcal/mol/Å; the La proposal reaches0.501266Å and its selected
residual is0.005612. Existing geometry job4038 found Glu OE1 shortening
2.322→2.151Å and Asn OD1 2.840→2.404Å along the Ca proposal. This identifies the
motion to test; a constrained nonzero gradient does not establish optimizer
failure or an excessive physical response. See the exact saved
[geometry report](../../pqq_inconclusive_review_20261014/geometry/REPORT.md).

## Existing DFT is not a substitute for this test

The complete September20 fold DFT archive contains both metals for target0 and
control1, but they are **80-atom origin cores with CPCM(Water)**, with different
XYZ hashes. Their inputs are singlepoints and their task directories contain no
`.engrad` files. The inventory retains their actual energy/output/execution pins
as explicitly incompatible evidence, never as reuse. Current193-atom exact XYZ
hashes are absent from the four inspected completed full-fold, torsion, PedE and
ExaA native-reference manifests. This is a bounded archive search, not a claim
that every file on the filesystem was searched.

The earlier MMOL sample1 Asp362±0.2rad profile used MACE and GFN2 only:
[report](../../accommodation_reference_profile_20260920/REPORT.md).
The completed native donor validation covered4MAE and two PLM extra-Asp paths,
explicitly excluding Glu/coupled-mode qualification:
[DFT report](../../accommodation_torsion_20260920/FINAL_DFT_REPORT.md).
Recent ExaA vacuum DFT covers different sources/coordinates and a La-generated
path: [final report](../../exaa_donor_response_20261014/FINAL_REPORT.md).
Neither validates MMOL sample0's coupled four-mode response. Saved SCF orbital
gradients are not nuclear-coordinate gradients; all12 current MACE force arrays
are available, but an exact native DFT atomic-gradient comparison is unavailable.

## Minimum reference and reuse directions

Use the existing native ORCA6.1.1 r2SCAN-3c/DefGrid3/TightSCF vacuum execution and
corrected TRAH collection machinery from
`diagnostics/exaa_donor_response_20261014/reference/`. A new finite adapter must
change its ExaA-specific identity/count assertions to the pinned target,
193 atoms/13 caps and six tasks; preserve source/state/paired-coordinate gates,
raw native warnings, actual convergence criteria, immutable attempts and
allocation-aware12MPI execution. The existing real192-atom parser fixture may
remain a parser fixture; it must not become the target atom-count assertion.
No CPCM or frozen-charge solvent term belongs in the initial vacuum-versus-vacuum
comparison. No old input, prepared geometry or production scorer is modified.

Inventory rows expose `sample`, `role`, `source_case_id`, `case_id`, `selection_id`,
`prepared_admission`, `state_signature`, `maps`, `source_preparation`,
`native_proposals`, `proposal_manifest`, origin `endpoints`, and
`matrix[Ca|La][origin|adaptive_Ca|adaptive_La]`. Each cell includes exact XYZ/MACE/
forces/native-result pins, charge, multiplicity, atom count and archived GFN
valence-electron count596. Native DFT all-electron/ECP counts must be checked by
the actual Slurm preflight;596 is not an ORCA basis-electron count.

The nearby192-atom ExaA eight-cell reference cost522,696 allocated CPU-seconds,
each endpoint using12MPI, with the last Ca proposal taking4,973 seconds. That is
an actual prior cost anchor, not an MMOL timing measurement. A six-cell round is
roughly0.4million allocated CPU-seconds if comparable; measure the two actual
origins before forecasting the four proposal cells. Shared resource policies
and live job ownership remain root's responsibility.

Control sample1 is the first consumed available Ca-conditioned peer by index.
It has the historical2.1016Å extra-Asp contact and both current proposals hit the
0.8Å bound. Its current accommodation shift is+17.014809kcal/mol and it is
La-called. Retain that source-specific ambiguity; do not replace it with a
more favorable peer. A control round is conditional on the target result and
is not part of the proposed initial six-cell execution. No new labels, bands,
optimized coordinates, biological preference claims or production promotion
follow from this inventory.
