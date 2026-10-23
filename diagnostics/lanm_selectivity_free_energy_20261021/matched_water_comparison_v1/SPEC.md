# Saved comparison of matched water removal

Declared 21 October 2026. One local numerical thread; zero molecular calls or
scheduler submissions. This collector owns no molecular execution or recovery.

## Inputs and admission

Use the immutable six-cell matched-water-removal preparation, its frozen
REFERENCE.py and MACE bindings, the admitted hydrated comparison and original
dry comparison. Root's new MACE result is under
`matched_water_removal_v1/run_v1/mace/execution/RESULT.json`; six native results
are under `run_v1/reference/tasks/CELL/execution/RESULT.json`. Later AMOEBA6
results are admitted only after their configuration is explicitly pinned.

Verify exact input/source/result hashes, XYZ atom identities/order/coordinates,
charges, physical versus effective multiplicities, and native admission. Check
the actual Tinker XYZ and raw energy/gradient outputs as well as source geometry.
Reconstruct ORCA forces from native engrad using minus gradient times
627.509474/0.529177210903; MACE energy/force conversion is 23.060547830619
kcal/mol per eV. Preserve MACE state metadata and ECP physical-state limits.
No model fitting, changed scientific gate or default promotion is included.

Missing files, unavailable rows and unbound optional configurations remain
explicit statuses without energies. A failed provenance or numeric-consistency
check fails collection; it is not converted into a missing row or zero energy.
Partial collections receive a fresh output directory and retain their snapshot
hashes. They may report complete individual contrasts while listing every
missing term. No result is overwritten when later endpoints arrive.

## Separate contrasts

For each donor/model, report hydrated and matched-dry La/Dy exchange using the
**same saved aquo8 endpoints**:

`C = E(Dy,donor) + E(La,aquo8) − E(La,donor) − E(Dy,aquo8)`.

Then report `Cwet−Cmatched_dry`. The aqueous terms cancel algebraically; this
equals the Dy-minus-La difference in adding the fixed six or seven waters to
the identical surviving metal–ligand geometry. Verify both direct expressions
and their element/charge balance. Free-water energies cancel between metals.
This is a total many-body energy response, not an isolated ligand–water term,
an equilibrium hydration free energy or a protein affinity.

In a separate section, report `Cmatched_dry−Coriginal_dry` only when the original
dry endpoints exist. Reassemble Coriginal_dry with the **same aquo8 endpoints**;
never subtract the older water-monomer-referenced summary. This quantifies the
separate dry donor-placement change. Missing historical reference cells do not
block the new hydration contrast and are not authorized as new calculations.

The exact rigid-ligand derivative is `−u·sum(F_i)` over all ligand atoms. Report
wet, matched-dry and their difference for each metal, and the Dy-minus-La change.
The ligand coordinate has units Å. Do not use the stored metal-only projection
or normalize this whole-ligand translation derivative by atom count.

## Interpretation and preparation tests

Report model-minus-frozen-f-reference contrasts only for jointly available
terms, retaining the unqualified reference/state/ensemble caveats. No success
threshold is chosen from these diagnostic results. Validate sign and missing
term propagation with synthetic balanced energies, unit/force parsing against
saved outputs, and rejection of wrong charge/geometry/provenance. These tests
make no native calls and do not establish physical accuracy.
