# Component field name in the frozen initial collection

The immutable `prepared_v2/collection_1202082.json` field named
`components_delta_R_kcal_mol` contains **expanded Ca-minus-La component contrasts**,
not expanded-minus-core changes. The final paired report reparses actual baseline
and expanded outputs and calculates the changes correctly under
`component_delta_R_kcal_mol` in `result_v1/result.json`.

The live collector now calls its expanded-only field `expanded_R_components_kcal_mol`.
Original scientific outputs, collection and frozen implementation are preserved.
No energy calculation or numeric result was changed to repair this field name.
The final actual-fixture check verifies SCF + dispersion + gCP component changes
close to the total change. CPCM dielectric is already inside SCF and is not added
a second time.
