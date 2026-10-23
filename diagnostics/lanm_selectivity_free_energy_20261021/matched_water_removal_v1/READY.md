# Six matched water-removal controls are prepared

The matching workspace's `prepared_v1/RESULT.json` records successful saved-only
preparation and validation. **No molecular calculation or submission occurred.**

| Donor | Retained atoms | Waters removed | Charge | La/Dy physical multiplicity |
| --- | ---: | ---: | ---: | --- |
| Monodentate acetate | 8 | 7 | +2 | 1 / 6 |
| Bidentate acetate | 8 | 6 | +2 | 1 / 6 |
| Acetamide | 10 | 7 | +3 | 1 / 6 |

Every surviving XYZ atom line is unchanged from its hydrated central parent.
All six ECP cells have 40 explicit electrons and effective multiplicity 1.
La/Dy paired coordinates, state bookkeeping and the six-cell denominator pass.
The frozen reference executor's `validate` function passed six geometries and
rendered six decks. The frozen MACE reader independently recovered identical
symbols and coordinates. Preparation used one thread and 0.339 s internal wall
time; there were zero native, MACE or GPU calls.

## Root-ready artifacts

Under `workspaces/lanm_selectivity_free_energy_20261021/matched_water_removal_v1/prepared_v1/`:

- `PANEL.json`: source records, deletion maps and pinned existing results.
- `NATIVE_PANEL.json`: six-cell binding to unchanged `REFERENCE.py`.
- `MACE_CONFIG.json`: six calls through unchanged MACE runner/helper/config.
- `REFERENCE_VALIDATION.json`: actual saved-only executor validation.
- `xyz/`: six exact coordinate subsets.
- `preview_inputs/`: exact method previews at four ranks and 3000 MiB MaxCore;
  these are previews, while runtime MaxCore comes from the actual memory grant.
- `RESULT.json`: evidence and complete prepared artifact pins.

The reference IDs are `ecp_pbe0__matched_dry_<donor>_central_<metal>`, for the
three donors and La/Dy. Root can use the existing reference `RUN.sbatch` with
the unchanged executor's `run --panel NATIVE_PANEL.json --cell ID --out ...`
arguments. The MACE binding uses the existing MACE launch wrapper and
`RUN.py --config MACE_CONFIG.json --output ...`. Do not use the initial
56-cell-specialized reference collector for this panel.

## Interpretation to preserve

Use the same saved aquo8 terms for hydrated and matched-dry contrasts. Their
difference is the metal-dependent cost of adding the fixed surrounding waters;
water and ligand geometry are otherwise matched. It is a total many-body
energy response, not an isolated ligand–water contribution or an affinity.
Native SCF, frozen-f, model accuracy and absent ensemble/entropy limitations
remain unchanged. Root review and submission are still pending.
