# Native scout and saved collection handoff

Frozen runner: `REFERENCE.py` ea95da20; panel
`workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1/prepared_v2/PANEL.json`
726c877f; independent review98542e00. READYdcb74a65 contains executable argv.
Root owns `run_v1`, releases, submission and completion watchers. This note makes
no native call and does not change a frozen artifact.

## Four actual scouts

Each result belongs at
`run_v1/tasks/<cell_id>/execution/{START.json,endpoint.inp,endpoint.out,endpoint.err,RESULT.json}`.

- `control__water__Dy__complex`
- `control__water__Dy__ion_ghost_water`
- `control__water__Dy__water_ghost_ion`
- `control__cc__acetamide__La__complex`

Inspect these actual outputs once root reports completion. Verify four-rank
allocation and grant-derived MaxCore; exact executed input; RHF/NoFrozenCore;
90/201 orbital functions and18/8/10/40 electrons; Dy369/La351 metal AuxC functions;
real/ghost ECP and nuclear inventories; SCF, localization, local-MP2, CC and
iterative-triples traces; energies/closure, warnings and actual costs.
Native admission remains distinct from the later18-cell approximation test.

## Remaining14 controls

The frozen `release_contract` requires stage `controls_remaining`, the exact
other14 ordered IDs from PANEL, and a pinned `control_scout_review` with:

```
status: passed_control_scout_review
panel: exact prepared_v2 PANEL pin
results: {each of the four scout cell IDs: exact native RESULT pin}
```

The gate re-admits all four scout raw outputs through the frozen parser. It
cannot accept pending or failed original RESULT files. If an actual native job
completed but the unseen output grammar fails parsing, preserve the original
source, RESULT and all raw files. A versioned saved-output recollection must
bind those authorities explicitly. Any subsequent execution/collection handoff
must then bind that reviewed parser recovery; do not edit the old RESULT, weaken
the numerical gate, or repeat the molecule. No speculative recovery is prepared.

## Full18 approximation qualification

Copy/fill `CONTROL_RESULT_INVENTORY_TEMPLATE.json` as root's
`run_v1/CONTROL_RESULT_INVENTORY.json`. It already lists the exact18 native result
paths. Missing files remain unavailable; the collector retains all18 slots.
READY's `control_collection_argv` writes a fresh
`workspaces/.../coordination_reference_uncertainty_v1/control_qualification_v1`.
Only `status=local_controls_passed`, `complete=18`, `all_controls_passed=true`
can support the next stage. There are35 individual HF/CP/AuxC/exchange checks.

## Cluster stages and full26 collection

`cluster_scout` releases only `cluster__tight__aquo9_central_Dy` and binds the
actual18-control qualification. The remaining7 require that qualification plus
`cluster_scout_review` containing status `passed_cluster_scout_review`, the exact
panel and native scout RESULT pins, and `cost_and_memory_accepted=true` from root.

`FULL_RESULT_INVENTORY_TEMPLATE.json` supplies all26 exact paths for the final
saved collector. READY's `full_collection_argv` writes a fresh
`workspaces/.../coordination_reference_uncertainty_v1/full_collection_v1`.
Both four-cell tier contrasts need all26 admitted records; no partial quartet,
mixed-method CP/B, fitted correction or target replacement follows this panel.

The saved tests include a positive18-row collector regression and a fully
missing26-row case. Those synthetic receipts prove reporting mechanics only;
they are not native/scientific qualification authorities.
