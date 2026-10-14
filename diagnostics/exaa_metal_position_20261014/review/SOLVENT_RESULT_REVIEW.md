# Placement changes the donor response; full-pool benefit remains unestablished

Independent read-only review of completed solvent job4132. **The reciprocal effect survives the qualified composite energy, but it does not yet justify implementing a joint optimizer.** Finish the saved full-three-to-five-pool comparison and metal-force diagnostics first. A large change in work along a selected path is not the same as lowering the best available energy beyond the existing candidates.

## Actual response

All eight new cells are available; eight original cells were reused. The energy is MACE vacuum plus native ALPB evaluated at each endpoint's GFN2 vacuum charges. For each placement separately, work is saved `adaptive_La` energy minus that placement's origin energy; the differential is Ca work minus La work.

| Host | Original differential | Transferred differential | Change | MACE part of change | Solvent part of change |
|---|---:|---:|---:|---:|---:|
| Sample0 control | +6.539282 | +11.261659 | +4.722377 | +3.373337 | +1.349040 |
| Sample3 target | +18.491708 | +8.597700 | −9.894008 | −9.328603 | −0.565406 |

Values are kcal/mol copied from the completed collection. The target-minus-control differential changes from +11.952426 to −2.663959. The solvent contribution preserves the reciprocal direction and is smaller than the MACE contribution. This supports sensitivity to the fixed source-conditioned metal point within the declared local model. It does not identify a uniquely correct point or show that either transplanted endpoint is an equilibrium geometry.

## Why this does not yet demonstrate a useful search improvement

The secondary four-geometry calculation measures minima relative to the **original** origin, so its reference differs appropriately from the primary placement-specific donor works. In the control, both metals retain the original `adaptive_La` geometry; adding the two transferred geometries changes no selected energy. In the target, only Ca changes to the transferred origin: its selected work is −9.039922 versus −7.898852 in the restricted original two-geometry pool. La retains original `adaptive_La` at −26.390560. The resulting differential improvement is only −1.141070 kcal/mol, despite the much larger primary path change.

That comparison omits the real historical `adaptive_Ca` candidate. The already completed packing comparison reports target original-context Ca work **−9.698819** at `adaptive_Ca`, below the transferred-origin value. Thus even the apparent 1.141070 benefit may disappear against the correct full baseline. This is a reason to complete the prepared saved-data comparison, not a substitute for its exact identity, state and numerical checks. No full-three-to-five result is claimed in this note.

The force diagnostic can establish whether a retained metal force points along a physically admissible displacement at a saved endpoint. Vacuum forces alone cannot certify descent of the solvent-added energy, convergence to a useful joint minimum, or source-independent classification. A local force can motivate a small subsequent physical test if the full-pool result leaves a concrete question; it should not automatically trigger optimizer engineering.

## Qualification, limitations and cost

The executed recipe retains the pinned standalone GFN2 parameters, independent component-v2 checks, native water-diagonal replay and informative vacuum-charge precision checks. Original absolute-point-v1 failures remain visible separately; they were not converted into passes. The corrected full ORCA JSON is provenance only for this route, which makes no ORCA call. Numerical agreement does not validate the dry oxidized-PQQ state or a biological affinity.

No thresholds, bands, class calls, new optimizer, DFT evaluation or source selection were introduced. The existing fixed geometric intervention preserves host nonmetal coordinates, state and caps. New transferred geometries still lack an independent DFT test. Both sources are consumed development examples, not held-out validation.

The stage receipt records **166.6373 s on seven allocated CPUs**, 1166.4611 allocated CPU-s, 32 started/terminal SCC calls and 16 native component evaluations. There are zero new MACE or DFT calls. These are this continuation's costs, not end-to-end new-site latency or measured utilization.

## Authoritative records

- `workspaces/exaa_metal_position_20261014/solvent_v1/execution/FINAL_COLLECTION.json`, SHA256 `81a4167f726e750ee89d22de5443406a1f6cc596b31240a2066f6f25a3a6d2d0`.
- Stage cost: the adjacent `EXECUTION.json`; frozen config SHA256 `0f162cb433686011436a4c9f04a71ebded1e6c3e8de112383032e7ec6933d708`.
- Existing full historical pool: `diagnostics/exaa_source_packing_20261014/packing/FINAL_REPORT.md` and its pinned actual comparison SHA256 `664a35cf2fc117f77acadcf45f5e2a8fa5409a83952f27874e987dafa010ac3c`.

This review reads completed artifacts only. It performs no scientific calculation on the login host, molecular call, refit or submission. Root owns the subsequent saved-data execution and scientific decision.
