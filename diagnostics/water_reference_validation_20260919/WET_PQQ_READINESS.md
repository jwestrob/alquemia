# Existing wet-PQQ input readiness

Read-only follow-up requested by the parent after parvalbumin completed. **No
existing preparation in the inspected manifest inventory qualifies as a labeled,
water-bearing PQQ validation input.** No new waters, microstates, calculations or
biological labels were created.

The bounded inventory examined 561 `carve_manifest.json`, `preparation.json` and
`repair_manifest.json` records under `diagnostics/` and `workspaces/`. Of these,
405 explicitly record PQQ atoms/fragments; none records a nonempty retained-water
inventory in the recognized water-inventory fields. These include many protocol
variants and duplicates, not 405 independent structures or biological controls.
The pinned record list is
`workspaces/water_reference_validation_20260919/wet_pqq_inventory_v1/inventory.json`.
This inventory does not claim an exhaustive search of all legacy files, raw
crystallographic waters or published metal-use experiments.

Specific existing records explain the practical gap:

- The canonical 25 PQQ cases and two crystal controls are dry at the scored core;
  their identity checks preserve current calls but cannot validate wet-PQQ
  preparation.
- Q46444/1KB0 uses the frozen dry fixed-core preparation. Its source has 1,016
  waters, but the existing [selector audit](../pqq_q46444_1kb0_external_validation_20260915/SELECTOR_AUDIT.md)
  reports none within 3.6 Å of its metal. Outer source waters are not retained
  coordinating water and cannot be relabeled as such.
- The newly prepared 8GY2 control explicitly has no crystal/synthetic water.
  It is Ca-associated, not a direct La/Ca-affinity label; its biological assembly
  also contains heme. See the [existing report](../benchmark_augmentation_20260918/REPORT.md).
- Other deposited candidates such as 4CVB, 7WMK, 4MH1 and 3DAS remain chemistry,
  assay-mapping or noncanonical-core gates in that report. Their source coordinates
  do not themselves supply a qualified water-bearing PQQ scoring preparation.

To qualify this branch, an input needs a real retained site water with O/H source
mapping and a compatible PQQ/cofactor state, an explicit full-context assembly and
exclusion policy, and evidence tied to the exact protein/state. The current wet
component needs a coherent PQQ/cofactor adapter. Existing dry-core decision bands
do not automatically calibrate a changed water-bearing protocol. Until those
inputs exist, wet-PQQ status remains unsupported; dry-PQQ regression is preserved.
