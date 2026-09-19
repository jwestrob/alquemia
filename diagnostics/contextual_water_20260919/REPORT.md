# Reusable contextual water preparation — complete replay

**The demonstrated water-H preparation is now usable with explicit prepared-case
inputs, without pilot case-name logic or a prior DFT calculation.** Production
defaults remain unchanged. This component preserves retained water identities,
oxygen coordinates, all nonwater coordinates, protonation and charge.

The implementation retains the established 3.5 Å polar context, source-backed
complete fragments, two native OMOL orientation starts and rigid-water optimizer.
Optional context groups explicitly share the fragment inventory across compatible
indexed protein structures. Source atoms outside the prepared protein assembly
are rejected. Only the resulting site-water H coordinates enter the original
scanner structure. Context outer waters are not silently added to it.

## Actual replay

Four completed native MACE proposals and eight completed masked MACE endpoint
energies were reused exactly. Two alpha-lactalbumin preparations reproduce the
previous scanner coordinates and four-state descriptors:

| Case | Ca−La interaction, model-kcal | Detached environment difference |
|---|---:|---:|
| ALPHA_1F6S | 64.786514287 | −2.113135965 |
| ALPHA_6IP9 | 64.945305994 | −2.663771489 |

The required expression is `(Ca_bound−Ca_detached)−(La_bound−La_detached)`.
Endpoint-specific H coordinates prevent canceling the detached terms. Components
close within 8e−15 model-kcal. These are existing uncalibrated descriptors, not
new accuracy evidence. The demonstrated alpha-versus-GGR comparison remains one
biological comparison; 6IP9 itself is alpha-lactalbumin.

All 30 zero-site-water inputs (25 canonical PQQ, two PQQ crystal controls and three
GGR structures) return their original endpoint records unchanged. They use the
original scorer/cache; this component does not invent a zero score.

Nine real-fixture tests pass, zero skips. They check exact replay, arbitrary case
identifiers, oxygen/nonwater/charge invariants, PQQ identity, rejected unconverged
proposals, an explicitly corrupted real oxygen coordinate, missing detached
receipts and optimizer cache invalidation. Frozen runner dry-runs pass in the CPU
driver and MACE environment. Both execute operations correctly return zero work
for the completed manifests. Prepare/collect/report were actually exercised.

**New scientific evaluations: zero. New allocations: zero.** Measured preparation
of 32 inputs took 35.67 wall seconds / 32.95 user CPU seconds / 1.65 system CPU
seconds on the login host; this is an implementation check, not a production
hardware benchmark. Earlier incomplete preparation directories are preserved.

## Use and limits

The component accepts existing whole-chain `preparation.json` records. Wet cases
require the source-backed peptide-amide v3 core, atom graph, protonated source,
topology and complete existing site waters. It is not a raw-CIF water-placement
tool. Missing waters, oxygen locations, inventory and occupancy remain upstream
scientific decisions. Water-bearing PQQ/cofactor contexts currently fail explicitly
until a compatible adapter exists. Existing dry PQQ inputs need no adapter.

New-case execution uses the existing pinned model/environment and runner through
an isolated dispatch snapshot; it has not been scientifically rerun on a new
protein in this task. The frozen proven optimizer is copied from the actual
qualified source manifest, avoiding replacement by later unrelated live edits.
No daemon or second scheduler workflow was installed.

New protocol IDs:

- `source_context_rigid_water_H_native_omol_preparation_v1`
- `masked_omol_context_prepared_water_four_state_v1`

Prepared/scored exports retain explicit unavailable states, all four endpoint
energies, source mappings and receipts. Occupancy probabilities and calibrated
decisions remain null. See [commands](COMMANDS.md), [preparation tables](preparation_export_v1/TABLES.md)
and [score tables](score_export_v1/TABLES.md).
