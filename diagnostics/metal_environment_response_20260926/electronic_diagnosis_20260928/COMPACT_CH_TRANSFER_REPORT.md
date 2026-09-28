# Targeted C–H repair transferred to all six compact endpoints

**Prepared, not electronically evaluated.** All three consumed sources now have matched La/Dy compact inputs that differ from their archived counterparts only at real C-bound hydrogens. No original input/result was overwritten and no electronic job was submitted by this preparation.

| Source | Changed real H atoms | Largest displacement Å | Previously severe cap–CB–H angle, repaired |
|---|---:|---:|---:|
| Hans8DQ2 |14|1.202240|109.011898° (was166.711397°)|
| Hans8FNR |14|1.078545|108.668109° (was164.927825°)|
| Mex8FNS |10|0.413971|not selected for a matching defect|

The large displacements reflect orientation repair, not a small thermal fluctuation. These are deterministic prepared states under the separately reviewed whole-source bonded/LJ surrogate; they are not experimental structures, full quantum minima, or equilibrium samples.

## Exact preparation and admission

`workspaces/metal_environment_response_20260926/compact_CH_repair_v1/manifest.json` records protocol `nikasha_compact_CH_repair_preparation_v1`, three `sources` entries, and per-metal `xyz` path/hash records. Each source retains original XYZ, source-map and reviewed-admission hashes, repaired indices, original/new coordinates, parent source IDs, displacements and explicit fixed/paired checks. Six new electronic endpoints are needed if root executes a comparison; the old Dy origin is not reusable after H movement.

Input admission comes from each source's **REVIEWED_ADMISSION.json** in `scaffold_H_repair_v1`. The original serialization failure, unavailable optimizer-success flag and generic signed-volume failure remain archived. The separate chemistry review excludes CH2 label-volume changes from stereocenter inversions and checks actual final stationarity/geometry. This transfer does not rewrite that history or weaken its other geometry gates.

Every moved H matches an actual amide-v3 source ID and has exactly one actual full-source graph neighbor, carbon. The new coordinate is copied from the admitted repaired source. Retained heavy coordinates match the repair frame; old compact coordinates match the archived map. Both metals receive identical repaired coordinates, with unchanged inventory/charge/state hypotheses. Synthetic caps, metal and every other nonselected XYZ atom line are unchanged exactly.

## What remains unrepaired

- Hans8DQ2/Hans8FNR retain original exchangeable Asn83 HD21/HD22 and backbone H90. Their old bond lengths/orientations remain, as do every synthetic cap and all heavy atoms.
- Mex8FNS retains original backbone H91 plus water320/331 H1/H2; water inventory and geometry are unchanged.
- This is **not** generic normalized-H preparation, solvent reorganization, protonation adjustment or whole-protein relaxation. It does not remove every source-strain mechanism.

Two real-artifact tests pass: complete source coverage with exact fixed/paired coordinates, and actual graph/admission membership for all38transferred H atoms. No fabricated energies or scientific executable mocks were used. Root owns any new electronic comparison and its interpretation; the existing source reversal remains the historical result.

Reproduce preparation into a fresh destination:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/metal_environment_compact_h_transfer.py --request diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/COMPACT_CH_TRANSFER_REQUEST.json --output workspaces/metal_environment_response_20260926/compact_CH_repair_reproduction_v1
```
