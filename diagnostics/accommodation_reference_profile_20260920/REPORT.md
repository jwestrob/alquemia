# Geometry-selected reference accommodation: completed

**The frozen three-point rule weakens a native MACE misclassification but does
not add a new correct decisive call.** The solvent-composite model was already
correct and gains separation from its old La boundary. Four new OMOL energies
and all eight new GFN2 singlepoints completed; no DFT was added. Three actual
fixture/replay tests pass, zero skips.

## Fixed result

Source: `mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-1`, a consumed
La-associated PQQ functional reference folded with Ca. Its extra Asp362 was the
sole <2.2 Å extra-Asp case among the 233 prepared reference folds. Selection came
from geometry; the native misclassification had already been inspected, so this
is development evidence. It is not a new direct-affinity label or blind test.

Both models and both metal endpoints select **−0.2 rad** from the predeclared
Asp362 chi2 grid {−0.2, 0, +0.2}. The positive displacement is retained and worsens
the La-like contrast; no favorable subset or band adjustment was made.

| Model | Ca lowering | La lowering | Change in Ca−La | Old-band transfer: q0 → grid |
|---|---:|---:|---:|---|
| Native OMOL | −1.034 | −7.319 | +6.285 | Ca → indeterminate |
| OMOL + GFN2 ALPB−vacuum | −1.715 | −7.518 | +5.803 | La → La |

Energies are kcal/mol. Native's grid value remains **0.4521** below the frozen
La boundary; no angle extension or threshold adjustment was performed to cross
it. Composite's La-side margin increases **3.4982→9.3012 kcal/mol**. This is one
consumed structural case and an explicit transfer test of old single-geometry
bands, not calibration or broad validation of the new descriptor.

The selected motion increases nearest Asp O–metal distance from 2.1016 to2.1968 Å;
the other O changes 3.3558→3.2874 Å. It remains below the existing 2.2 Å warning.
Both selected angles sit on the finite grid boundary. Neither a physical minimum,
a warning-free structure nor a thermodynamic population has been established.

## Meaning and limits

`R_grid = min_q E_Ca(q) − min_q E_La(q)` is a bounded local-accommodation descriptor.
All six endpoint/grid combinations are available; each model requires its full
grid before evaluating the descriptor. Metal-specific energy lowering is the
selection rule, not agreement with a biological label. Source Ca conditioning,
scaffold, charge, PQQ, all other coordinates and water inventory are retained.
The same whole-carboxylate coordinate measure is used for both endpoints.

The response moves in the favorable direction for this reference, but the
composite already succeeded before it. The earlier 16-task native DFT validation
on the independent frozen 4MAE/PLM contexts is still separate and pending; these
new cheap values do not qualify the energy surface. No free-energy or entropy
claim, default change or new high-level calculation occurred.

## Actual execution and cost

- Exact q0 reuse: two OMOL receipts from `folds_v1/mace_v3` and four successful
  GFN2 receipts from executed `solvent_shards_v2/shard_2`, with source manifest,
  task identity, actual coordinates, charge, recipe and native-state checks.
- New OMOL job 1203879: 4/4 values;10 allocated seconds ×32 CPUs, 10 GPU-seconds.
  Warm process 2.362 s including 1.171 s model load; peak allocated VRAM
  2,944,835,584 bytes.
- New GFN2 job 1203880: 8/8 values;28 allocated seconds ×64 CPUs, zero GPUs.
- Total incremental allocation: **2,112 core-seconds and 10 GPU-seconds**.
  q0 campaign cost belongs to the root campaign; local preparation/analysis/test
  cost is not included. Zero Slurm batch MaxRSS records are unavailable.

[PLAN.md](PLAN.md) records the frozen rule. [TESTS.txt](TESTS.txt) separates actual
artifact/algebra checks from the 12 executed model calls. Complete unrounded
results are `workspaces/accommodation_reference_profile_20260920/result_v1.json`;
all inputs/receipts are under `prepared_v1/`. Production remains unchanged.
