# Matched fixed-FF rescore of the closed one-water panel

Declared 16 October 2026 before new solvent calls. Root owns submission.
Question: does the observed source-conditioned water-position response persist
under the now-matched energy model? This is a finite model-dependence check of
the existing consumed panel, independent of the Mex continuation scout.

## Exact existing denominator

Freeze the complete sixteen-slot panel in
`lbt3_water_remaining_v1/execution/RESULT.json`, Slurm12670, SHA256
`aa4dc92c0a62d7579f84c4c2f0f6e559f3ab275f8c13417f1b1f07eaea6fa63a`.
Its two reused scout rows come from12664, SHA256
`26aa78c269cd954b63f48ad034771058d72ccd1aae35fd60376a2715466ee18f`.
Preparation12565 is SHA256
`4379954911dea5cf6a36e72c0c80d744f4a7bc9137895e0a7c61a593e2b04cd7`.

Order: sources7CCO,7CCN; poses inner_plane0,inner_plane90,outer_plane0,
outer_plane90; metals La,Dy. Preserve all eight7CCO cells and four7CCN outer
cells. The four7CCN inner cells remain geometry-rejected and unscored. No
replacement, new placement, new water, occupancy change, changed contact rule,
source selection, optimization or geometry calculation.

Every admitted state retains its exact original256 physical records and the
same modeled O/H/H at indices256/257/258, IDs W/1//O,H1,H2. Metal index255,
whole259 charge−2, physical multiplicities La1/Dy6. Source-map, graph, original
termini and peptide protonation remain pinned. Preserve the archived water
template coordinates; applying TIP3P charges does not imply TIP3P optimization.

## Energy and charge model

E = saved native MACE-OMOL vacuum(whole259) + native fixed-FF ALPB(whole259).
The whole peptide, target and one water are Q. External-field and MM ledger
terms are partition-defined zero. Do not add water/peptide MM vacuum terms or
any old self-consistent GFN solvent offsets. No new MACE, GFN SCC, MM or DFT.

Reuse the exact original255 ff19SB peptide charges plus formal+3 target from
matched static12715. Append O−0.834,H1+0.417,H2+0.417, read from the existing
Mex-pinned OpenMM `tip3p.xml`, SHA256
`607f0fc9566c3770db2d9eb579fed68c4157578445eed1ec0aa7dccc23d57a6c`.
The original Mex native atom export independently records this same triplet.
No charge fitting, normalization, population projection or coordinate repair.
Water charge0 and full charge−2; serialize all charges with round-trip precision.

Reuse `LBT_MATCHED_MODEL_V2.inputs` for original FF/terminal/model/native-build
admission. The native wrapper differs from the qualified12715 source only in
the declared atom-count assertion256→259. Compile with the same pinned compiler,
static libraries and runtime libraries. Keep ALPB water/gsolv,230 grid points,
298.15K,zero ionic strength,p16 and exact declared-FF charge-rounding bookkeeping
dq0 unchanged. This is not a claim of zero physical charge-model uncertainty.

Before the first call, verify all twelve raw MACE/CELL/TASK/XYZ/physical/map/graph
pins, energy scalar replay, parent model/checkpoint/versions and La1/Dy6 state.
The saved MACE receipts are energy-only and contain no XYZ field: bind the exact
adjacent TASK.json to the parent row.task, then its XYZ/physical-source pins.
Require all original256 physical records exactly unchanged, paired-metal
coordinate equality, O/H/H identity, and original archived geometry-admission
flags. Do not rerun placement or geometry optimization. No256-atom energy or
solvent receipt substitutes for a259-atom state.

## Calls, availability and accounting

Exactly twelve new259-atom native fixed-FF ALPB attempts at most, one per admitted
slot in the frozen order, one CPU and zeroGPU. Archive native gradients from
the unchanged wrapper, but use energies only. Each slot gets at most one attempt;
failure remains unavailable without retry or substitution. Preserve all sixteen
rows and mark every dependent comparison unavailable if any required cell fails.
Full completion means twelve scored plus four explicitly rejected slots, never
sixteen available cells. Record actual attempts, build/runtime pins, raw receipts,
application timing and scheduler cost. Root submits separately after review.

## Frozen comparisons

For each source and each plane0/90, retain the paired inner-minus-outer works and
H = (E_Dy,inner − E_Dy,outer) − (E_La,inner − E_La,outer).
Compute balanced component differences in native eV/hartree before conversion
using the unchanged23.06054783061903 and627.509474 constants. Negative H means
the inner transfer is relatively more Dy-favorable. Retain the earlier1kcal/mol
interpretive resolution; it is not a statistical uncertainty or classifier band.
Report MACE and fixed-FF solvent components, total H, and the old saved H as a
separately named comparison. Never insert the old solvent into the new score.

Both7CCO plane comparisons can be available. Both7CCN transfers remain
unavailable because their inner constituents were rejected; report all its
outer raw/component rows. Retain the original complete four-pose menu per source,
available only if all eight source cells exist. Both metals receive the same
menu. Preserve old selections alongside the new result. No partial favorable
selector, global wet pool, source vote, or imputed missing cell.

The outer position is not bulk water. All transfers retain one explicit water;
there is no water chemical potential, hydration population, binding affinity or
raw minimum across dry/wet inventories. These fixed original sources do not
hydrate the later accommodated peptide pool or establish a robust Mex–LBT
preference. Retain terminal-charge, source-conditioned NMR, spin/frozen-f and
fixed-charge limitations. No further campaign follows automatically.

CLI: `python LBT_MATCHED_WATER.py --config CONFIG --output FRESH`.
Config pins implementation,plan,matched_helper,matched_origin_result,
water_preparation_result,water_scout_result,water_panel_result,tip3p_xml,
Mex water-charge provenance. No execution is authorized for the child author.
