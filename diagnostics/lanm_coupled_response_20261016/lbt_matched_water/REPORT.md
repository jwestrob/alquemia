# Fixed-water response under matched FF solvent

**The negative numerical sign persists, but the earlier roughly −2.5 kcal/mol response does not persist at its prior magnitude.** Both matched-FF contrasts are within the predeclared 1 kcal/mol near-zero range. The result shows substantial solvent-model dependence for this fixed source and water-placement question.

All 16 declared slots are preserved: 12 scored successfully and the same four 7CCN inner slots remain geometry-rejected. No new geometry, water occupancy, MACE or GFN calculation was introduced.

## Paired water-transfer response

`H = (E_Dy,inner − E_Dy,outer) − (E_La,inner − E_La,outer)`. Negative H means the inner transfer is relatively more Dy-favorable; it does not by itself mean the inner position lowers either absolute energy. All entries are kcal/mol.

| Source / plane | MACE H | Old self-consistent solvent H | New fixed-FF solvent H | Old total H | Matched-FF total H | Frozen interpretation |
|---|---:|---:|---:|---:|---:|---|
| 7CCO / 0° | -2.108176 | -0.434413 | +1.321202 | -2.542589 | -0.786974 | near zero |
| 7CCO / 90° | -2.099925 | -0.375789 | +1.360345 | -2.475714 | -0.739580 | near zero |
| 7CCN / 0° | — | — | — | unavailable | unavailable | inner cells rejected |
| 7CCN / 90° | — | — | — | unavailable | unavailable | inner cells rejected |

The saved vacuum differential stays near −2.1. The solvent contribution changes sign and opposes that response under fixed FF charges. The old resolved negative result is therefore not robust to this declared solvent-model change. The 1 kcal/mol range is an interpretive resolution, not an error bar or evidence of an exactly zero effect.

| 7CCO plane | La inner−outer | Dy inner−outer |
|---|---:|---:|
| 0° | +8.002234 | +7.215260 |
| 90° | -0.798329 | -1.537909 |

Plane0 transfer is uphill for both metals; plane90 transfer is downhill for both. In the complete four-pose 7CCO menu, both metals select inner_plane90, unchanged from the old model and without exact ties. The complete 7CCN menu remains unavailable. No partial outer-only selector or global wet pool is substituted.

## Complete sixteen-slot record

Native energy components below are copied from the saved cells; elemental absolute offsets are not affinity measures. The four rejected rows retain their original status.

| Source | Pose | Metal | Status | Saved MACE (eV) | New fixed-FF ALPB (hartree) |
|---|---|---|---|---:|---:|
| 7CCO | inner_plane0 | La | complete | -186315.834765647 | -0.529407776382 |
| 7CCO | inner_plane0 | Dy | complete | -209539.117139605 | -0.543127377053 |
| 7CCO | inner_plane90 | La | complete | -186316.516233774 | -0.525277876797 |
| 7CCO | inner_plane90 | Dy | complete | -209539.798509505 | -0.538927457288 |
| 7CCO | outer_plane0 | La | complete | -186315.952982097 | -0.537815773780 |
| 7CCO | outer_plane0 | Dy | complete | -209539.143936886 | -0.553640843275 |
| 7CCO | outer_plane90 | La | complete | -186316.170498466 | -0.536711197355 |
| 7CCO | outer_plane90 | Dy | complete | -209539.361712829 | -0.552528625532 |
| 7CCN | inner_plane0 | La | geometry_rejected | — | — |
| 7CCN | inner_plane0 | Dy | geometry_rejected | — | — |
| 7CCN | inner_plane90 | La | geometry_rejected | — | — |
| 7CCN | inner_plane90 | Dy | geometry_rejected | — | — |
| 7CCN | outer_plane0 | La | complete | -186313.118291302 | -0.476455189963 |
| 7CCN | outer_plane0 | Dy | complete | -209539.443749128 | -0.482573228006 |
| 7CCN | outer_plane90 | La | complete | -186313.124191901 | -0.476983500832 |
| 7CCN | outer_plane90 | Dy | complete | -209539.447942838 | -0.483080148653 |

Full precision, the original old-model summary, all source/task/physical/graph pins and every native receipt remain in [the execution result](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_coupled_response_20261016/lbt_matched_water_v1/execution/RESULT.json). The original twelve-cell energy panel and four preparation exclusions are unchanged.

## Scope and limits

The whole259-atom peptide/metal/water system is the electronic region. Its energy is archived MACE vacuum plus new native fixed-FF ALPB; external-field and MM partition terms are zero. The solvent vector is the existing255 peptide FF charges, formal+3 target, and TIP3P O−0.834/H+0.417/H+0.417. Geometry, charge−2, La1/Dy6 state, charged termini and original source maps are unchanged. The native functional/build lineage is unchanged except its declared259-atom assertion.

This is one modeled water on the original fixed sources. It neither samples hydration occupancy nor tests water on the later accommodated peptide pool. Outer water is not bulk water; no dry/wet raw minimum, missing-cell imputation, source-robust preference, protein Kd or affinity is inferred. The 7CCN omissions prevent a complete second-source transfer test. Terminal-charge, source-conditioned NMR, spin/frozen-f and fixed-charge model limitations remain.

## Calls and cost

Slurm12732 COMPLETED, exit0:0: **12 native ALPB calls, 9 scheduler CPU-seconds, zero GPU**. All12 MACE cells were reused. No new GFN SCC, MM, DFT, geometry or optimization. Application duration was7.882003559963778s, including admission/build; prior preparation/scoring costs are excluded. [ACCOUNTING.json](ACCOUNTING.json) preserves the terminal record and exact counters.

- Result SHA256: `8d9a0ea1e8876a39b41e5858eba7d064b6b5b2236a6b09e71a2ce3338cba43ee`.
- Source/config/plan: `2b1c2ae2` / `aaa94fdf` / `6c92404c`.
- Original complete panel12670: `aa4dc92c`; reused scout12664: `26aa78c2`; preparation12565: `43799549`.

This report transcribes and interprets completed saved results and performs accounting arithmetic only. No new scientific computation, submission or follow-on campaign is authorized by this report.
