# Physical response in the six consumed PLM sources

Completed 2026-09-23, existing artifacts only. **The score-spread reduction accompanies relief of the previously identified short extra-Asp contacts, predominantly through greater native-MACE stabilization of La.** This is a descriptive association in two unlabeled proteins, not a causal energy decomposition or validation of composite energetics/biological accuracy. No new molecular calculations, proposals, selection rules or scoring thresholds were used.

## Actual source and geometry accounting

The authoritative six-source execution is `workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json` (SHA256 `06e2c0b7727ff0aa53e6085b144a307c2eb6ca71cd0cc8f06112e3b8ec9e4be7`). All declared seed-101 samples 0/1/2 are retained for each protein; all were La-conditioned structures. Their metal-use labels remain unknown. The 07ab and 8344 abbreviations below mean `PQQSEQ_07ab500e3df76b30d71c` and `PQQSEQ_83440678cbbd658047c9`.

All twelve selected geometries are the corresponding metal's own four-angle proposal, confirmed from the actual **composite** common-pool matrix, not inferred from lower vacuum energy. Mathematical and operational selections coincide. The 213-atom (07ab) and 201-atom (8344) scored contexts retain charge −3 Ca/−2 La, unchanged deprotonated acidic roles and the archived dry inventory. PQQ and metal coordinates remain fixed. Source-heavy atoms move through intact bonded side-chain rotations; cap positions follow their existing physical map, with no independent cap coordinates.

The role maps identify Glu203/Asn287/catalytic Asp329/extra Asp331 in 07ab and Glu193/Asn263/catalytic Asp305/extra Asp307 in 8344. Anchor Glu and extra Asp move in all six sources. Asn moves only in 07ab sample2 and 8344 sample1. Catalytic Asp, PQQ and the remaining scaffold are fixed. Exact active mode IDs, atom identities, full angular coordinates, source/core/context pins and proposal receipts are in the pinned `GEOMETRY.json`.

## Distances and constraints

Extra-Asp distances are OD2 / OD1 in Å, retaining source atom identity. Contact counts use the existing chemically typed 3.1 Å convention; they describe geometry, not experimentally established bonds. The older 2.2 Å line is a La-structure warning, not a new threshold or a declaration that every shorter Ca contact is defective.

| Source | Original | Selected Ca | Selected La | Maximum heavy displacement Ca / La (Å) |
|---|---:|---:|---:|---:|
| 07ab sample0 | 1.847 / 3.714 | 2.282 / 4.142 | 2.350 / 4.214 | 0.683 / **0.800** |
| 07ab sample1 | 1.735 / 3.675 | 2.219 / 4.331 | 2.219 / 4.331 | **0.800 / 0.800** |
| 07ab sample2 | 2.120 / 3.604 | 2.168 / 3.559 | 2.317 / 3.442 | 0.521 / 0.690 |
| 8344 sample0 | 1.766 / 3.538 | 2.411 / 3.362 | 2.451 / 3.304 | 0.761 / **0.800** |
| 8344 sample1 | 2.055 / 3.248 | 2.028 / 3.263 | 2.451 / 2.974 | 0.301 / **0.800** |
| 8344 sample2 | 1.700 / 3.562 | 2.304 / 3.346 | 2.304 / 3.321 | **0.800 / 0.800** |

All selected La short-O distances rise above 2.2 Å. Two selected Ca contacts remain shorter: 07ab sample2 (2.168 Å) and 8344 sample1 (2.028 Å); the latter shortens slightly. Seven of twelve selected geometries reach the declared heavy-displacement boundary (five La, two Ca). These are bounded candidates, not unconstrained minima.

Extra Asp remains geometrically monodentate in five sources; only 8344 sample1 La gains a second contact within 3.1 Å. Anchor Glu in 07ab stays one-contact: its selected La second-O distances are 3.134–3.149 Å, just outside the existing cutoff. Anchor Glu in 8344 was already two-contact and remains so, with more equal La–O distances. Thus the common effect is relief/reorganization of a short contact, **not a universal monodentate-to-bidentate switch**. [DISTANCES.csv](DISTANCES.csv) retains all 126 role-atom observations, including actual source IDs and per-atom displacement.

## Actual energy response

Work is each selected candidate minus its own source origin. The score shift is `ΔR = work(Ca) − work(La)`; positive values here reflect greater La stabilization. The component table is in model kcal/mol, with no free-energy/affinity interpretation.

| Source | Native-MACE differential work | ALPB−vacuum differential work | Composite ΔR |
|---|---:|---:|---:|
| 07ab sample0 | +38.384 | −1.528 | +36.855 |
| 07ab sample1 | +66.918 | −3.936 | +62.983 |
| 07ab sample2 | +18.912 | +2.652 | +21.564 |
| 8344 sample0 | +66.168 | −2.692 | +63.475 |
| 8344 sample1 | +28.155 | −4.569 | +23.587 |
| 8344 sample2 | +82.414 | −3.641 | +78.773 |

The most compressed source in each protein receives the largest positive contrast shift. Within the **same envelope representation**, the three-source ranges fall from 40.138 to 2.307 (07ab) and 72.416 to 17.230 (8344) kcal/mol. Large source-dependent native-MACE work dominates this compression of the ranges; solvent subtraction modifies it, usually opposing it. Exact per-metal work and unrounded original/selected scores are in [COMPONENT_WORKS.csv](COMPONENT_WORKS.csv).

This analysis does not isolate the energetic contribution of extra Asp from simultaneous Glu/Asn motion. It also does not isolate the effect of enlarging the environmental context: the before/after comparison holds this context fixed. The compression was identified in earlier PLM geometry work; this report traces how the completed candidates respond to it. Unknown labels, boundary-limited motions and residual 8344 source spread prevent treating this as established improved biological accuracy.

## Figure and verification

Editable [SVG](../../workspaces/plm_envelope_response_20260923/figures_v2/plm_geometry_response.svg) and [PDF](../../workspaces/plm_envelope_response_20260923/figures_v2/plm_geometry_response.pdf) show both extra-Asp oxygens, all six sources, selected boundaries and the actual differential components. Connecting lines are visual guides between discrete states, not a sampled transition path. Left panels show distances; right panels show candidate-minus-origin Ca−La work. Open squares flag the 0.8 Å boundary. All experimental labels are unknown. Figures v1 are preserved; v2 only repositions the legend to uncover data.

Read-only verification replayed all 18 geometries against both metal score cells (36 XYZ comparisons, tolerance 1e−12 Å), confirmed identical paired maps and charges, fixed PQQ/metal coordinates, retained acidic states and exact component/score algebra (1e−7 kcal/mol check). CSV export retains all 126 atom rows and six response rows. `ARTIFACTS.json` pins the result, full geometry, CSVs, script and figures. The original analysis-script bytes are preserved beside `GEOMETRY.json`; subsequent source change only moves the figure legend. No new molecular CPU/GPU cost or scientific tests were introduced.
