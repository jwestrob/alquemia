# Matched-FF rescore of the complete LBT candidate pool

**Complete: 18 geometries, 36 La/Dy cells, all 15 declared menu summaries available.** The earlier motion reduces the large 7CCN exchange magnitude, while the candidate/union global result stays close to its static value. Both sources retain the expected sign, but substantial source dependence remains.

## Comparison and full denominator

`R = (E_Mex,Dy − E_Mex,La) − (E_LBT,Dy − E_LBT,La)`, in kcal/mol. Positive R favors Mex–La + LBT–Dy. Mex remains the admitted static12685 baseline in every row; latest Mex12717 accommodation is reserved for the separate root-owned saved comparison.

Both metals select independently from each identical, complete menu. The old source-local prior/control/candidate menus are retained, with explicit origin and full-union rows. Control proposals used seven coordinates; candidate proposals used ten. Counts below are geometries (each has both metal states); indices abbreviate `geometry_XX`.

| Source scope | Menu | Geometries | Selected La / Dy | R (kcal/mol) | Change from origin |
|---|---|---:|---|---:|---:|
| 7CCO | origin | 1 | 00 / 00 | +2.532554 | +0.000000 |
| 7CCO | prior | 5 | 06 / 03 | +1.580410 | -0.952144 |
| 7CCO | control | 7 | 10 / 12 | +1.579011 | -0.953542 |
| 7CCO | candidate | 7 | 11 / 13 | +2.509905 | -0.022648 |
| 7CCO | union | 9 | 11 / 13 | +2.509905 | -0.022648 |
| 7CCN | origin | 1 | 01 / 01 | +68.101238 | +0.000000 |
| 7CCN | prior | 5 | 08 / 08 | +44.723749 | -23.377488 |
| 7CCN | control | 7 | 14 / 14 | +40.013871 | -28.087367 |
| 7CCN | candidate | 7 | 15 / 17 | +30.772099 | -37.329139 |
| 7CCN | union | 9 | 15 / 17 | +30.772099 | -37.329139 |
| global | origin | 2 | 00 / 00 | +2.532554 | +0.000000 |
| global | prior | 10 | 06 / 03 | +1.580410 | -0.952144 |
| global | control | 14 | 10 / 12 | +1.579011 | -0.953542 |
| global | candidate | 14 | 11 / 13 | +2.509905 | -0.022648 |
| global | union | 18 | 11 / 13 | +2.509905 | -0.022648 |

All selections have singleton exact-tie sets. Full precision, all 36 cell records, native component energies, exact candidate memberships, provenance, and original selections are preserved in [RESULT.json](RESULT.json), copied byte-for-byte from the execution result.

## What the motion changes

- The 7CCN exchange decreases from +68.101238 at the origin to +30.772099 with candidate/union selection (saved change −37.329139). This narrows the source disagreement without removing it.
- The 7CCO exchange moves from +2.532554 to +2.509905 (saved change −0.022648). The global menu makes those same selections and has that same result; it does not establish a larger separation after accommodation.
- The prior and control menus remain visible: their 7CCO/global exchanges are +1.580410 and +1.579011. The near-unchanged global result refers specifically to candidate/union selection.
- Candidate and union select the same endpoints in both source contexts. The union adds no selected minimum beyond the candidate menu. The candidate endpoints lower the selected own-metal energies relative to the control menu in this finite comparison.
- All 15 exchange signs are positive. There is no source sign reversal. A global-only summary would conceal the remaining difference in magnitude between the two source contexts.
- Under matched fixed-FF solvent, every global menu selects both metals from 7CCO-derived geometries. The old self-consistent-solvent global prior/control/candidate summaries selected Dy from 7CCN (09/14/17); the new selections are 7CCO (03/12/13). This is model-dependent reselection, not independent evidence that the deposited 7CCO source is the physical equilibrium state.

## Component evidence

Candidate and union have identical component exchanges. The global row duplicates 7CCO. Entries below are copied from the saved component results (kcal/mol).

| Source scope | MACE vacuum | GFN field difference | MM ledger | Fixed-FF ALPB | R |
|---|---:|---:|---:|---:|---:|
| 7CCO | +3.856884 | -0.041319 | +0.066701 | -1.372361 | +2.509905 |
| 7CCN | +32.531601 | -0.041319 | +0.066701 | -1.784884 | +30.772099 |
| global | +3.856884 | -0.041319 | +0.066701 | -1.372361 | +2.509905 |

Selected candidate/union own-metal component work from the corresponding origin is also saved directly. These describe the LBT geometry changes; its field and MM partition terms are zero.

| Source | Metal | MACE work | Fixed-FF ALPB work |
|---|---|---:|---:|
| 7CCO | La | -3.813829 | -0.121633 |
| 7CCO | Dy | -8.815869 | +4.903056 |
| 7CCN | La | -55.180827 | -1.606119 |
| 7CCN | Dy | -16.808325 | -2.649482 |

The 7CCO Dy vacuum improvement is partly opposed by its solvent change; the 7CCN La vacuum response is much larger than its Dy response. These saved components explain why improving own-metal energies does not imply increasing the balanced exchange.

## Model and interpretation limits

This is a rescore of a fixed consumed candidate set. Its geometries were generated under the earlier MACE plus self-consistent GFN solvent objective. No matched-FF optimization, new stationarity check, changed bound, retry, extra geometry, or favorable source selection was performed. Old solvent energies were not inserted into the new total. The finite common pool is not a population.

The matched energy construction is MACE vacuum + local GFN embedded-minus-isolated field response + MM ledger + global fixed-FF ALPB. For whole-Q LBT, field/MM terms are partition-defined zero; no whole-peptide FF vacuum term is added. Mex retains its capped local-Q partition, fixed EF3 target context and three Nd spectators. Different electronic partitions and dry LBT versus hydrated Mex remain limitations.

LBT retains the existing 255 ff19SB charges plus formal +3 target, free charged termini and holo charge −2; the unresolved primary-source terminal-charge convention remains explicit. La1/Dy6 MACE metadata, frozen-f/effective-spin qualifications, and the fixed-charge approximation remain. These consumed peptide/protein contexts do not establish a protein Kd, classifier, independent validation, or equilibrium relative affinity.

## Execution and provenance

- Slurm12720: COMPLETED, exit0:0; 13 scheduler seconds on one CPU, zero GPU. Application elapsed 12.1840522789862 seconds. Incremental scheduler cost: 13 CPU-seconds.
- 32 new native 256-atom fixed-FF ALPB calls, four exact12715 origin ALPB reuses, and all36 saved MACE reuses. Zero new MACE, GFN SCC, MM, DFT or optimizer work. Earlier search/static costs are excluded from this incremental receipt.
- Execution result SHA256: `ee867343321293d83c8845aee7ad6b7f5e6e4661a738f1c7c70c015ae1262928`.
- Implementation `c9856d48`; config `a330a901`; plan `0f752f7f`; inventory `aa81b392`. Native helper/build/charge/model admission and all raw saved-MACE evidence passed before calls. Original12560 outcomes remain in the result.
- [ACCOUNTING.json](ACCOUNTING.json) preserves the observed terminal scheduler record separately from application timing. [The frozen plan](../LBT_MATCHED_POOL_PLAN.md) defines the full denominator and limits.

This report transcribes and interprets completed Slurm outputs. No molecular calculation or new scientific numerical analysis was run while preparing it. Root owns integration and any separate comparison with Mex12717.
