# Matched-objective LBT response: finite experiment complete

**Matching the search objective changes the 7CCN source materially but does not
resolve source dependence.** The 7CCO and global paired contrasts move only
**+0.03675883 kcal/mol**, while 7CCN moves **+4.17796159 kcal/mol**. The absolute
source gap remains **24.12099058 kcal/mol**, compared with the previously saved
28.26219334 kcal/mol. These are conditional finite-menu responses, not affinities
or an accuracy improvement.

All four prescribed starts were retained. Jobs **12762/12765** completed their
own-state searches; **12767** completed all four opposite-state scores. The final
pool contains **22 geometries and 44 available metal cells**. The preceding
18 geometries, 36 cells and 15 menu columns remain unchanged and separately named.

## Same model and geometry rules

All four starts are the exact old ten-mode endpoints. Their new starting gradients
were reconstructed without molecular calls from raw MACE forces plus archived
12720 fixed-FF ALPB gradients. The whole 256-atom peptide, charge −2, physical
multiplicities La1/Dy6, native checkpoint, FF charge vector, source reference and
ten coordinates remain unchanged. The original 0.40 Å source, 0.20 Å target and
0.20-radian torsion envelopes, 0.05 Å trial limit and prescribed fallback remain.
Own-energy admission uses the exact matched warm start; work from the original
source is also archived. No old SCC-solvent composite gradient was reused.

## Four own-state searches

Values below are copied from the completed search receipts. Each molecular
request comprises one MACE and one fixed-charge ALPB evaluation.

| Start | Accepted steps / fallback | New requests | MACE work | ALPB work | Composite work from warm start |
|---|---:|---:|---:|---:|---:|
| 7CCO/La | 8 / 7 | 12 | +0.128323 | −1.006573 | **−0.878250** |
| 7CCO/Dy | 5 / 4 | 8 | +0.203641 | −1.045132 | **−0.841491** |
| 7CCN/La | 11 / 11 | 15 | −4.518667 | +0.040375 | **−4.478292** |
| 7CCN/Dy | 5 / 5 | 9 | +0.124089 | −0.424419 | **−0.300330** |

Work units are kcal/mol. Three incremental decreases are solvent-driven; the
larger 7CCN/La decrease comes from MACE vacuum energy.

Every search ends at `no_admissible_decreasing_trial`, preserving its last accepted
point. At each final iteration, all four normal trials fail the stated bounds;
then both allowed, geometrically admitted fallback trials fail the energy test:

| Start | Final normal-trial obstruction | Last fallback coordinates | Actual fallback composite changes, kcal/mol |
|---|---|---|---:|
| 7CCO/La | Torsion bound | target x; A10 N–CA | +0.215660; +0.015116 |
| 7CCO/Dy | Source bound; largest trial also target bound | target y; target x | +0.094622; +0.022068 |
| 7CCN/La | Source bound | target x; target z | +0.035155; +0.540043 |
| 7CCN/Dy | Source and torsion bounds | target x; target z | +0.099623; +0.259572 |

The four retained endpoints pass covalent, chirality, protected-group and all-atom
contact checks. None has a new/worsened severe contact. Near-bound positions and
nonzero force remain:

| Start | Source extent, Å | Target extent, Å | Maximum torsion, rad | Active-span maximum force: start → terminal, kcal/mol/Å |
|---|---:|---:|---:|---:|
| 7CCO/La | 0.299969 | 0.160214 | 0.199915 | 12.7572 → 2.4763 |
| 7CCO/Dy | 0.399931 | 0.180537 | 0.191205 | 13.6593 → 7.3705 |
| 7CCN/La | 0.399944 | 0.141658 | 0.175249 | 25.0831 → 22.6571 |
| 7CCN/Dy | 0.399944 | 0.167206 | 0.199624 | 3.7747 → 5.5306 |

All low-force gates remain false. A decreasing energy need not decrease the
projected force, as 7CCN/Dy illustrates. These stops are not stationary minima or
proof that every smaller feasible move is exhausted.

## Shared paired pool

Both metals see the same menu within each scope. Rejected trials do not become
candidates. Every newly selected endpoint is its corresponding own-state terminal.

| Scope | Selected La | Selected Dy | Change in paired Dy−La from old union, kcal/mol |
|---|---|---|---:|
| 7CCO | `geometry_18`: new 7CCO/La | `geometry_19`: new 7CCO/Dy | **+0.03675883** |
| 7CCN | `geometry_20`: new 7CCN/La | `geometry_21`: new 7CCN/Dy | **+4.17796159** |
| Global | `geometry_18`: new 7CCO/La | `geometry_19`: new 7CCO/Dy | **+0.03675883** |

The source gap remains large despite some reduction. A global menu choosing 7CCO
for both metals does not establish source robustness. The 7CCN shift exceeds the
predeclared 1 kcal/mol interpretive resolution and establishes material sampling
sensitivity to the objective. This finite experiment corrects the objective used
for these proposals; it does not establish improved biological accuracy. The
balanced comparison with fixed Mex common5 is a separate saved-data report.
No automatic extension of the bounds, modes or optimization follows from this result.

## Actual cost and evidence

| Job | Purpose | Scheduler seconds | CPUs / GPUs | Allocated CPU-s / GPU-s | New MACE / ALPB |
|---|---|---:|---:|---:|---:|
| 12762 | Scout | 27 | 12 / 1 | 324 / 27 | 12 / 12 |
| 12765 | Remaining three | 48 | 12 / 1 | 576 / 48 | 32 / 32 |
| 12767 | Opposite-state pool | 17 | 12 / 1 | 204 / 17 | 4 / 4 |
| Total | This finite experiment | — | — | **1,104 / 92** | **48 / 48** |

No new GFN2, MM or DFT calls. These are allocation costs, not utilization, and
exclude earlier model/force qualification and later saved comparison jobs.
[FINAL_ACCOUNTING.json](FINAL_ACCOUNTING.json) pins all implementations,
configurations, actual results, searches, selected points and completion snapshots;
[FINAL_SACCT.psv](FINAL_SACCT.psv) retains the raw scheduler records.

Actual result SHA256 prefixes: scout `a47736fe`, remaining `e6ad83e9`, pool
`9ab5a556`. The prior source gap is read from the saved
[12736 comparison](../../../workspaces/lanm_coupled_response_20261016/monomer_continued_comparison_v1/execution/RESULT.json),
not recomputed here. This report reads saved scientific outputs and scheduler
receipts only; no scientific calculation, geometry change or submission was run.
