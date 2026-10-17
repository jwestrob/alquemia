# Mex seven-coordinate continuation and complete common5 comparison

**Complete: both own searches and all five geometries / ten paired cells are available.** Each metal selects its own terminal geometry. Dy gains more accommodation energy: the saved Dy−La accommodation change is −1.490527 kcal/mol from the original source, versus −0.102379 for the original common3. The additional common3→common5 change is −1.388149 kcal/mol. This is not evidence of improved discrimination.

Both searches accepted 12 steps from their exact pilot starts and stopped at the predeclared accepted-step limit. They remain nonstationary. The common pool adds exactly two opposite-metal evaluations and preserves the original common3 result.

## Complete shared menu

Work is relative to each metal at the same original physical source (kcal/mol). Both metals receive all five geometries; no failed or unfavorable slot is dropped.

| Geometry | La work | Dy work | Dy−La work |
|---|---:|---:|---:|
| origin | +0.000000 | +0.000000 | +0.000000 |
| La_pilot | -0.252489 | -0.260656 | -0.008167 |
| Dy_pilot | -0.233258 | -0.354867 | -0.121609 |
| La_terminal | -3.588288 | -1.665229 | +1.923059 |
| Dy_terminal | -0.672681 | -5.078816 | -4.406135 |

The La-terminal geometry is relatively La-directed (+1.923059 differential work), whereas the Dy-terminal geometry is relatively Dy-directed (−4.406135). The common5 selector allows each metal to choose independently from the same menu and selects La_terminal/Dy_terminal. Its common selected work is −4.333552 kcal/mol. All five physical geometries are distinct; declared aliases and the ten raw paired cells remain in the execution result.

## Own response and remaining forces

| Metal | Work beyond pilot | Work from original | Active-span maximum atom gradient | Active-span Cartesian RMS | Source displacement | Target displacement | Max total torsion |
|---|---:|---:|---:|---:|---:|---:|---:|
| La | -3.335800 | -3.588288 | 14.056838 | 0.194325 | 0.181529 | 0.102301 | 0.077454 |
| Dy | -4.723948 | -5.078816 | 5.122077 | 0.087399 | 0.198122 | 0.081955 | 0.085673 |

Energy units are kcal/mol; gradient units kcal/mol/Å; displacement Å; torsion radians. Both active metrics remain above the frozen 0.1 maximum / 0.03 RMS stopping thresholds. The saved force outside the seven-coordinate span has L2 norm 1141.255790 (La) and 1143.673083 (Dy) kcal/mol/Å; the limited search does not establish full Cartesian stationarity.

All original-source guards pass: fixed complement exact, no newly worsened severe contacts, no recorded chirality or peptide-basin failures, and protected bond/angle checks preserved. Both metrics retain rank seven. Source/target envelopes are 0.20Å and total torsion limit is 0.20rad; Dy approaches the source envelope at 0.198122Å, but the recorded termination is the step limit, not a bound or convergence. Final self-check fields for incremental movement are not interpreted as the largest step taken.

## Own terminal component work from original source

| Metal | Local MACE vacuum | GFN field difference | MM ledger | Fixed-FF ALPB | Total |
|---|---:|---:|---:|---:|---:|
| La | -1.740423 | +0.246059 | -2.440617 | +0.346693 | -3.588288 |
| Dy | -4.630797 | +0.885451 | -2.734612 | +1.401142 | -5.078816 |

Vacuum and MM terms favor each own endpoint; field response and fixed-charge solvent oppose part of that gain. The paired cross scores are included above, so the conclusion does not rely on own-only energy decreases.

## Scope and limits

The Hamiltonian is local MACE vacuum + local GFN embedded-minus-isolated response + the existing MM ledger + global native fixed-FF ALPB. The physical system has 2049 atoms, a 209-atom capped local region and 1840 shifted field sites. Local physical multiplicities remain La1/Dy6; GFN uses the separately declared effective UHF0/frozen-f model. The seven coordinates are target xyz, Leu83 N–CA hinge and Asp84/Asp86/Glu95 CA–CB torsions.

This is conditional EF3 response in one Mex source with three fixed Nd spectators. It does not establish equilibrium affinity, protein Kd, source robustness or a stationary minimum. Comparison against LBT menus belongs to the separately declared saved-data comparator; this report does not combine those scores or authorize more search.

## Calls and actual costs

| Job | Work | Composite / MACE | Local GFN2 | MM configurations / groups | ALPB | Scheduler seconds | CPU-s | GPU-s |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 12727 | La continuation | 12 / 12 | 24 | 12 / 48 | 12 | 121 | 1452 | 121 |
| 12729 | Dy continuation | 12 / 12 | 24 | 12 / 48 | 12 | 121 | 1452 | 121 |
| 12731 | Two terminal cross cells | 2 / 2 | 4 | 2 / 8 | 2 | 32 | 384 | 32 |

Incremental three-job allocation cost: **3288 CPU-seconds and 274 GPU-seconds**. Each used 12 CPUs and one H200, with one model load per stage. No DFT or whole-protein SCC. Earlier pilot/static/qualification costs are excluded. Application wall times are 119.218240s, 119.164630s and 30.473602s; these differ from the scheduler intervals above. [ACCOUNTING.json](ACCOUNTING.json) preserves raw timing/counter fields, terminal records and root-review pins. Allocation cost is not measured process utilization.

## Saved evidence

- Common5: `workspaces/lanm_coupled_response_20261016/mex_continuation_pool_v1/execution/RESULT.json`, SHA256 `c93f7df93c0c58bebd99d7f2281804c528c596cf30d8cb3ff6572a05a2ad746f`.
- La result/search: `mex_continuation_scout_v1`, result `3fac70c7`, search `c969f90a`.
- Dy result/search: `mex_continuation_Dy_v1`, result `99c1026a`, search `083be154`.
- Frozen scout/base `c571f20e`, remaining wrapper `dddc4d25`, common5 adapter `50cbd762`; plan `88dc7df6`. Root admission records exist for all three jobs.

This report transcribes saved scientific values and performs only accounting arithmetic. It adds no molecular call or scientific numerical analysis. Existing SCOUT_REPORT.md and SCOUT_ACCOUNTING.json remain unchanged.
