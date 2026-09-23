# Nikasha: parallel search and basin experiments

**All three approved pilots are complete. None establishes additional classifier
improvement over the existing adaptive scorer.** The basin calculation works on
two of four curves, but its breadth contribution is small and electronic-solver
continuity remains a concrete concern. Production, its reference bands, explicit
DFT access and all historical outputs remain unchanged. The earlier adaptive improvement remains
separate: 200 to 202 correct on 204 matched structural sources, with coverage and
cost limitations in the [previous report](../nikasha_next_phase_20260922/REPORT.md).

## What these experiments test

Jacob approved three parallel experiments in [PLAN.md](PLAN.md). They target
different explanations for residual structural sensitivity: an optimizer starting
in the wrong basin, missing geometries favored by solvent, and differing basin
breadth at comparable minimum energy. The common eight sources and four-source
integration subset are pinned in [INPUTS.json](INPUTS.json). These are selected,
previously inspected development examples, not independent biological validation.

All three use the same source atoms, chemical states, source-connected physical
motions and existing native float64 MACE plus native GFN2 solvent transfer:

`E_M(q) = E_OMOL,vac,M(q) + E_GFN2,ALPB,M(q) - E_GFN2,vac,M(q)`

`R = min_q E_Ca(q) - min_q E_La(q)` over one shared finite geometry pool.

Existing operational selection tolerance and exact mathematical minima remain
separate. A changed search has no newly validated threshold. Both released-band
and adaptive-band transfers are reported without recalibration or padded edges.
No unavailable cell is replaced with a successful baseline result.

## Completed findings

| Branch | Actual result | Implication |
|---|---|---|
| Reference-derived alternative starts | Four of 16 starts fit the fixed domain; all eight resulting metal searches return the old basins. Four difficult folds admit neither template. | This bounded starting-point rule adds no useful solution and does not test a wider collective rearrangement. |
| Solvent-guided finite probes | All eight scores complete; seven correct and one wrong under transferred adaptive bands, unchanged. | Solvent favors some new coordinates, but the tested local rule repairs no decision and does not reduce the selected structural spread. |
| Direct donor-torsion basin breadth | Two of four curves pass the frozen integration gates. Width changes the contrast by only about 0.08–0.11 kcal/mol there. The expanded finite pool retains the old minima on all three fully scored sources. | No additional classifier gain. Keep the conditional integration primitive, but omit this grid from routine scoring. |

The alternative-start Q9Z4J7 result crosses an exact calibration boundary by
0.00000377 model kcal/mol. Its literal inconclusive status is retained; this is
not a meaningful loss of chemical discrimination. A separate, much larger Q88JH5
change has a numerical concern, described below. Full denominators, individual
works, selected states and both band transfers are in the shared comparison;
see the [start report](../structure_informed_starts_20260922/REPORT.md),
[solvent-probe report](../solvent_guided_20260922/REPORT.md), and
[complete shared comparison](../../workspaces/nikasha_parallel_pilots_20260922/COMPARISON_v1.md).

## What the basin experiment adds

| Source | Integration status | Minimum contribution to ΔR | Width contribution to ΔR |
|---|---|---:|---:|
| 1H4I | Unavailable: two La vacuum points fail native SCF | unavailable | unavailable |
| 4MAE | Frozen integration gates pass | −0.389541 | −0.106738 |
| Q9Z4J7 | Fails quadrature and domain checks | not qualified | not qualified |
| A0ACD6B9F2 Ca-conditioned sample 4 | Frozen integration gates pass | +16.161591 | −0.083615 |

All contributions are model kcal/mol, relative to each metal's own original
geometry on the same one-coordinate curve. For the large A0ACD6B9F2 response,
the minimum shift dominates; breadth nearly cancels between metals. That supports
prioritizing accommodation energy over a routine entropy calculation for this
tested coordinate. It does not establish that entropy is negligible for other
motions, chemical states or water inventories.

Q9Z4J7 has a 0.645 kcal/mol paired quadrature discrepancy and 68.7% of the La weight
outside the fixed inner domain. No post-result extension or interpolation was
used. The two 1H4I failures reached 500 SCF cycles; both remain visible. In the
expanded common pool, 4MAE, Q9Z4J7 and A0ACD6B9F2 select exactly their existing
adaptive winners, with zero contrast change. 1H4I's new pool is unavailable;
its previous score stays separately available.

Passing integration gates establishes conditional quadrature/domain convergence,
not a unique or accurate electronic surface. The unsmoothed solvent component
has large adjacent-point jumps on Q9Z4J7 and A0ACD6B9F2 (about 5.60 and 24.71 kcal/mol
for La). The latter lies outside the strongly weighted well. These finite-step
jumps are retained without attributing a unique cause; they do not justify a
physical entropy or predictive-improvement claim. See the
[basin report](../local_basin_breadth_20260922/REPORT.md) for full curves, gates,
component accounting and editable figures.

## A concrete numerical limitation

Two Q88JH5 La geometries differ by only 6.46e-7 Å and have indistinguishable MACE
energy, yet native vacuum GFN2 differs by 4.805 kcal/mol. ALPB nearly agrees, so
the subtraction carries that jump into the composite. Eight fresh exact-input
repetitions reproduce both outcomes, excluding host/scheduling noise at this
scale. The numerical discontinuity is established; its unique cause is not.

This limits claims about small geometry improvements in the affected calculation.
It does not demonstrate that every static prediction is wrong. Generic printed
SCF density thresholds do not certify native-mixer convergence; a separate native
accuracy control was not found in the documented interface. Actual evidence and
the installed standalone-xTB fallback are recorded in the
[numerical-control note](../solvent_guided_20260922/NATIVE_ACCURACY_NOTE.md).
The fallback remains untested and is not mixed into the current solvent transfer.

Saved native charge/multipole states support an eight-call self/cross-restart
diagnostic, [proposed here](../structure_informed_starts_20260922/RESTART_PROPOSED.md).
It has not run. No tightening, repair or threshold change is implied by normal
termination of the completed calls.

## Implementation and accounting

`nikasha_finite_candidates.py` admits declared source-connected geometries,
checks states and mappings, deduplicates exact coordinates, and reuses existing
pool runners and compatible cells. `nikasha_parallel_compare.py` retains static,
adaptive and challenger fields side by side. The new finite-pool protocol is
`nikasha_declared_finite_geometry_native_OMOL_GFN2_ALPB_v1`; the separate integral
protocol is `physical_single_rotor_conditional_breadth_v1`.

The basin integral uses only its 65 prescribed points, with fixed spectators and
the common physical measure dq/(2π) at 300 K. Its minimum and width contributions
belong to that same curve. The shared scorer additionally retains old adaptive
pool candidates; its expanded-pool minimum is a different quantity. No width
term is added to an unrelated adaptive minimum. This experiment cannot establish
whole-pocket entropy, water occupancy, binding free energy or an equilibrium
population of chemical states.

| Branch | New MACE calls | GFN2 attempts | Allocated core-seconds | Requested GPU-seconds |
|---|---:|---:|---:|---:|
| Alternative starts, including numerical repeats | 295 | 40 | 11,712 | 124 |
| Solvent-guided probes | 54 | 108 | 18,208 | 35 |
| Basin grid | 512 | 1,024 | 235,552 | 221 |
| **Total** | **861** | **1,172** | **265,472** | **380** |

All MACE calls succeeded; 1,170 GFN2 attempts converged and two failed. There were
zero new DFT calls. Allocation totals include worker startup, failed calculations
and collection time. Four zero-runtime pending basin allocations were replaced
after a verified inherited node-exclusivity mismatch; no chemistry was duplicated.
The solver jobs ran eight 8-rank tasks concurrently. Their startup/collection
overhead was 86–138 allocated seconds per shard; MACE startup/collection outside
the worker phase was about 132 seconds. Local preparation, tests and reporting add
unmetered work; initial basin geometry preparation alone took 264.04 wall seconds.
These are measured development costs, not source-inclusive production latency
or active CPU utilization. The dense grid has not earned a place in the scanner.

Shared tests: five admission tests, three comparison tests and 30 existing pool
regressions passed on real pinned artifacts and explicitly corrupted copies.
The branch suites pass 12 alternative-start, 9 solvent-probe and 9 basin tests. The
basin JSON writer required a boolean serialization fix; its partial v1 is retained
and the successful algebra/report is v2, with no molecular recalculation.
These implementation tests are distinct from the executed molecular calls. No new
DFT, folds, MD, model training, production rescore or default change occurred.

## Next decision

Close the two tested search rules and omit routine basin-width integration.
The practical numerical target is the demonstrated solvent-energy discontinuity:
the proposed eight-call restart check can test electronic initialization without
launching another broad campaign. A higher-upside geometry question remains
collective scaffold response: can neighboring residues move together so the same donor displacement is useful?
The [matched-target proposal](../collective_scaffold_20260922/PROPOSED.md) separates
that effect from ordinary force-field relaxation. It is proposed-only and cannot
be promoted as a complete scaffold-energy correction. Its implementation can
progress alongside the targeted numerical investigation; small composite-energy
improvements need that continuity check before they can justify model selection.
The restart, standalone-xTB fallback, scaffold and redox proposals remain
unexecuted. No production promotion or remote push occurred in this phase.

Exact report-only replay commands are in [COMMANDS.md](COMMANDS.md). Machine-readable
outcomes, component/result pins and measured costs are in [RESULT.json](RESULT.json).
