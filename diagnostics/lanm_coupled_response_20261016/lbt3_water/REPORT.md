# One-water response: conditional Dy-directed effect in 7CCO

**The finite panel is complete: 12 of 16 cells are available.** In the fixed
7CCO peptide, moving the modeled water inward gives a Dy-minus-La response of
about −2.5 kcal/mol in both declared orientations. Both exceed the predeclared
1 kcal/mol interpretation resolution. All four 7CCN inner-water cells remain
geometry-rejected, so their transfer contrasts are unavailable. This does not
resolve source dependence or establish a binding preference.

## Same-composition result

Every cell contains the original fixed 256 atoms plus one intact neutral water:
259 atoms, total charge −2, physical La multiplicity 1 or Dy multiplicity 6.
The water oxygen is 2.6 Å from the target in the inner pose and 4.5 Å in the
outer pose. The two orientations share the predeclared, geometry-selected source
ray. All energies use the unchanged native MACE OMOL vacuum plus standalone
GFN2 ALPB(Water) minus GFN2 vacuum recipe; GFN2 retains effective UHF 0.

`H = (E_Dy,inner − E_Dy,outer) − (E_La,inner − E_La,outer)`.

| Source | Orientation | La inner−outer | Dy inner−outer | H |
|---|---:|---:|---:|---:|
| 7CCO | 0° | +5.023217 | +2.480628 | −2.542589 |
| 7CCO | 90° | −4.151979 | −6.627693 | −2.475714 |
| 7CCN | 0° | unavailable | unavailable | unavailable |
| 7CCN | 90° | unavailable | unavailable | unavailable |

Values are kcal/mol, copied from the completed Slurm summary. At 0°, inward
placement is uphill for both metals, with a smaller penalty for Dy. At 90°, it
is downhill for both, with a greater gain for Dy. The relative effect has the
same sign across these orientations despite their different absolute transfer
signs. MACE contributes −2.108176/−2.099925 kcal/mol to H; the solvent correction
contributes −0.434413/−0.375789 kcal/mol for 0°/90°, respectively.

Both metals select `inner_plane90` from the complete shared four-pose 7CCO menu.
The 7CCN menu remains unavailable: its four outer cells scored successfully, but
the four inner cells failed the original contact gate and were never evaluated.
No alternate placement replaced them. Geometry rejection is evidence about this
frozen construct and placement rule, not proof that the molecule cannot bind water.

## Meaning and stopping decision

The admitted 7CCO poses demonstrate a conditional energetic response to explicit
water placement that favors Dy relative to La within this model. The peptide and
Asn5 did not move, so this does not establish water–Asn displacement, equilibrium
occupancy or a coupled structural mechanism. The modeled water is not deposited
evidence. The outer pose is not a bulk chemical-potential reference; no raw energy
minimum was taken across dry and wet inventories. Missing 7CCN transfers prevent
a two-source conclusion. No affinity, population or natural-protein classifier
follows from these cells.

**Close the declared panel.** No new placement search, optimization, altered
contact threshold or automatic water campaign follows from this result.

## Execution and preserved evidence

| Stage | Job | Terminal state | Scheduler seconds | Allocation | New MACE / GFN2 calls |
|---|---:|---|---:|---|---:|
| Geometry | 12565 | COMPLETED, 0:0 | 2 | 1 CPU | 0 / 0 |
| Paired scout | 12664 | COMPLETED, 0:0 | 34 | 12 CPUs, 1 H200 | 2 / 4 |
| Remaining admitted cells | 12670 | COMPLETED, 0:0 | 119 | 12 CPUs, 1 H200 | 10 / 20 |

Actual panel total: **12 MACE and 24 GFN2 evaluations**, no gradient calls or
optimization. The remaining stage reuses the two scout cells. No old 256-atom
energy was reused. Application wall times were 33.158 and 117.907 seconds for
scout and remainder. Recorded GFN2 CPU times were 33.906 and 160.818 seconds;
MACE forwards took 1.002 and 1.732 seconds, with one model load per allocation.
These timings distinguish measured component work from reserved resources.

- [Compact result and all 16 cell statuses](RESULT.json), SHA256 `cd8f8deea79c119e0e5ad04d3091f47fdaccc6c47f65069151955788aca6e6f1`.
- [Terminal scheduler accounting](ACCOUNTING.txt).
- [Geometry plan](../LBT3_WATER_PLAN.md) and [score plan](../LBT3_WATER_SCORE_PLAN.md).
- [Actual final result](../../../workspaces/lanm_coupled_response_20261016/lbt3_water_remaining_v1/execution/RESULT.json), SHA256 `aa4dc92c0a62d7579f84c4c2f0f6e559f3ab275f8c13417f1b1f07eaea6fa63a`.
- [Actual scout result](../../../workspaces/lanm_coupled_response_20261016/lbt3_water_scout_v2/execution/RESULT.json), SHA256 `26aa78c269cd954b63f48ad034771058d72ccd1aae35fd60376a2715466ee18f`.

This report extracts already saved results and receipts; it performs no new
energy arithmetic, molecular evaluation or placement selection. Production and
the original preparation exclusions remain unchanged.
