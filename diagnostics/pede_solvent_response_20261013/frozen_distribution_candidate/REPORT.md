# Fixed-charge solvent candidate: real PedE gain, incomplete transfer

14 October 2026 UTC. **Both wrong two-Ca PedE calls become correct, and all five
two-Ca PedE sources now classify correctly.** However, five previously correct
three-source summaries become inconclusive. Pursue this candidate as an explicit
alternative; it has not earned blanket replacement of the current scorer.
Production, old calibration and every original failure remain unchanged.

## What changed and why

The completed DFT comparison supported a real source-dependent solvent response
but found a substantially smaller contrast than self-consistent GFN2/ALPB.
Native component replay localized most GFN2 sensitivity to additional solute
charge redistribution, including its intrinsic energy compensation. Neither
result established that all polarization should be removed.

This finite experiment tested a different, declared approximation:

`E(M,q) = E_OMOL,vac(M,q) + G_ALPB(q_GFN2,vac(M,q); q)`.

The maintained xTB 6.7.1 native ALPB routines supply electrostatics, hydrogen
bonding, SASA and reference terms. The reaction field responds to fixed solute
charges; the solute density does not re-equilibrate. This is a perturbative
descriptor, not a self-consistent solvent treatment or binding free energy.
MACE weights and all source states, waters, protons and geometries are unchanged.
Every source retains the same complete common candidate pool for both metals.
Both rows are selected again under the new expression, with the unchanged
0.1 kcal/mol operational rule. No response multiplier or target fitting was used.

All charge/MACE inputs were reused from the pinned corrected-parameter lineage.
This wrapper invokes only the native solvent functional, not an ORCA or xTB SCC
calculation. It does not alter the corrected La orbital parameter assignments.

## Numerical implementation and calibration

All 1,542 distinct geometry/metal inputs are available. The 3,076 new component
evaluations and four reused input pairs pass the unchanged per-component
water-charge reproduction checks and charge-rounding limits. Maximum water
component residual is 3.73413e−8 Eh. Cavity and Born matrices are identical
between the two charge distributions for each input. This validates the saved-
charge replay, not physical adequacy of freezing the density.

The exact 25 designated canonical members were evaluated first. All remained
separated under the existing calibration rule; operational gap 8.4503234703
kcal/mol. New bands were frozen in `CANONICAL_COLLECTION.json`, SHA256
`ee00fa03aba64c23e6a737fb2b3b4a9f40dc139e58bdb9225e8ed10c9af2b58d`,
before any transfer evaluation. The old gap was 3.4721710470 kcal/mol; increasing
this training gap alone is not predictive evidence. Raw R remains Ca minus La,
larger more La-like; its large model-specific offset is not an affinity scale.

The scorer independently reproduced all archived selections/calibrations over
278 records, all 25 canonical inputs and all 100 grouped comparisons. Missing-
input checks used explicitly corrupted copies of those real fixtures. Static
review was independent of both implementation authors. Actual molecular-component
evaluations, parser/algebra tests and static review are distinct evidence.

## Fixed-band transfer

Counts are **correct / wrong / inconclusive / unavailable**. Populations overlap;
they must not be summed as independent proteins.

| Population | Corrected baseline | Fixed-charge candidate |
|---|---:|---:|
| Canonical 25, calibration | 25 / 0 / 0 / 0 | 25 / 0 / 0 / 0 |
| Expanded 135 source panel | 134 / 1 / 0 / 0 | 129 / 0 / 6 / 0 |
| Original 125 Ca-source inputs | 105 / 1 / 1 / 18 | 105 / 0 / 2 / 18 |
| Same 125 with documented recoveries | 119 / 3 / 3 / 0 | 121 / 0 / 4 / 0 |
| All 15 two-Ca sources | 12 / 2 / 1 / 0 | 13 / 0 / 2 / 0 |
| Existing 100 three-source summaries | 94 / 0 / 0 / 6 | 89 / 0 / 5 / 6 |

The exact original source denominator and recovery identities are retained.
There is no new voting policy for the Ca-source set. All new replay cells are
available; the unavailable original inputs and six excluded triples remain
pre-existing preparation exclusions. These are consumed development structures,
not fresh blind tests or independent biological replicates.

### What improved

- Two-Ca PedE samples 1 and 2 change from wrong La calls to correct Ca calls.
  All five are now correct. Their accommodated contrast range drops from
  19.9420 to 5.78134 kcal/mol. Static origins also give five correct calls under
  the candidate's separately calibrated origin bands; new proposals were not
  needed for that result.
- Repaired MMOL_2048 changes from inconclusive to correct La. Its five-source
  range drops from 13.5653 to 5.92343 kcal/mol.
- A0ACD6B9F2 Ca-conditioned sample 4 changes from wrong Ca to inconclusive.
  This is an abstention, not a recovered correct prediction. No formerly correct
  biological La-class call regresses in this comparison.

### What did not improve

- La-conditioned PedE samples 3/4 and ExaA Q9Z4J7 samples 2/3/4 change from
  correct Ca to inconclusive. These are Ca-class proteins folded with La;
  source conditioning and biological reference class are distinct.
- Two-Ca ExaA sample 3 also changes from correct to inconclusive. Its source
  range slightly increases, 8.56476 to 8.92447 kcal/mol. Structural robustness
  therefore does not improve uniformly.
- Existing PedE triple_078 and ExaA triple_096–099 become inconclusive despite
  complete input coverage. The previous all-available-triples fidelity is lost.
- MMOL_1770 Ca sample 0 and ExaF C5AXV8 two-Ca sample 4 remain inconclusive.

The zero-wrong count combines genuine PedE corrections with abstention. Old-band
transfer is retained per source and policy in the result/table; it is not silently
reported as success under the old scale. Operational and mathematical decisions
agree here. Origin-only results remain separate: among the 100 historical triples,
they give 93 correct, one inconclusive and six unavailable. Selecting whichever
policy works for each protein is not part of this model.

## Cost and execution

| Job | Work | Wall s | Allocated CPUs | Allocated CPU-s |
|---|---|---:|---:|---:|
| 4006 | Real archived regression | 1 | 1 | 1 |
| 4007 | Thin compilation + canonical components | 8 | 12 | 96 |
| 4008 | Canonical calibration | 2 | 1 | 2 |
| 4009 | Transfer components, reused build | 51 | 12 | 612 |
| 4010 | Fixed-band comparison | 8 | 1 | 8 |
| 4014 | Saved-result table/group/component summary | 11 | 1 | 11 |

Total 730 allocated CPU-seconds; zero GPU, new MACE, SCC or DFT evaluations.
Native child processes report 173.9 measured CPU-seconds and at most 10,752 KiB
RSS per child; this is not whole-job peak memory or measured allocation utilization.
These timings reuse expensive archived charges/proposals. They do not establish
new-source end-to-end latency or a production speedup.

The prior DFT/replay phase cost remains separate and fully recorded. All new jobs
are terminal and acknowledged; all completion wakes were armed. Older unrelated
reporting jobs remain preserved under Jacob's no-cancellation instruction.

## Decision and next question

**Pursue the fixed-charge candidate, without blanket promotion.** It gives a
concrete, affordable PedE correction and improved source stability, while exposing
an unresolved transfer tradeoff. This supports investigating the solvent-response
model rather than attributing the PedE failures solely to bad folding.

The [saved-matrix diagnosis](selection_diagnosis/REPORT.md) is also complete as
job4016: the new abstentions combine changed solvent contributions and a stricter
canonical Ca boundary. Reselection does not cause them; where it changes a
regressing source, it helps. PedE sample1 is corrected with identical selections;
sample2 additionally benefits from a different saved Ca candidate and has only
a 0.23276 kcal/mol correct-class margin. No band was retuned.

That diagnostic adds one allocatedCPU-second (731 for this candidate round).
The next physical question is how much electronic response is reliable across
different pockets. Do not automatically repeat DFT, refold these proteins,
switch models by desired label or interpolate response to force correct calls.

Results/pins: [EVIDENCE.json](EVIDENCE.json). Executable records and restart
instructions: [COMMANDS.md](COMMANDS.md). Full machine-readable per-source contrasts,
decisions, selected geometries, works and missing reasons remain in the workspace
`FINAL_COLLECTION.json` and `summary_v1/SOURCES.tsv`.
