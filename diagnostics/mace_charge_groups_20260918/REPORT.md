# Whole-protein charge-group POLAR: PQQ improves; broader screen fails

The declared `intact_POLAR_medium_chemical_group_charge_frozen_OBC2_v1`
model passes **1/7 directions**, all 30 numerical checks and all three grouping
checks. The PQQ pair improves from -17.628168 to +19.957810 model kcal; all six
alpha/GGR directions remain wrong. These are four consumed biological groups,
not seven independent or blind observations. Production baseline/default unchanged.

## Matched comparisons

Larger R = E(Ca)-E(La) is more La-like. Entries are differences between sites,
in model kcal; the common metal-energy offset cancels. No absolute class, aquo
reference, universal zero, or free-energy interpretation is supplied.

| Expected higher / lower | Previous POLAR + GB | Grouped vacuum component | GB contribution | Declared grouped total |
|---|---:|---:|---:|---:|
| PQQ_4MAE / PQQ_1H4I | -17.628168 | 51.626321 | -31.668511 | 19.957810 |
| ALPHA_1F6S / GGR_1GLG | -21.259552 | 104.518476 | -129.853476 | -25.335000 |
| ALPHA_1F6S / GGR_2FW0 | unavailable | 45.818337 | -91.695198 | -45.876861 |
| ALPHA_1F6S / GGR_2FVY | unavailable | 54.973349 | -89.268989 | -34.295640 |
| ALPHA_6IP9 / GGR_1GLG | -35.012731 | 92.758266 | -134.096976 | -41.338711 |
| ALPHA_6IP9 / GGR_2FW0 | unavailable | 34.058127 | -95.938699 | -61.880572 |
| ALPHA_6IP9 / GGR_2FVY | unavailable | 43.213139 | -93.512490 | -50.299351 |

The unavailable native scores were not recomputed or filled with another model.
PQQ labels are functional class evidence. Alpha/GGR compares a qualified
cross-study alpha strong-site direction with GGR's direct same-assay direction;
it is not a matched experimental free-energy difference. Alpha 6IP9 is already
La-bound. Keep both alpha and all three GGR structures in the denominator.

## What changed and what was learned

The intact protein, source coordinates, protonation, cofactor/water inventory,
learned weights and all energy readouts were retained. The three POLAR charge
restoration stages now conserve charge within source-defined chemical groups.
The same changed endpoint density feeds the full model and frozen-charge OBC-II.
No carve energy, new DFT, geometry search or fitted label correction is added.

This demonstrates a PQQ ordering gain over the previous whole-protein POLAR
model. It does not establish an improvement over the production DFT baseline,
which already separates those PQQ examples. The tested solvent model still
fails the broader objective and should not be expanded as a validated scorer.

The vacuum component gives all seven expected directions, while GB reverses
all six alpha/GGR directions. This locates the reversal algebraically; it does
not prove the solvent physics is wrong or justify omitting it from an affinity
calculation. A separately frozen vacuum compatibility descriptor is a reasonable
next accuracy test, with these inspected cases designated method development.
It would not retroactively convert this declared total from 1/7 to 7/7.

## Physical and numerical checks

The one-group limit reproduces both native energies exactly at saved precision.
Maximum native density difference is 4.802e-15 e; native force difference
1.510e-13 eV/A. Grouped force rotation error is at most 1.420e-7 eV/A.
The 30 declared energy/rotation/group-charge checks pass; maximum reported
energy error is 1.9974e-6 model kcal. Total charge error is at most 4.139e-13 e;
minimum group Fukui denominator ratio is 0.62321, far from the frozen 1e-10 gate.

| GGR structure | Vacuum connected-minus-primary | Total connected-minus-primary | Total within +/-2 |
|---|---:|---:|---|
| GGR_1GLG | -2.513855 | -0.077926 | yes |
| GGR_2FVY | -3.275712 | -0.083309 | yes |
| GGR_2FW0 | -3.116626 | 0.021714 | yes |

All vacuum grouping changes exceed the existing 2-model-kcal criterion.
Preserve this failure if pursuing the vacuum component separately. The solvent
term cancels most of that grouping sensitivity; its successful representation
check did not imply predictive accuracy. Formal group charges are an unvalidated
constraint on learned density coefficients. Neither an exact charge partition
nor self-consistent solvation is claimed. Combined solvent gradients and
mechanical corrections remain unavailable.

## Execution and costs

24 MACE and 24 OBC-II evaluations completed (1201508 and 1201513). Preserved
setup failures: 1201505 failed dispatch before a forward; 1201507 attempted one
forward and ran out of GPU memory after bypassing the old unused-grid wrapper.
V3 restores wrapper composition; all scientific inputs/settings match V2.

Including both failures: **895 allocated GPU seconds (14 min 55 s), 14,320
allocated core-seconds, 987.336 reported CPU-seconds**. Successful MACE inference
sums 660.037016 s; solvent inference 1.400081 s. Per-pair model time: alpha
17.38–18.10 s, GGR 49.03–49.39 s, PQQ 113.35–117.72 s. Peak successful MACE GPU
allocation 10,207,527,424 bytes; failed attempt 16,392,278,528 bytes. OpenMM GPU
memory was not measured. Source preparation 11.245473 wall / 11.105136 CPU s;
local tests/preflights/reporting and historical reuse costs are additional.
This is development cost, not an established production-overhead comparison.

Six real-fixture tests pass, zero skips, 45.957 s: physical group inventory,
corrupted charge/weight rejection, native energy/density/force identity,
actual full report replay, sign/conversion algebra, explicit missing-solvent
behavior, saved charge traces and preservation of the execution wrappers.
These parser/replay checks reuse actual scientific outputs; they do not create
new molecular evaluations. Earlier twelve runner regression tests also passed.

## Judgment and next operation

- Numerical consistency: passes the declared checks; physical approximation
  limits remain.
- Predictive usefulness: PQQ gain relative to native whole-protein POLAR;
  no broad gain in the declared seven-direction screen.
- Affordability: practical A5000 execution; accuracy remains the limiting issue.

Retain the production baseline. Preserve the group adapter as a useful research
building block. Next investigate a separately declared vacuum descriptor on
additional biological groups, retaining the observed grouping failure, and
compare environment terms without score-dependent parameter changes. No new
pilot is launched by this report. Exact replay commands are in COMMANDS.md;
full endpoints, receipts, snapshots and figures remain in the workspace.
