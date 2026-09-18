# Simultaneous substitution: numerically sound, no accuracy gain

All8 actual forwards completed and all23 numerical/algebra checks pass.
Changing every modeled metal simultaneously leaves the author-domain challenge
at **6/9 comparisons** and makes RTX less La-like by10.038875model kcal per
metal. This static correction does not repair the RTX failure.

| Protein | Sites | Single-substitution mean | Simultaneous per metal | Change |
|---|---:|---:|---:|---:|
| A0A7 | 6 | 49.234436 | 48.019569 | -1.214867 |
| HEW5 | 8 | 70.120351 | 69.505120 | -0.615231 |
| RTX | 8 | 2.777993 | -7.260882 | -10.038875 |
| PARV_4CPV | 2 | 45.889415 | 45.376832 | -0.512583 |
| AEQ_1SL8 | 3 | 32.595055 | 32.124324 | -0.470731 |

Units are atom-referenced descriptor model kcal, not binding free energies.
All ordered original site scores are retained in the full result. A0A7/HEW5
still exceed all three GGR structures; RTX fails all three. Parvalbumin's joint
mean passes2/3 supporting GGR comparisons, as did the old mean; this is distinct
from the original4/6 site-wise comparison. Aequorin remains an ordered vector
plus this joint descriptor, without individual labels or an absolute decision.

## What this establishes

The fixed-coordinate non-additivity is large for RTX and opposes the desired
ordering. It is smaller in the other four proteins. This rules out this tested
static all-site substitution as a remedy. It does not show that experimental
cooperativity is absent: the calculation retains Ca-conditioned geometry and
full occupancy and contains neither folding nor a binding/titration ensemble.

The three author domains have La ITC affinities and Ca CD folding thresholds
under different conditions; the latter are not Ca Kd. The benchmark remains
qualified protein-level cross-readout evidence. GGR structural replicas form
one biological group; the nine comparisons are not nine independent tests.
Parvalbumin is supporting cross-study evidence; aequorin is weak protein-level
directional evidence with unresolved sites. Source MPVP-scar, pH, modeled-site
count and cofactor/water exclusions remain as recorded in the parent studies.

## Implementation and checks

New opt-in protocol `masked_omol_whole_protein_collective_metal_substitution_v1`.
Same actual complete proteins, coordinates, protonation and native masked OMOL
checkpoint/evaluator. The old one-selected-metal guard remains unchanged.
Collective tasks explicitly map every metal and increase total charge by the
site count, preserving even electron parity. No force-field La charges or
global neutralization was invented. No baseline or historical score changed.

The [plan](PLAN.md) records the subtraction algebra and acceptance rules.
All-Ca and single-La states are reused only from matching actual receipts.
Five all-La endpoints plus RTX Ca replay, La rotation and La permutation ran.
Readout closure, rigid rotation/permutation, all-Ca replay, one-site GGR identity
and non-additivity algebra pass. This is a model descriptor, not quantum/PB
energy, an aquo-reference score or a reproduction of assay cooperativity.

## Cost and reproducibility

Job1201453: **115GPU-allocation seconds,1840allocated core-seconds,
140.033reported CPU-seconds** on oneA5000/16CPU/64474MiB. Eight model forwards
took43.179467s; peak GPU allocation3939399680bytes. Zero newDFT/solver/trajectory
or fitting calls; zero scientific execution failures. Successful preparation
8.861036wall/6.571481CPU-s. One initial dictionary-layout preparation error
occurred before output creation/inference and is recorded; technical recovery
preserved the same inputs. Other local work and historical reused calls remain
additional costs, not zero.

Three real-fixture tests pass23.181s, zero skips, including actual report replay
and rejection of corrupted metal inventory. [Compact results](RESULT.json),
[commands](COMMANDS.md), [PDF](../../workspaces/mace_collective_20260918/figure_v2/comparison.pdf),
[SVG](../../workspaces/mace_collective_20260918/figure_v2/comparison.svg).
Full results, source pins and execution receipts live in the matching workspace.

**Recommendation:** retain the baseline and close this collective substitution
as an RTX accuracy fix. It is affordable and numerically coherent, but it did
not add predictive value on the declared comparisons. The broad MACE goal
remains active; neither this result nor the failed coupled response completes it.
