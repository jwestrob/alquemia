# Joint metal/donor response: real force relief, no demonstrated classifier gain

Allowing the metal to move removes the diagnosed native-MACE force in the seven
successful searches. It does not improve the available class decisions. Keep
this as research support, not a routine correction or a promoted classifier.

## Actual outcome

The eight original starts used the exact prepared seven-coordinate rule.
Seven produce valid interior candidates. Metal displacements are0.137–0.380Å;
final maximum heavy-atom displacements are0.257–0.743Å. Metal-gradient norms fall
to0.000145–0.000869kcal/mol/Å. These are stationary responses of the **vacuum-MACE
proposal energy** within the selected coordinates, not solvent-composite minima.
Omitted angular loads remain, especially PLM07ab La's Glu chi1.

PLM8344 La fails at an oversized trial: the metal/context-atom37 separation moves
from3.013Å to0.937Å, violating the unchanged1Å severe-overlap guard. Its maximum
trial heavy displacement is1.961Å. No MACE energy is evaluated at this rejected
geometry. The original input is valid; this is an optimizer-trial failure, not
evidence that the protein is chemically invalid. No fallback or rescue is used.

As specified, intermediate SLSQP requests may exceed the final displacement
limit. There are32 completed infeasible requests; the largest is3.804Å. All seven
admitted final candidates pass the actual source/cap/bond/overlap checks and
leave nonselected atoms fixed. The failure remains in the four-source denominator.

Each successful pair competes over the same seven real geometries, retaining
all five prior candidates and their energies. Both endpoints select their own
joint proposal. All six new cross-MACE and24GFN2 cells succeed.

| Source | ΔR versus angular pool, model kcal/mol | Old-band transfer |
|---|---:|---|
| 1H4I | +0.762495060 | Ca-supported; unchanged |
| 4MAE | −1.646585986 | La-supported; unchanged |
| PLM8344 | unavailable | No joint decision; prior score separate |
| PLM07ab | +0.751242586 | Inconclusive; unchanged, truth unknown |

There is no joint calibration. In particular, PLM07ab's proximity to an old
boundary is not a reason to alter that boundary or call its class known. The two
crystal contrasts move closer together. Together with the wider angular result,
this does not justify another full-fold campaign simply to enlarge the search.
The failed boundary case remains a concrete limitation of this implementation.

## Cost and reproducibility

Jobs1209901/1209902/1209903 are complete. Proposal execution made152 actual MACE
calls, with15 cache reuses and one geometry-rejected request. Six additional
cross evaluations and24 native GFN2 calls completed the pool. Total incremental
allocation: **6,880 core-seconds and77 GPU-seconds**, including the failed search;
no new DFT. Earlier origins, angular searches and their pools are reused, so this
is not end-to-end source latency. Historical production timing is separate.

Protocol: `nikasha_joint_metal_angular_common_geometry_native_OMOL_GFN2_ALPB_v1`.
Manifest: `workspaces/adaptive_metal_20260922/common_pool_v1/manifest.json`, SHA256
`6f5c747813ef21b2262f38b5f3e79ab36c98d58bbee986998d93eb7761f830df`.
Actual result: `after_solvent_0_1209903.json`; exact components, choices, forces,
old-reference decisions and cost receipts: `summary.json` in the same directory.
The generator is [summarize.py](summarize.py). [Execution scope](EXECUTION_PLAN.md)
retains the prior preparation-only manifest and its separate launch agreement.

Four real joint-pool tests pass with zero skips. They check retained prior
energies, identical geometry pools, exact reused endpoint states, unavailable
failed searches and actual row-minimum algebra. The preceding22 pool/comparison
regression tests also pass. These are implementation checks; the molecular
evaluations above provide the physical evidence.

No production default, old result, reference or PLM cohort was changed.
