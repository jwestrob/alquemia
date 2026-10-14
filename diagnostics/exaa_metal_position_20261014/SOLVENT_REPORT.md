# Ion-placement response survives solvent; score benefit remains unproved

All eight new solvent cells pass the existing component-v2 and native replay
admission. The reciprocal donor-response effect survives the frozen-charge
solvent descriptor. However, the secondary limited pool gains only 1.141070
kcal/mol for the target and zero for the control. This does not yet show a gain
over the actual historical three-candidate scorer: that comparison must also
include its saved adaptive_Ca geometry.

## Actual response, kcal/mol

W_M is saved adaptive_La minus origin at one ion placement. D=W_Ca−W_La.
The metal point is the same for both donor arrangements at each placement.

| Host | Quantity | Original ion | Transferred ion |
|---|---|---:|---:|
| Control0 | Ca total work | −6.395494 | −11.761451 |
| Control0 | La total work | −12.934776 | −23.023109 |
| Control0 | D total | +6.539282 | +11.261659 |
| Target3 | Ca total work | −7.898852 | +1.426089 |
| Target3 | La total work | −26.390560 | −7.171611 |
| Target3 | D total | +18.491708 | +8.597700 |
| Target minus control | D total | +11.952426 | −2.663959 |

Control delta_D is +4.722377: +3.373337 MACE plus +1.349040 solvent. Target
delta_D is −9.894008: −9.328603 MACE plus −0.565406 solvent. Both meet the frozen
reciprocal interpretation rule. All native electrostatic, hydrogen-bond, SASA
and reference components are retained in the unrounded result; none was fitted.

## Why the large response change is not a large score improvement

The primary diagnostic changes the energy difference along the same donor path.
Its origin energy changes too. Giving each metal access to every compatible
geometry is the relevant separate test for scoring.

In the secondary four-geometry pool, both control metals retain their original
adaptive_La geometry. For the target, La also retains its original adaptive_La
geometry; only Ca selects the transferred ion with the original donor arrangement.
The target's paired contrast shifts −1.141070 kcal/mol against the two-candidate
baseline. That is much smaller than the change in donor work. It is not a class
transition or evidence of better predictions.

The actual prior scorer also contains adaptive_Ca, omitted by design from the
donor-response diagnostic. Its exact cached results are available. Root is
checking the full historical three-candidate pool versus all five compatible
geometries, plus saved physical forces, before deciding whether a new search is
worthwhile. These follow-ups use saved outputs only; no new molecular evaluations.
Do not call the two-to-four comparison a gain over the released/current scorer.

## Model, cost and limits

E = MACE_OMOL vacuum energy + native ALPB evaluated at the corresponding
standalone GFN2 vacuum charges. All 192 atoms, proton/cofactor state, water
inventory, caps and metal charges remain as declared. These are modeled
perturbations of consumed ExaA folds, not equilibrium structures or affinities.
The native factory parameters and published standalone parameter file are pinned;
the corrected ORCA full JSON is provenance only, not an executed parameter input.
Historical absolute xTB/tblite parity failure remains separate and failed.

Job4132 completed32 SCC calls and16 native component replays, with zero new MACE,
DFT or optimizer calls. Eight original components and all16 MACE values were
reused. Scheduler wall166s,7CPUs:1162 allocatedCPU-s; TotalCPU928.873s,
batchMaxRSS714460K. The enclosing Python receipt measures166.637s/1166.461
allocatedCPU-s; integer scheduler timing explains the scope difference.

Full result: `workspaces/exaa_metal_position_20261014/solvent_v1/execution/FINAL_COLLECTION.json`.
Exact config, sources, execution and per-cell receipts are pinned by that record
and indexed in SOLVENT_RESULT.json. No bands, calibration or production default
changed. All molecular jobs in this test are terminal; delayed callbacks must
not trigger reruns.
