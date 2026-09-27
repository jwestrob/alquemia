# Classical components: derivatives pass; whole-protein relaxation remains unsupported

All 32 real configurations completed and all 200 frozen implementation checks pass.
This supplies the missing classical mechanics component of the finite additive model.
It does not qualify the electronic component, solve the absent-solvent problem, or
demonstrate improved La/Ca classification.

## What changed physically

The 1H4I Thr159 A-to-B movement changes retained bonded energy by -0.1394726824
and MM electrostatics by +0.5574093753 kcal/mol. Both changes are identical for Ca
and La; the classical contribution to this particular response double difference
is exactly zero in the recorded arithmetic. Thus the classical terms do not erase
the previously observed electronic response. This is specific to the fixed core
and this hydrogen-only movement, not a general cancellation of protein mechanics.

At source A, real-QM/MM Lennard–Jones energies are -57.1460190077 (Ca) and
-57.5284098364 (La) kcal/mol, a Ca-minus-La contribution of +0.3823908287.
These are explicitly unqualified pure-12-6 cross parameters, not affinity evidence.

## Numerical evidence

| Test | Largest absolute residual | Units |
|---|---:|---|
| Repeat energy | 0 | kcal/mol |
| Rigid energy | 7.40328687243e-10 | kcal/mol |
| Rigid force | 5.06406649947e-09 | kcal/mol/Å |
| Directional finite difference | 7.84644959708e-05 | kcal/mol/Å |
| Step refinement | 5.88552211411e-05 | kcal/mol/Å |

All tolerances remain those in COMPONENT_CHECK_PLAN.md. Native electronic rigid
rotation failures remain failures. A combined-force comparison is being assembled
from exact compatible saved geometries; no additional quantum calls are implied.

## Preparation warning

The large whole-protein bonded/LJ totals require interpretation before movement.
Read-only diagnosis finds a severe contact in the rebuilt, experimentally missing
A596 tail. A source-strain report is in progress. Do not relax this dry model or
reinterpret fixed-source contrasts as equilibrium preferences.

## Execution and next work

Scout1219495: 11 s ×1 CPU. Remaining31 evaluations1219497: 21 s ×32 CPUs;
the exact Ca_A cache was reused. Collectors1219496/1219498: 1 and2 s ×1 CPU.
Total allocation for this classical stage including collectors: 686 CPU-seconds,
zero GPU-seconds. Maximum recorded per-process lifetime peak RSS: 305812 KiB;
this is not summed node peak memory. Both allocations used --mem=0.

Results: `workspaces/metal_environment_response_20260926/component_checks_v1/FINAL_COLLECTION_1219497.json`.
The actual completion event was received and acknowledged. All four jobs are terminal.
Next: combine compatible saved electronic/classical force checks, finish source-strain
diagnosis, and establish honest La/Dy electronic-state capability before any LanM scout.
No production change, full-protein optimization, new quantum job or library rescore.
