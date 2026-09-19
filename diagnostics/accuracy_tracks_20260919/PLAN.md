# Parallel accuracy development — 2026-09-19

## User agreement and objective

Jacob rejected training a surrogate of the existing DFT method as the main
accuracy-development direction: “The point is to improve, not make one method
a faster version of the other.” The target is **metal specificity**, including
noncatalytic binders such as LanM, not catalytic competence.

After discussing chemical-state competition, protein-environment effects and
errors in electronic treatment, Jacob explicitly requested: “Explore all those
in parallel with subagents. Remember we can do lots of compute and the whole
point is to help us climb this hill.” This authorizes contained investigations
and actual calculations on all three tracks. Existing baseline, scientific
records, other agents' work and scheduler rules remain protected. Finite task
manifests make experiments inspectable; there is no arbitrary project CPU/time
budget. Preserve failed attempts and report measured execution costs.

## Ownership and concrete work

| Track | Agent | Initial work |
|---|---|---|
| Competing chemical states | water_basins | Finish reusable contextual water preparation, then test physically justified proton/water-state competition with balanced reference accounting. Reuse the occupancy calculations; the failed wider-basin approximation cannot supply missing entropy. |
| Protein environment and constraints | second_shell | Execute the frozen ten-endpoint context-supported coordination pilot, then use its result and earlier environmental failures to select a coherent next physical model if warranted. |
| Electronic-treatment errors | khoury_benchmark | Choose and execute a compact independent electronic-method test on matched real cores; distinguish method sensitivity from evidence of improved accuracy. |

The parent reviews the interpretation and coordinates follow-ons. Agents own
separate scripts, tests, diagnostics and workspaces. They send short scientific
scope/scale notes before new executions for coordination, without reopening
already granted user authorization. Each track writes a vault note.

## Common scientific comparison

- Begin with consumed real references and their archived preparations. Preserve
  experimental evidence strata and group structural replicates by protein.
- Separate chemistry, geometry, environment and electronic-method changes in
  matched comparisons. State changes must have their real energetic accounting.
- Record the physical target and selection rule before seeing pilot outputs.
  A failed result is not permission to select favorable states or methods.
- Preserve PQQ discrimination while seeking improvement on known difficult
  comparisons. Lower endpoint energies, a larger uncalibrated model scale, or
  agreement with current DFT alone do not demonstrate improved classification.
- Changes do not inherit baseline bands. Report relative site contrasts when
  no compatible absolute reference/calibration exists.
- Keep missing terms unavailable; do not fill occupancy, entropy or environment
  failures with zero or silently substitute successful baseline values.
- No production/default change, broad PLM rescore, push or deployment.

## Immediate launched scope

The parent reviewed and approved
`diagnostics/coordination_preparation_20260919/PLAN.md`: five real structures,
Ca/La pairs, ten native MACE constrained context optimizations, ten original-core
MACE evaluations and up to ten original-core native CPCM-DFT single points.
The agent is to execute and report actual discrimination, not just preparation.
Other tracks will record their exact methods/manifests in their own diagnostic
directories before execution; this file does not invent their results or inputs.
