# Joint donor response with the passing water model frozen

## Question and actual trigger

Can the previously declared eight donor pair parameters reproduce both metals'
wet and dry ligand responses on the consistent PBE0-D4 development surface?
Water fit17605 passes all existing consumed water controls. Frozen transfer17612
passes its ten Cartesian derivative checks but fails51 of57 hydrated-donor
controls, including all27 rigid-body force controls. The largest hydrated balanced
exchange error is14.162637 kcal/mol. All nine dry rigid-body controls fail exactly
as before. Improving water alone therefore does not remove the represented donor
errors. This does not prove that changing these eight parameters will be adequate.

This declaration activates `paired_donor_development_v1/DESIGN.md` and its additive
`REFERENCE_ELIGIBILITY_AND_EXECUTION_OUTLINE.md`, with a concrete budget and the
actual passing water candidate. Their historical water-candidate references are
not used as this run's baseline. Root already accepted the consistent28 PBE0-D4
rows as provisional development targets. The limited central-amide derivative
qualification and all electronic-state/CC/counterpoise limitations remain.

## Frozen inputs and parameter family

Use the unchanged28-row training result99ea27f8, body helperbf65a20a and symmetric
objective86e85fcb from `paired_donor_development_v1/objective_v1/OBJECTIVE.py`.
Use actual17612 results as the starting model:18 hydrated donors, six matched-dry
donors and four water/aquo centers. Preserve the other four radial-water records
from17605 as invariant controls outside the donor objective. Admit exact native
provenance and input identities before any new calculation.

Freeze all six selected17605 water parameters:
`[3.906136716492464,0.28782327521600637,3.5865909437359376,
0.2837564269122712,0.20765190997435956,0.3173248056679709]`.

Eight variables are radius/damping for each of four pairs, in the order below.
Mono- and bidentate acetate share a single carboxylate pair for each metal.

| Pair | Initial radius, Å | Initial damping | Fixed epsilon, kcal/mol |
|---|---:|---:|---:|
| La–carboxylate O |3.6318|0.1518|0.2280|
| La–amide O |3.3000|0.2987|0.3019|
| Dy–carboxylate O |3.347237658604727|0.139|0.3795688850865599|
| Dy–amide O |3.2657502966294456|0.139|0.3795688850865599|

Keep original resolved-radius bounds0.85–1.15 times these initial radii and
damping0.05–0.50. These are development limits, not uncertainty intervals.
Preserve the existing small-system class/type mapping, fixed donor epsilons,
La amide-C override, Dy amide-C mixing, water-H, multipoles, polarizabilities,
bonded terms, coordinates, chemical states and native numerical settings.

## Objective, algorithm and finite execution

Do not modify the prepared objective: four equal energy/La/Dy/paired blocks;
energy subweights0.50 aquo-referenced exchange,0.30 matched wet-minus-dry,
0.20 individual radial work. Keep both metals and every declared fixture.
Energy controls remain≤1 kcal/mol; per-fixture rigid RMS≤2 and maximum component≤5
kcal/mol/Å. Paired radial-work controls remain required although not added as a
separate loss term. Water energies/forces must remain exactly unchanged.

Use the existing dimension-general coupled direction solver: central probes±0.005
of each original normalized interval, with actual clipping/widths recorded;
relative least-squares cutoff1e-6; normalized component trust box±0.025; four
half/full scaled Gauss–Newton and bounded-least-squares trials in existing order.
Complete each batch before selection. Accept only actual objective improvement;
the current point wins ties. Derivative probes cannot become selected candidates.

Run at most24 steps,480 new vectors (481 including the initial vector), and9216
native programs. Sixteen derivative vectors cost at most192 programs: eight
carboxylate probes affecting eight cells each and eight amide probes affecting
four cells each, analyze plus analytic testgrad. Four full trials affect24 donor
cells each, at most192 more programs. Thus384×24=9216 before exact reuse.
Initial28 values and all water controls must replay with zero new calls.

Use one eight-CPU standard_expand task, up to eight single-thread workers and
default memory/time. Cache identity includes the effective donor and frozen water
parameters, exact inputs and native provenance. Reuse unchanged-metal and
unchanged-donor-family cells; preserve all old results/ledgers. Stop new starts on
native/admission failure, drain siblings and retain the full unavailable denominator.

Stop on all original controls passing, a complete trial batch with no actual
improvement, the finite budget, or native/admission failure. No automatic extension,
epsilon addition, bound widening, retry or threshold change. Record actual costs,
rank/conditioning at the evaluated origin and proximity to bounds; a local stop
does not prove global family impossibility. Root owns release/submission after
actual zero-call preparation and a targeted independent review.

## Limits

This is development on consumed geometries and references. No reserved24 outcomes,
new DFT, geometric search, protein mapping adjustment, molecular dynamics,
production promotion or exchange-affinity claim is part of this experiment.
Any passing selected model still requires separately declared numerical and
independent validation before full-system sampling.
