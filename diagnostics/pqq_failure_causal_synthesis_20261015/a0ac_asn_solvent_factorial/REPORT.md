# Geometry gives the solvent benefit; charge redistribution offsets it

**Within the current native ALPB model, changing the complete saved geometry is
favorable at either saved charge distribution. Replacing the old vacuum charges
with the flip's vacuum charges offsets that benefit.** This behavior is shared
by both sample4 sources. A large mixed term makes the magnitude assigned to
each step depend on the decomposition path.

Slurm4813 completed all8 new off-diagonal evaluations. Eight genuine archived
diagonals were revalidated and reused, yielding all16 cells and all four complete
source/metal factorials. No SCC, MACE, DFT, geometry optimization or new diagonal
evaluation was performed.

## Both paths, with the full denominator

Let x0 be the complete saved `adaptive_La` geometry and x1 the complete
origin-derived `amide_flip`; q0 and q1 are their respective actual vacuum GFN2
charges. The geometry step therefore changes more than an isolated Asn atom or
bond: it compares the two complete204-atom arrangements used in the prior result.

All values below are actual saved totals in kcal/mol, rounded for display.
Unrounded four-term values, matrix cells and bounds are in [RESULT.json](RESULT.json).

| Conditioning / evaluated metal | Diagonal G11−G00 | Geometry at q0 | Charge at x1 | Charge at x0 | Geometry at q1 | Mixed term |
|---|---:|---:|---:|---:|---:|---:|
| Ca4 / Ca | −7.763898 | −22.924209 | +15.160311 | +2.375675 | −10.139573 | +12.784636 |
| Ca4 / La | −9.514018 | −25.073250 | +15.559232 | +1.850556 | −11.364574 | +13.708676 |
| La4 / Ca | −6.747875 | −22.955000 | +16.207125 | +3.190905 | −9.938780 | +13.016219 |
| La4 / La | −9.556529 | −25.479253 | +15.922724 | +2.250236 | −11.806765 | +13.672488 |

The geometry-first path is `(G10−G00)+(G11−G10)`; the charge-first path is
`(G01−G00)+(G11−G01)`. Every component closes exactly in the saved arithmetic.
The mixed term `G11−G10−G01+G00` is retained instead of choosing one path as a
unique attribution. In particular, the roughly9.5kcal/mol La solvent benefit is
not created by a favorable change of vacuum charge vector: that charge change
is unfavorable at both geometries. It attenuates a larger favorable geometric
reaction-field contribution under this model.

The Ca rows use the same `adaptive_La`→flip comparison as the La rows. They do
not represent replacement of the Ca winner, which remained `adaptive_Ca` in
the completed common4 selection. This factorial is not a new pool score.

## All four terms retained

The actual diagonal changes are predominantly Gelec; Ghb and Gsasa contribute
smaller favorable terms, while Gshift is unchanged:

| Conditioning / evaluated metal | Gelec | Ghb | Gsasa | Gshift |
|---|---:|---:|---:|---:|
| Ca4 / Ca | −7.538561 | −0.175130 | −0.050207 | 0 |
| Ca4 / La | −9.280011 | −0.186932 | −0.047075 | 0 |
| La4 / Ca | −6.543128 | −0.168731 | −0.036016 | 0 |
| La4 / La | −9.301311 | −0.214263 | −0.040956 | 0 |

For target La, the geometry-first Gelec change is−24.96752660567937kcal/mol
and the following charge step is+15.687515921262941. The reverse path gives
+2.026656944711151 then−11.30666762912758. The Gelec mixed term is
+13.66085897655179; the Ghb mixed term is+0.04781679309358578. Gsasa and
Gshift have zero mixed terms. Both other-source and Ca-metal decompositions
remain fully recorded, without selecting a favorable path or dropping a term.

## Scientific qualification and precision limits

Before dispatch, all eight original diagonal executions passed existing
component-v2 and native qualification. Their actual command/binary, XYZ/charge
hashes and raw printed terms/bounds matched the saved records. The crossed
inputs retain exact atom order, physical maps, state and charge-source identity.
The new cavity and Born-matrix hashes match their same-geometry diagonal;
charge-independent terms and printed charge sums pass the unchanged gates.

All16 individual cells satisfy the original per-term precision requirements.
All diagonal and two-cell path differences have informative aggregate bounds.
The two **Ca-metal mixed terms do not pass** the unchanged1e−6Eh aggregate
precision ceiling: their conservative bounds are1.084790250692721e−6Eh
(Ca4/Ca) and1.0880723902754505e−6Eh (La4/Ca); their Gelec informative flags are
also false. Both La-metal mixed terms remain informative, with bounds
8.216517974053357e−7 and8.260501808537134e−7Eh. These flags and all outcomes
are retained; no precision threshold was changed. Historical absolute-v1
vacuum/water parity failures remain distinct from accepted component checks.

## What this establishes—and what remains conditional

The result narrows the mechanism within the existing ALPB expression: the
favorable response survives either fixed charge vector, while charge adaptation
opposes it and strongly couples to geometry. “Geometry” includes interatomic
distances, Born/reaction geometry and geometric HB/SASA weights; this is not
an isolated cavity-volume test. Crossed distributions are diagnostic vectors,
not new self-consistent electronic states.

The calculation does not validate this reaction field for the intact protein,
identify a unique atom's causal contribution, prove a relaxed Asn basin, or
establish biological preference. The dry oxidized-PQQ/capped representation
and the included unrelaxed sample4 H contacts remain unchanged. No new bands,
class decisions, fitted corrections or production changes are introduced.

## Actual cost and evidence

Job4813 completed0:0 on `node-64-768g-2`: **4 scheduler seconds,1 CPU**, TotalCPU
00:02.623, batch MaxRSS177,188KiB. The application receipt measures
4.417345321038738s and2.018821777 parent CPU-s; child-native costs are separately
pinned. The application and scheduler timers use different boundaries/resolution.
Actual count:8 new native evaluations,8 reused diagonals and zero new diagonal
calls. This cost excludes all earlier SCC/ML and diagonal generation.

The authoritative result is
`workspaces/pqq_failure_causal_synthesis_20261015/asn_solvent_factorial_v1/execution/FINAL_COLLECTION.json`.
[RESULT.json](RESULT.json) pins it, the manifest, actual preflight, partial record,
execution and submission. [READY.md](READY.md) records the completed command and
independent review. Reporting copied existing Slurm values without new scientific
calculation or reruns.
