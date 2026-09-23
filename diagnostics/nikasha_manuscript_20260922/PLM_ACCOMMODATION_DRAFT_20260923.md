# Accommodation contribution to the integrated PLM paper

This is concise current manuscript material. The released static workflow and
archived DFT results remain available; candidate production promotion is separate.

## What the machine-learning component contributes

Nikasha uses a pretrained MACE-OMOL molecular potential to propose bounded,
chemically connected movements of metal-coordinating side chains. This lets the
score account for how each metal accommodates a supplied structure. The tested
improvement is reduced sensitivity to some unfavorable donor arrangements in
predicted structures. The model was not trained on the PLM proteins' presumed
metal preferences, and those preferences remain unknown.

For each fixed chemical state and source context, four independent angular
motions are selected from source-mapped side-chain torsions. Selection uses both
the Ca−La differential force and individual-metal forces, normalized by actual
heavy-atom displacement; native MACE supplies these proposal forces. Both metals
use the same active coordinates, with independent energy-lowering searches.
Angles are bounded by±0.8rad and source-heavy displacements by0.8Å. Caps follow
their physical mappings. The scored PQQ, metal position, proton and water
inventories remain fixed in this branch.

Both metals then evaluate the same finite geometry pool: origin, Ca proposal,
La proposal. At every candidate geometry,

`E_M = E_MACE-OMOL,vac,M + E_native-GFN2,ALPB(water),M − E_native-GFN2,vac,M`.

The paired contrast is `R = selected_E_Ca − selected_E_La`. Mathematical row
minima and operational selections are both retained; operational selection keeps
the origin when its improvement is below0.1modelkcal/mol. These are finite
energy-selected candidates, not equilibrium populations or demonstrated minima
of the composite energy. No entropy, water-reservoir or affinity term is added.

NativeGFN2 scalar inputs use ORCA6.1.1, freshNoAutostart, native mixer,
TolE1e−10Hartree,SmearTemp300 andMaxIter500, one rank per scalar task. Native
MACE-OMOL uses the existing float64 checkpoint, SHA256
`9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`.
The SLSQP objective is shifted from the origin and scaled by0.03674932217934773;
its operative tolerance is1e−8 in that Hartree-equivalent objective, at most200
iterations. Each component is converted once using the retained factors
23.06054783061903kcal/mol per eV and627.509474 per Hartree.

## Supported benchmark statement

Across225 previously inspected alternative structures from25 reference protein
groups, the strongest completed consistent-pocket/accommodation protocol gives
207correct, zero wrong, one inconclusive and17unavailable results. On207 sources
shared with the released scorer, correct calls increase203→206 and both former
wrong calls are corrected. All94 available predefined three-source summaries
remain correct. The original25 calibration inputs and three consumed crystal
controls retain their classes under the candidate's own reference.

A same-solver, same-pocket ablation gives199correct/1wrong/8inconclusive on the
same208 supported sources before accommodation, versus207/0/1 afterward. Each
of the eight repaired source decisions requires an actual contrast change under
the accommodated bands; changing the reference alone does not explain them.
This isolates useful local accommodation beyond consistent pocket preparation.
One C5AXV8 source becomes inconclusive relative to release, and17 preparations
remain unavailable. Replicas are grouped by their actual proteins: these are
developmental structure-transfer observations, not225 independent biological
experiments or prospective validation.

The strong benchmark above fixes fragment membership across the declared ten
reference folds. The practical three-source version instead freezes the union
of complete fragments within a4.3Å polar-contact envelope: the existing3.5Å
contact distance plus the0.8Å admitted displacement. The envelope is a bound
for moving original donor anchors versus omitted fixed atoms, not convergence of
all environment effects. It has its own canonical-only reference. Its completed
primary100-triple transfer gives91correct/9unavailable; one failed solvent cell
blocks three groups in addition to six prior exclusions. A separate, qualified
two-start restart sensitivity restores94correct/6unavailable, with no band or
geometry change. This matches released group fidelity. Individual-source
abstentions increase (four of104prepared pairs versus one with the strict-tenfold
comparator); the separate A0A3Ca3 probe also becomes inconclusive. On the same91
primary complete triples, accommodation repairs two same-context origin
abstentions, but range medians increase7.180→8.941kcal/mol (39decrease/52increase).
Do not describe universal spread reduction or a group-accuracy gain over release.
The original failure and recovery remain separate, and automatic restart is not
part of the generic execution interface.

## Real PLM application example

The practical route completed the exact three existing La-conditioned AF3
sources for each of two consumed, unlabeled PLM proteins. Within the same
prepared context, accommodation reduces raw-score ranges40.138→2.307 and
72.416→17.230modelkcal/mol. Both median predictions remainCa-supported;
one source of the second protein remains inconclusive. The score changes
accompany relief of short extra-Asp contacts: selectedLa distances increase
from1.700–2.120Å to2.219–2.451Å. TwoCa contacts remain short, and7/12 selected
proposals reach the movement boundary. This describes a physical response and
improved source consistency without assigning biochemical truth to these cases.

Archived nativeDFT calculations support the direction and approximate magnitude
of nativeMACE response along the previously specified±0.2rad extra-Asp motions.
Those limited checks do not validate every later optimized candidate. Repeating
the solvent cells with strict stopping changes differential response by at most
0.081kcal/mol; the extra solvent contribution still worsensDFT agreement on the
two tested PLM paths. This numerical explanation is closed. The composite remains
an approximate classifier, and force/entropy qualification is not implied.

The actual two-protein/six-source integration took251allocation seconds on one
H200 with32CPUs:183MACE evaluations and72nativeGFN2 scalar calculations. Separate
preparation functions took16.955s and reused archived protonation. Report this as
batch timing with those exclusions; it is not an isolated single-protein latency
or a matched speed comparison againstDFT.

The candidate export appends results for only these two proteins to the existing
176-protein PLM table, preserving every originalDFT, sequence, gene, phylogenetic,
motif and transcript field. The other174 proteins remain explicitly unscored by
this candidate. Existing normalizations and missing expression data are retained.
OriginalDFT aquo-referencedS and composite rawR have different gauges and must
not be subtracted. Neither score establishesKd, physiological occupancy,+catalytic substrate or a causal environmental effect.

## Exact supporting material

- [Strongest full structural comparison](../strict_native_transfer_20260923/REPORT.md)
  and [matched ablation](../strict_static_ablation_20260923/REPORT.md).
- [Actual practical integration](../pqq_three_source_envelope_execution_20260923/REPORT.md)
  and [reusable experimental command](../pqq_three_source_envelope_api_20260923/COMMANDS.md).
- [All six PLM geometries, components and editable figure](../plm_envelope_response_20260923/REPORT.md).
- [Strict numerical check against archivedDFT](../strict_donor_response_20260923/REPORT.md).
- [Joinable protein export and unchanged biological evidence](../plm_candidate_overlay_20260923/REPORT.md).
- [Full100 primary transfer](../motion_envelope_transfer_20260923/REPORT.md)
  and [separate recovery sensitivity](../native_failed_cell_recovery_20260923/POOLED_SENSITIVITY.md).
