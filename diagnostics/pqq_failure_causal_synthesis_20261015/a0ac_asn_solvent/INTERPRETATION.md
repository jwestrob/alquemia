# What the completed Asn flip establishes

**The local search missed an oxygen-facing configuration that lowers the current
model's La energy. It has not established a relaxed physical basin or a corrected
classification.** The solvent contribution amplifies this configuration change,
but its large contribution is shared by both sample4 sources.

## Actual result and the essential counterexample

Slurm4623 completed all20 new endpoints,80 SCC and40 native component calls;
all60 originals were reused. The full common3→4 comparison retains every source,
and mathematical and operational selections agree. Ca keeps its old adaptive_Ca
winner in both sample4 sources; La changes from adaptive_La to the fixed flip.
Values below are copied from the saved collection, in kcal/mol:

| Source conditioning/sample | La flip minus old La winner: MACE | Frozen ALPB | Total | Common-pool Ca−La change |
|---|---:|---:|---:|---:|
| Ca4, target | −1.460981 | −9.514018 | −10.975000 | +10.975000 |
| La4, counterexample | +8.883471 | −9.556529 | −0.673058 | +0.673058 |

The other eight pools do not change. The different gains in the two sample4
pools therefore do not arise from an exceptional target-only solvent bonus:
their La solvent benefits are almost the same. The important difference is the
vacuum energy already reached by each old adaptive search. Relative to origin,
both sample4 La flips are strongly favorable in MACE (−27.140067 and−27.277089),
but much of that stabilization was already available in the old pools. Those
origin changes must not be presented as new gains over accommodation.

The target's solvent advantage over the old La winner is chiefly electrostatic:
Gelec−9.280011, Ghb−0.186932, Gsasa−0.047075 and Gshift0. The counterexample shows
the same pattern. Passing component parity, diagonal replay and charge-precision
gates establishes faithful numerical evaluation of this model; it does not
independently validate that reaction field or charge distribution for the protein.

## Physical reading and remaining limits

Saved geometry analysis4626 confirms an actual intact neutral-amide rotation:
target OD1–metal4.595274→2.562018Å and ND2–metal2.446070→4.718440Å.
La-conditioned4 changes similarly. This is direct evidence for a donor-orientation
mechanism worth investigating, rather than an atom-name correction or a missing
Asn atom. Earlier Asn modes converged inside their old local domain, which never
included this separately declared180° conformer.

The fixed target also introduces an included HD22–Gly292 backbone-H contact at
1.382191Å; the counterexample has the analogous1.374477Å contact. These are real
represented nonbonded contacts, distinct from the unchanged covalent ND2–H nearest
record. Full-source checks show no new omitted severe-contact obstruction, but
the flipped geometry has not been relaxed or shown to be a minimum. Unfavorable
flips in the other eight sources cannot rule out all alternate relaxed rotamers.

The result remains conditional on the dry oxidized-PQQ3−/neutral-Asn state,
204-atom capped representation, vacuum-MACE Hamiltonian and native ALPB model.
No cofactor-state intervention was tested. In particular, **frozen charge means
frozen against solvent-induced density response, not identical charges across
geometries**: each geometry has its own GFN2 vacuum charge vector. The observed
ALPB change currently combines altered reaction geometry and altered vacuum
charge distribution. Neither lower energy nor improved peer agreement proves a
biological preference or licenses a new band fit.

## One highest-value bounded next test

Use the existing native fixed-charge replay to cross the saved charge vectors
and geometries for **Ca4 and La4, both exchanged metals, old adaptive_La versus
amide_flip**. These two complete geometries are shared by both metals; adaptive_La
is the previous La winner whose roughly9.5kcal solvent difference needs explaining.
Reuse all eight accepted diagonal native results and evaluate only the eight
off-diagonal geometry/charge combinations. No new SCC, MACE, DFT, optimization,
coordinates or state change is needed. This is a diagnostic, not a new scorer.

For each source/metal retain the full two-by-two table G(x_i,q_j), all four native
terms and exact atom-identity mapping. Report both decompositions of the observed
diagonal difference: change geometry at q0 then change charges at x1, and change
charges at x0 then geometry at q1. Retain their mixed term rather than choosing a
path that gives the preferred attribution. “Geometry” here includes interatomic
distances, Born/cavity response and geometric HB/SASA weights; it is not pure
cavity volume. Keep the qualified factory, binary and existing precision limits;
propagate stored-charge bounds and retain unavailable outcomes.

This finite test distinguishes a predominantly geometric reaction-field benefit
from one dependent on geometry-induced vacuum charge redistribution. A substantial
mixed term would show that neither isolated explanation suffices. It will not
establish either contribution as physically correct or prove an Asn basin, but it
directly isolates the dominant untested mechanism without repeating completed
chemistry. It is proposed only; nothing has been prepared or executed here.

## Evidence

- Solvent: `workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_solvent_v1/execution/FINAL_COLLECTION.json`,
  SHA256 `212f8d0ffc82b070a3a7c78ef47e2441acaa76ce061825ffd3662c98f058202f`.
- Geometry: `workspaces/pqq_failure_causal_synthesis_20261015/asn_geometry_response_v1/execution/RESULT.json`,
  SHA256 `d8c1db63f86aef65911cc999460e78ccaa1249410f30db9429423053bbc5d251`.
- Earlier search limits and exact peer inventory: `../a0ac_asn_inventory/REPORT.md`;
  vacuum comparison: `../a0ac_asn_scoring/REPORT.md`.

This interpretation copies saved results only. No numerical analysis, molecular
call, submission, production change or recalibration was performed.
