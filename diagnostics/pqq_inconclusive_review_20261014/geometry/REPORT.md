# Inconclusive sources have different geometric limitations

The saved structures do **not** support one shared preparation failure, low-confidence
filter, or optimizer-bound rule that explains all nine inconclusive contexts. The
clearest distinct issues are a displaced Asn in XoxF1 Ca4, limited/uneven donor
accommodation in MMOL1770 Ca0 and several ExaA sources, and unusual metal–PQQ
placement in ExaA two-Ca3. ExaF two-Ca4 and PedE La3/4 have no comparably distinctive
gross defect in these descriptors. This is diagnosis, not a demonstrated remedy.

## Actual execution and coverage

Slurm4038 completed all48 declared source/context cases and their144 saved local
geometries, including all nine targets. Nine exact previously measured raw-source
and confidence rows were reused. Source conditioning/occupancy remain separate;
these48 rows are not48 independent proteins. No new molecular energies or folds.

- [Frozen input](../../../workspaces/pqq_inconclusive_review_20261014/geometry_v1/manifest.json), SHA256 `af1fbc96f08b9d57f0a9eb656ce68415f77f1e2ef4a25b912bbff236f982a2e5`.
- [Actual result](../../../workspaces/pqq_inconclusive_review_20261014/geometry_v1/execution/RESULT.json), SHA256 `b116a1076f14e2d0303b89f31feb42a42311de195e21b985278688e2c30d5ed2`.
- [Receipt](../../../workspaces/pqq_inconclusive_review_20261014/geometry_v1/execution/EXECUTION.json):21.070592 wall seconds,22.401105 process CPU seconds,266680KiB peak RSS, one allocated CPU.

The physical-source reconstruction reproduces the saved target coordinates to
floating-point precision. The targets' archived proposals report successful
optimization and admissible covalent geometry; that does not make them stationary
on the composite scoring energy. Their forces/residuals are from the vacuum MACE
proposal Hamiltonian. Metal/PQQ positions stay fixed in these proposals.

## Target-versus-peer findings

Distances are Å. CN is the existing typed-donor count at3.1Å; a threshold crossing
is a descriptor, not proof of a discrete physical coordination state.

| Target | Actual geometry and motion | Comparison and interpretation |
|---|---|---|
| XoxF1 A0ACD6B9F2 Ca4 | Asn OD1 is4.595 from metal; CN6. Selected Ca/La proposals leave Asn at4.535/4.673, with no recorded boundary flag. | Other Ca samples have Asn2.552–2.672 and CN8. However correctly classified La4 also has distant Asn4.558; missing Asn contact is not a sufficient cause. Metal pLDDT95.99 and metal→PQQ PAE1.50 resemble Ca peers. Current small donor motions cannot recreate the alternate contact basin. |
| XoxF4-2 MMOL1770 Ca0 | Ca proposal shortens Glu OE1 from2.322 to2.151 and Asn OD1 from2.840 to2.404; Glu motion reaches0.8. Selected-mode MACE residual3.401kcal/mol/Å remains. La proposal has no bound and Glu2.386. | Same-Ca peers also show bounds, and the special2.10Å extra-Asp source is **Ca1**, not this target. No unique low confidence or source covalent failure. This is a concrete asymmetric accommodation/model-response question, not a reason to discard the origin. |
| ExaF C5AXV8 two-Ca4 | Both metals select the La proposal; CN7, Glu OE1 remains2.590; no proposal bound, max movements0.328/0.396. | Same-state peers have similar donor geometry and local confidence. Metal pLDDT81.70 and PAE4.683 overlap those peers. No particular geometric failure explains this boundary-near inconclusive. |
| PedE Q88JH5 La3 | CN7 remains7, Glu OE1 increases2.366→2.467/2.581 in Ca/La proposals; no bounds. | PQQ/metal confidence and donor inventory overlap correct La peers. No demonstrated gross defect or failed search termination. |
| PedE Q88JH5 La4 | CN7 remains7; Glu OE1 increases2.304→2.433/2.527. Only La reaches0.8Å bound, residual0.111kcal/mol/Å. | Correct La0 and canonical La2 also have both proposal bounds. Therefore a bound flag does not distinguish this inconclusive; do not enlarge bounds only for failures. |
| ExaA Q9Z4J7 La2 | Both proposals reach0.8Å. Glu OE1 changes2.283→2.813/2.802, with retained CN7. Selected-mode residuals3.681/1.811kcal/mol/Å. | Strong reorientation toward a more balanced carboxylate contact; source confidence overlaps correct La peers. Finite search remains constrained, without demonstrating which further movement would help discrimination. |
| ExaA Q9Z4J7 La3 | Glu OE1 is2.167 at origin; OE2=3.114 leaves CN6 under3.1Å convention. Both proposals reach0.8Å and bring OE1/OE2 near2.65–2.70, CN7. | This is the exact previously flagged short-Glu source. Correct La1 also has short OE1=2.184 and no final boundary flag. Short origin distance alone is not a failure label. |
| ExaA Q9Z4J7 La4 | Both proposals reach0.8Å; Glu OE1 changes2.342→2.853/2.843. Both scoring rows select the La proposal. | Correct canonical La0 also has both bounds. The result cannot be explained just by assigning each metal its own proposal or by identifying a unique malformed source. |
| ExaA Q9Z4J7 two-Ca3 | Unique audited active metal is chainD. Glu OE1=2.099; PQQ N2=2.733 and bound carboxyl O=2.845. Ca/La proposals relieve Glu to2.408/2.493; La hits0.8Å. | Other two-Ca peers have N2=2.429–2.460 and bound carboxyl O=2.562–2.590. Metal pLDDT84.21 and PAE7.283 contrast with peers near91 and2.6. Sidechain-only search leaves this unusual metal/PQQ placement untouched, so it cannot settle that structural hypothesis. |

The target PQQ rings remain nearly planar and pass the previously declared broad
connectivity gate. PQQ itself is fixed during accommodation. Neither this gate
nor predicted cofactor bond lengths independently establishes redox chemistry;
the prior cofactor-state limitation remains. No blanket statement that all these
folds are wrong or that high AF3 confidence guarantees correct ion geometry follows.

## Contact descriptor correction

Do **not** use4038's `closest_mapped_nonbonded_heavy_pairs` or `nonbonded_scope`:
the fixed PQQ internal bonds are absent from the kinematic bond graph and genuine
C–O bonds entered that list. This is an analysis-descriptor defect, not evidence
of a new molecular clash or a change to the scoring calculation. Other fields
above do not use that descriptor. [The separate correction](CONTACT_CORRECTION.md)
preserves4038 and measures inter-residue and omitted-source contacts only.

Slurm4041 completed all48 contact rows with no energy calls:24.905696 wall seconds,
25.811470 process CPU seconds,100120KiB peak RSS. The
[corrected result](../../../workspaces/pqq_inconclusive_review_20261014/geometry_contacts_v2/execution/RESULT.json)
has SHA256 `ac0a8679605b29289d6427a437288bea7fd281ba4cb8d9493d38b47b2c517f11`;
[PEERS.tsv](../../../workspaces/pqq_inconclusive_review_20261014/geometry_contacts_v2/execution/PEERS.tsv)
provides compact recorded motion/confidence/bound descriptors. All mapped-source
coordinate-frame checks permit the omitted-source comparison.

There is a real omitted-packing lead, but it is shared by correct peers. The
Glu213 CG–omitted Leu485 CD2 distances in the **same two-Ca ExaA source state** are:

| Sample | Origin | Ca proposal | La proposal |
|---|---:|---:|---:|
|0, correct|3.8154|3.5971|3.5714|
|1, correct|3.8810|3.6727|3.6463|
|2, correct|3.8322|3.6175|3.6154|
|3, inconclusive|3.8389|3.5280|3.5124|
|4, correct|3.8372|3.5989|3.5707|

The target approaches Leu somewhat more, but both metals make nearly the same
approach. An exterior packing contribution could constrain both movements and
largely cancel from their contrast. Meanwhile its Asp350 OD2 moves **away** from
omitted Phe442 CD2 (3.7160→3.9174/4.0464Å), illustrating why a single nearest
distance cannot substitute for a complete interaction model. La-conditioned
ExaA targets2/3/4 also approach omitted Leu485 in both proposals. This supports a
bounded context-sensitivity test with matched controls, not a claim of a unique
clash, a demonstrated missing energy, or a known sign for a scaffold correction.

Many shortest omitted-source distances are unchanged sidechain CB–neighboring
backbone contacts; their stored displacement is zero. They are not evidence of
accommodation-induced strain. Lists retain the twelve nearest pairs; absence
from that list does not establish absence of an exterior interaction.

## Most informative next test

Use the structural evidence to separate hypotheses rather than add freedom to
every site or tune bounds to these labels. ExaA two-Ca3 is the strongest concrete
case for a small **joint metal/pocket placement** test, with all same-state peers
retained as controls: the present frozen metal/PQQ coordinates prevent the donor
search from addressing its distinctive displacement. XoxF1 Ca4 suggests a
different, source-supported Asn rotamer/contact-basin question, with correct La4
as an essential negative control for the simplistic missing-contact explanation.

Before a broad accommodation expansion, independently qualify the compressed
Ca-favoring donor response in MMOL1770 Ca0 and evaluate any omitted-packing lead
with the full matched-peer comparison above. Bounds and residuals motivate tests;
they do not establish the sign of a missing scaffold energy. ExaF two-Ca4 does not
justify a bespoke intervention merely because it lies close to a class boundary.
No new experiment or production change is executed by this report.
