# Frozen-f physical force projections: a local direction merits testing

All four actual La/Dy A/B endpoints assemble onto the same 1,891 physical atoms. Among the three prescribed donor/scaffold rotations, **Glu91 carboxylate rotation has the clearest differential load**. This specific 814-atom N83 hinge mainly carries common strain. This does not reject other scaffold motions or demonstrate improved selectivity.

## Actual normalized gradient loads

Units kcal/mol/Å, displacement normalized by the total Euclidean norm over physical atoms. These are gradients, so forces have the opposite sign. Delta means Dy−La. A/B retain their original fixed proton/metal/water inventory; B is the archived +2° Asp85 rotation.

| Physical direction | A La | A Dy | A delta | B La | B Dy | B delta |
|---|---:|---:|---:|---:|---:|---:|
| Metal toward nearest real oxygen |29.15743|20.56186|−8.59557|29.79727|21.44964|−8.34763|
| Asp85 CA–CB sidechain rotation |−1.44761|−1.36051|0.08710|0.75641|1.74747|0.99106|
| Glu91 CB–CG carboxylate rotation |0.94886|3.69257|2.74371|0.88253|3.65044|2.76791|
| N83–CA83 downstream scaffold rotation |−1.91175|−1.97897|−0.06722|−1.90677|−1.97133|−0.06455|

Rotational tangent norms are 3.46104, 3.65200 and 317.12003 Å/radian at A, respectively. Raw generalized loads, electronic/classical splits, moving atom counts, overlaps and exact source identities are in `FROZEN_F_PROJECTION_RESULT.json`; geometry-defined vectors are pinned in the workspace `PROJECTION_DESIGN.json`. Largest normalized overlap is 0.0181, between Glu91 and the scaffold motion. Synthetic caps are never independent motion coordinates. Classical differential contributions are zero for both donor rotations, approximately −2.6e−6 for the hinge, and +0.00617 for metal displacement; the differential signals here are predominantly electronic.

## Common hydrogen strain versus differential chemistry

At A, total common-gradient norms for C-bound H / other H / heavy atoms are 1062.31 / 508.12 / 1495.21 kcal/mol/Å; corresponding differential norms are 0.550 / 0.397 / 41.461. B is similar. Per-Cartesian common RMS values are 24.44 / 15.70 / 28.59, versus differential RMS 0.0126 / 0.0123 / 0.7927. Large hydrogen forces are overwhelmingly common, but heavy-atom common forces are also large: these data do **not** establish that hydrogen alone dominates total strain or causes metal specificity. Shared H-bonded distortions also load heavy anchors. Original known C–H angular defects remain in these inputs.

The next useful development test is an explicitly source-derived small Glu91 rotation after the separately prepared H repair, retaining both endpoint surfaces and the fixed chemistry. The observed loads select this hypothesis; it is not a blind test, affinity result, or optimization prescription.

## Execution and checks

Reused original La A/B classical components after exact coordinate/topology/charge/membership checks. Job1220330 evaluated only missing Dy A/B on OpenMM8.5.1 Reference: two configurations, eight component queries, zero quantum calls, 1 CPU, 2 scheduler seconds = **2 allocated CPU-seconds**, measured evaluator wall1.10488s. Slurm TotalCPU1.420s; reported MaxRSS0 is unavailable evidence, not zero memory. Classical ledger/metal parameters and all receipts remain pinned under `force_assembly_frozen_f_v1`. Earlier La component job cost1 allocated CPU-second; quantum costs belong to the source collection and are not charged as new work here.

All real-input geometry/field/receipt hashes passed assembly checks. Cap forces include both retained/omitted anchor Jacobians; native point-charge gradients scatter to physical source atoms. No additional QM–MM Coulomb, cap FF particles, or C4 terms are introduced. Translation mapping residual≤1.08e−13; total force-sum norms5.44–5.65e−6 kcal/mol/Å. Three additional real-artifact algebra checks passed by direct Python execution. Pytest is absent from this environment; the attempted pytest command did not run tests. These checks establish mapping/accounting, not a numerical derivative gate.

## Limits retained

The finite dry-system energy has no solvent-consistent global response; metal LJ cross parameters remain an unqualified candidate. Earlier classical finite-difference gates failed and are not rewritten. Full hybrid derivatives remain unqualified; no optimization or affinity prediction occurred. Restricted frozen-f valence state is distinct from physical Dy(III)'s sextet. Native-r2SCAN and frozen-f-PBE0 matched La environmental works differ by0.13312 kcal/mol, exceeding the small0.10779 kcal/mol Dy−La Asp85 response; that response is not method-independent evidence. No native-Dy matched result exists. Larger Glu91 load is a reason to investigate a real displacement, not proof the new method is accurate.
