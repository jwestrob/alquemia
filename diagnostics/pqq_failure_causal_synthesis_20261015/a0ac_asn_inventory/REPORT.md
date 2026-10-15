# XoxF1 sample4 has a distinct terminal-amide orientation

**The evidence supports a finite Asn290 amide-rotamer question, not a missing
donor or missing search mode.** Both Ca-conditioned sample4 and La-conditioned
sample4 point the neutral amide **ND2** toward the metal while OD1 is distant.
The other eight same-protein sources have the opposite contact pattern. Current
searches already vary both Asn angles and converge inside their local domain;
they do not reach the alternative orientation. This is not proof that AF3
misassigned atoms, that the amide state is wrong, or that changing it improves
classification. Correctly called La4 is an essential counterexample.

This inventory reads saved artifacts only: no geometry calculation, construction,
scientific imports, molecular calls or submissions. [evidence pointer](EVIDENCE_POINTER.json)
pins all ten sources, all sixty current common-pool cells and their saved Cartesian
forces, all sixty accepted frozen-charge native-ALPB results, maps/states and
the earlier tests. SHA256: `c91b152b1fb8b548395c41aac72fb2370b5ba73d015abbfce14767ac2787d4d8`.

## Exact peer evidence

These are copied descriptors from actual Slurm4038, not new measurements:
[RESULT.json](../../../workspaces/pqq_inconclusive_review_20261014/geometry_v1/execution/RESULT.json),
SHA `b116a1076f14e2d0303b89f31feb42a42311de195e21b985278688e2c30d5ed2`.
All ten declared same-protein sources are included; conditioning is not biological
class. Distances are Å; chi2 uses the saved CA–CB–CG–OD1 convention.

| Conditioning/sample | Metal–OD1 | Metal–ND2 | Asn chi2, degrees |
|---|---:|---:|---:|
| Ca0 | 2.568 | 4.027 | −62.077 |
| Ca1 | 2.612 | 4.073 | −62.343 |
| Ca2 | 2.552 | 4.024 | −62.621 |
| Ca3 | 2.672 | 4.102 | −61.500 |
| **Ca4, target** | **4.595** | **2.446** | **+78.781** |
| La0 | 2.561 | 4.004 | −63.681 |
| La1 | 2.597 | 4.064 | −64.557 |
| La2 | 2.542 | 4.004 | −65.032 |
| La3 | 2.637 | 4.094 | −64.471 |
| **La4, correct counterexample** | **4.558** | **2.414** | **+75.716** |

The saved source confidence is high even for Ca4; it does not establish chemical
correctness. The typed coordination count excludes amide ND2 as a donor, so a
short ND2 distance must not be relabeled as recovered oxygen coordination.

## Current representation and motion are explicit

Every current peer has the same 204-atom envelope identity
`a0acd6b9f2-pqq-la_model__envelope_2eda778295ffed58`: dry oxidized-PQQ3−,
neutral Asn, Ca−4/La−3 singlets, one source metal. This differs from the historical
180-atom search pilots. Asn belongs to the original80-atom core, not merely an
added environment fragment. The expanded source graph contains its complete
backbone, CB–CG, CG–OD1/ND2 and both ND2-bound hydrogens, with Gly289 and Pro291
peptide connectivity. It is not omitted by a distance cutoff.

For target Ca4, current QM indices are CB58, CG59, HD21/HD22=64/65, ND2=67,
OD1=69 (zero based). Physical-kinematic indices differ: the actual chi2 axis is
CB64–CG65 and its moving group is `[65,70,71,73,75]`. Use the named source map,
not these numbers as a universal recipe. Current caps terminate elsewhere; no
synthetic Asn cap needs independent motion. The state explicitly binds both H
atoms to ND2. Full states/maps and exact per-peer index mappings are in EVIDENCE.

The current pilot34 task selects Asn290 chi2, Asn290 chi1, Glu204 chi1 and Asp347
chi1. At q0 the stored normalized Asn chi2 gradients are −40.663/−62.512
kcal/mol/Å (Ca/La). After search, Asn chi1/chi2 residuals are approximately
−0.0070/+0.0027 for Ca and +0.00049/−0.00138 for La. Maximum heavy displacements
are 0.317892/0.487031Å, with no boundary flag. OD1 remains 4.535/4.673Å away.
These are vacuum-MACE gradients, not composite gradients or error estimates.

Exact current proposal manifest:
`workspaces/motion_envelope_20260923/pilot34_searches_v1/manifest.json`.
Exact target source preparation:
`workspaces/motion_envelope_20260923/run_v1/prepared/a0acd6b9f2-pqq-la_model__envelope_2eda778295ffed58/cases/a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4/context_preparation.json`,
SHA `fa4eb04da783468433f98642c84497dbff7086c34b881d998e547a75914a8859`.

## What has already been tried

| Actual earlier test on this source | Result and limit |
|---|---|
| Asn290 chi2 basin scan | 65 nodes over exactly ±0.7048068578544397rad, every other coordinate fixed. Width contribution −0.083615kcal/mol; adding the grid retained the existing adaptive winners. This local scan did not sample the alternative amide orientation. |
| Reference-derived starts | Target Asn chi2 changes −2.463862/−2.544326rad, with maximum heavy motion2.423669/2.291545Å. Both were rejected before energy for exceeding the frozen domain. The analogous La4 starts were also rejected. Therefore no alternative target basin was scored by this branch. |
| Six-angle warm expansion | Actual Ca4 contrast change only +0.101593kcal/mol; extra modes were Glu204 chi1/chi2. The branch is closed. It does not demonstrate an Asn flip. |
| Generic collective scaffold | Both sample4 sources were tested; all searches reached the physical guard, and protein-only relaxation damaged relative ordering. Do not repeat this generic proposal. |
| Solvent-guided probes | Both sample4 sources retained their prior pooled calls; no new useful basin established. |

These tests used their own older representations and native-ORCA GFN2
ALPB-minus-vacuum expression. Their energies cannot populate the current204-atom
standalone-xTB/tblite-qualified, frozen-vacuum-charge candidate by relabeling.
Exact raw collections and reports are pinned in EVIDENCE. In particular:
`workspaces/local_basin_breadth_20260922/grid_v1/design.json`,
`workspaces/local_basin_breadth_20260922/RESULT_v2.json`, and
`workspaces/structure_informed_starts_20260922/prepared_v3/PREPARATION.json`.

## One justified next test, not yet prepared

First admit a **separately declared finite180° terminal-amide rotation about each
host's real CB–CG axis**, carrying OD1, ND2 and ND2's two attached physical H atoms
together. Keep CG on the axis and every other coordinate, atom identity, proton,
cap and chemical state unchanged. This physically exchanges the amide's spatial
orientation; it is **not atom renaming or an O/N coordinate-label swap**. The
angle is a fixed amide-flip hypothesis, not fitted to a desired score or to the
observed approximately140° peer separation.

The exact graph supports this construction. The full protonated target source
already exists at `workspaces/accommodation_goal_20260920/folds_v1/source_preparation/cases/a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4/protonated.pdb`,
SHA `0a399c6378e6612d697c1a2b41d94eb65e5550ccf2bdf26d4e2dbbdb2b69d67b`.
All peers have matching source-to-context identities. Slurm4041 also established
the mapped full-source frame for omitted-protein contacts. There is **no prepared
full-source flip, admitted alternate geometry or scored flip** in this inventory.
Full-source contacts, covalent/rigid-amide invariants and unchanged cap geometry
must be checked on Slurm before any energy; source-heavy contact maps alone do
not certify hydrogen contacts.

Retain all ten origins in the geometric feasibility denominator. If root declares
the finite energy stage, both metals should see the same admitted flipped
geometry, with the three existing candidates retained separately in the common
pool. La4 and the first same-conditioning peer Ca0 must remain visible, alongside
the other peers. The180° test lies outside the old0.8rad/0.8Å local-search domain;
record a new discrete-rotamer protocol, **do not loosen or retroactively pass the
old gate**. Geometry failure closes that endpoint. Even a favorable shift would
only support an orientation sensitivity, not establish the observed rotamer,
thermodynamic populations or classifier improvement.
