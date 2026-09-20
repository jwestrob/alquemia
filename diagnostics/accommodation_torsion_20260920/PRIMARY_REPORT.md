# Primary torsion profiles: DFT validation running

The affordable surface predicts a strong metal-dependent response to the two
previously identified compressed PLM donors. All 64 native OMOL scalar evaluations
completed; 126/128 GFN2 vacuum/ALPB calls converged, yielding 63/64 complete
composite endpoints. Both control origins reproduce archived MACE and solvent
energies exactly. This is development evidence, not a corrected classifier.

## Fixed extra-Asp positive rotation

Work relative to the same metal's original source, kcal/mol; angle +0.2 rad.

| Structure | Ca work | La work | Ca−La work |
|---|---:|---:|---:|
| 4MAE XoxF control | +1.893 | +0.945 | +0.948 |
| PLM 83440678cbbd658047c9 | −23.019 | −49.976 | +26.957 |
| PLM 07ab500e3df76b30d71c | −12.619 | −30.293 | +17.674 |

The PLM changes relieve their pre-existing short metal–oxygen contacts. The
reverse rotations compress them further and cost substantial energy. Both
PLM profiles continue downhill toward +0.4 rad; no minimum is established.
No profile angle is selected as a new production structure or prediction.

This large response is chiefly already present in vacuum OMOL. Its Ca−La work
is +28.535/+18.917 in the two PLM structures; the solvent-transfer contributions
are −1.578/−1.242. Thus this is evidence for a targeted physical accommodation
hypothesis, not evidence that ALPB alone creates a useful signal. Full component
and signed profiles, including anchor-Glu motions, are retained in
`workspaces/accommodation_torsion_20260920/primary_result_v1.json`.

The solvent term is not smooth everywhere in this small sample: 4MAE La's
anchor-Glu −0.4 rad work contains a −3.392 kcal/mol solvent-transfer change,
compared with −0.208 at −0.2, while OMOL work is smooth (+0.923/+0.252).
This is a component observation, not an established diagnosis of a new basin
or solver branch. It remains outside the frozen native validation selection and
is not used to define an optimized score or harmonic curvature model.

## Explicit failure and pending test

The 1H4I La anchor-Glu −0.2 rad endpoint fails native GFN2 SCF in both vacuum
and ALPB. Both actual failure outputs remain preserved; no last-iteration energy
or baseline substitute is used. Unrelated MPI cleanup messages follow the
explicit SCF nonconvergence. No identical retry or parameter change was launched.
The 1H4I profile remains incomplete, with the full denominator retained.

All 16 **preselected** DFT validation endpoints have valid cheap counterparts;
the failed 1H4I point was never part of that frozen selection. Job 1203771 is
running the original scope: extra-Asp ±0.2 in 4MAE and both PLM cases, both
metals, plus four PLM origins. Two exact archived 4MAE context origins are reused.
No accuracy claim awaits the desired sign: native work and model errors will be
reported regardless of outcome. PLM truth labels remain unknown.

## Cost so far

MACE job 1203745: 13 allocated seconds, 32 CPUs and one H200 = 416 core-seconds
and 13 GPU-seconds. Warm process time was 5.696 s for 64 evaluations, including
1.175 s model loading; peak allocated VRAM 3,681,531,904 bytes. All four archived
control energies reproduce exactly, validating this warm scalar execution.

GFN2 job 1203746: 388 allocated seconds on 64 CPUs = 24,832 core-seconds,
including the two failures. Recorded batch MaxRSS 7,061,924 KiB. Total primary
allocation cost: **25,248 core-seconds and 13 GPU-seconds**. Local preparation,
analysis and tests are outside these allocation receipts. DFT cost is pending.

Original DFT/default scoring and all prior experiments remain preserved.
