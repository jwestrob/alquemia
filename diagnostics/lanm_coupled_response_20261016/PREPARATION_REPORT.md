# Coupled motions are physically admissible; finite electronic benefit untested

Slurm12460 completed on16October2026. Both own-metal proposals for both Hans
sources pass the unchanged geometry gates. The complete existing classical
ledger was evaluated for all eight candidate/metal combinations. No quantum
calculation, optimization trajectory or new preference score was produced.

The eight-coordinate span has rank8 in both sources. It combines Asp85/Asp87
whole-donor torsions, the native N83-CA83 backbone hinge and target translation.
The largest realized atomic displacements are0.0123–0.0140angstrom, below the
declared0.02angstrom limit; torsions are below0.09degree. Bond/angle and signed
volume changes are at floating-point roundoff, with no new severe-contact or
inversion failure. Physical, cap, field and omitted-anchor maps reconstruct.

| Source | Proposal selected by | Own-gradient linear work | Actual own-metal classical work |
|---|---|---:|---:|
| Hans8DQ2 | La | -0.431567 | -0.000694 |
| Hans8DQ2 | Dy | -0.376314 | -0.001024 |
| Hans8FNR | La | -0.624898 | -0.000224 |
| Hans8FNR | Dy | -0.491166 | -0.000591 |

Work units are kcal/mol. The first column of energies is a **local gradient
prediction**, not an evaluated finite electronic/total energy. At these tiny
displacements, changes in the evaluated classical ledger are small. This does
not establish that scaffold resistance stays small over a larger accommodation,
that the projected electronic work is accurate, or that the source reversal is
resolved. Individual target translation loads exceed the backbone-hinge load;
the current span is a limited controlled movement, not whole-protein relaxation.

## Actual cost and records

13 scheduler seconds on1CPU (13 allocated CPU-seconds),10.346 measured CPU-seconds;
batch MaxRSS324328KiB. Application elapsed12.777207seconds. Eight classical
energy/gradient cells,32 component queries, four reused classical origins,
zero new QM/MACE calls, zero integration steps. All four proposals are retained.

Workspace: `workspaces/lanm_coupled_response_20261016/preparation_v1/`.
`execution/RESULT.json` SHA256
`08d8b0fb445107c7b17964e22281af85a923462bbd9bad15f64e5db26f72a2bc`;
`PREPARE.py` SHA256
`12eec670c34ef60004ec3c63784eb12ee5769b213e0903883901ec2d28fd3933`.
Each source's PROPOSALS/MODES files link complete geometry, kinematics, candidate
inputs and both evaluated classical rows. SUBMISSION.json records the exact
command and completion watcher. No baseline or old output changed.

## Next scientific action

Prioritize the user's suggested monomer benchmark, keeping Hans as a consumed
source/assembly-sensitive control. Recover exact Mex and LBT3 evidence and usable
whole-scaffold representations before declaring its finite response comparison.
LBT3's opposite preference makes it a useful peptide control, not a substitute
for validating LanM. Do not spend eight new QM evaluations on the present tiny
Hans motions merely because their geometry passes. Keep them as reusable
diagnostic candidates; the overall profiler goal remains active and unmet.
