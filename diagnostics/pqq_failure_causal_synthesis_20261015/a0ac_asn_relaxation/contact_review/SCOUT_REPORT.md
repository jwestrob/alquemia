# Actual saved-contact scout 4859

Job 4859 completed successfully (`0:0`), using one CPU for one scheduler second.
The analyzer reports 0.6623908559558913 s wall, 0.37234734399999997 s process CPU and
44544 KiB peak RSS on `node-64-768g-1`. There were **zero new energy calls**.

All four real-artifact checks passed: Ca0's four saved geometries, preservation
of nine genuinely unexecuted source rows, named identities/paired coordinates,
and exact replay of job 4626's baseline donor distances. All ten sources and
40 candidate rows remain represented. `all_candidate_geometries_available=false`
is correct: the nine unexecuted rows retain both baselines and unavailable
relaxed proposals. It is not a failed interface check.

## Ca-conditioned sample 0 control

Distances below are copied from the saved RESULT, in Å. Both relaxed proposals
have `boundary_flag=false`.

| Saved geometry | OD1–metal | ND2–metal | HD22–Gly292 H | HD21–Gly292 H | HD21–Asn290 H |
|---|---:|---:|---:|---:|---:|
| Origin | 2.5677090956726385 | 4.0271059087140975 | 3.2214360151957067 | 2.813435977590391 | 4.675397523205914 |
| Fixed flip | 3.945817019756787 | 2.543135917236011 | 4.401529720716435 | 4.725961700292389 | 1.1666730115841222 |
| Relaxed Ca | 4.2559038833330005 | 2.6107944297602184 | 4.184329338407598 | 4.649236835195698 | 1.458125039010169 |
| Relaxed La | 4.234907296701167 | 2.8310258224160085 | 3.9386520353350303 | 4.37531006395547 | 1.3598249708327521 |

The control's short fixed-flip **HD21–own-backbone H** contact lengthens under
both saved relaxations. Its shortest H–H pair changes identity to
**HD22–Asp203 HB2**, at 1.3949795862578644 Å for relaxed Ca and
1.1934907303973479 Å for relaxed La (fixed-flip distance 1.5553200223248966 Å).
This illustrates why a generic nearest-contact value cannot establish relief
of the named Gly292 pair. No additional contact threshold is applied.

Ca0 is a control, not sample 4. These results do **not** establish relief of
sample 4's previously identified approximately 1.38 Å HD22–Gly292 H contact.
That question awaits the actual terminal remaining collection and the separate
all-ten analysis. Contact changes alone do not establish an energy mechanism
or an unconstrained physical minimum.

## Immutable evidence

Workspace: `workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_contact_scout_v1/`.

- `CONFIG.json`: `f7b1ddd29a75ae3ae2d278cf7b29309d1bf34805c22529627060c5e762ecd8b5`.
- `execution/RESULT.json`: `e2c95fee91f9ff79e55f47084a4b8e0d15e1c97e855b9caf03d482d9928346a7`.
- `execution/EXECUTION.json`: `4167ef3833bebc0a9d20f9bf2b8938bdd2668b6fbc9e8b7943f9b571c45fc556`.
- `execution/TABLE.csv`: `619b0f08255753f1773f58a0ad527aad32a12fe434aab6cac591583e0cea7363`.
- `completion_event.json` records actual scheduler completion and allocation.

This report reads saved values only. No extra analysis or molecular execution
occurred on the login host; no new job was submitted by this agent.
