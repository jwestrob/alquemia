# Physical donor response: saved-gradient result

The source-defined donor motions are implemented and pass geometry/gradient
checks on all 16 actual endpoint states. They reveal substantial donor
sensitivity in the failed GGR transfer cases. **No improved score is claimed:**
this step measures derivatives, with no new energy evaluation or curvature.

## Results

Values below are first-order changes for a common probe with maximum
heavy-atom motion 0.02 A, in kcal-equivalent. They are neither evaluated
energy changes nor relaxation corrections. The column is the largest absolute
value across all declared modes of that type. All cases are consumed
development, representing two biological groups.

| Representation | Metal | Donor |
|---|---:|---:|
| ALPHA_1F6S | 1.222681 | 0.229227 |
| ALPHA_6IP9 | 1.463305 | 1.047183 |
| GGR_1GLG_connected | 0.485616 | 0.422372 |
| GGR_1GLG_extended | 0.499820 | 0.487979 |
| GGR_2FVY_connected | 0.138241 | 0.545992 |
| GGR_2FVY_extended | 0.541307 | 0.573787 |
| GGR_2FW0_connected | 0.050467 | 0.426887 |
| GGR_2FW0_extended | 0.682592 | 0.474097 |

The largest donor mode is Gln142 chi1 in every GGR representation, and Asp82
chi1 in both alpha-lactalbumin structures. In connected 2FW0/2FVY, donor
sensitivity exceeds metal sensitivity. This makes coupled donor/metal response
a useful next question, without establishing the sign of a relaxation benefit.
Both endpoint gradients and their curvature matter.

Maximum differences between connected/extended first-order probes are
0.207731 (1GLG), 0.709610 (2FW0), and 0.679548 (2FVY). These are derivative
representation differences, **not** finite-energy partition tests.

## Method and actual checks

Protocol `matched_hybrid_physical_donor_coordinates_v1`; 16 coordinates per
GGR and 11 per alpha site. Every declared donor chi bond and backbone peptide
crankshaft is included. All explicit waters and exterior atoms remain fixed.
Both GGR representations and both metals share the same full physical map.
Source atoms and cap anchors retain their actual identities. See [plan](PLAN.md).

The projected hybrid gradient is Jcore^T(gDFT-gTcore)+Jfull^T gTfull, with
Ca-minus-La sign. Metal derivatives retain kcal/A and rotations kcal/radian;
the displayed probe scales avoid comparing raw mixed-unit norms.

All finite geometry Jacobian errors are <=7.052e-9 A per coordinate unit
(tolerance 1e-7); bond-length changes <=2.088e-14 A (tolerance 1e-9).
These finite geometry checks are not numerical DFT gradients. Actual gradient
projection agrees with the independent source/cap chain rule. Three real-fixture
tests passed in 18.098 s, zero skips, including combined motions and rejection
of an explicitly corrupted real source bond.

## Cost and next work

Zero new DFT/MACE/solver calls. Local preparation: 28.044262 s wall /
26.487931 s CPU; report: 0.052886 s wall. Tests are additional; earlier gradient
calculations retain their original costs. No Slurm job was launched.

Curvature, relaxation, entropy and calibrated scores remain unavailable
(`response_model_not_validated`). The next model must couple physical motions
and test actual displaced native energies across all three GGR structures,
rather than selecting the previously successful 1GLG case. Production unchanged.

## Reproduce

Run from the repository root, choosing a fresh output directory:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_site_response.py prepare --config workspaces/mace_site_response_20260918/config.json --agreement diagnostics/mace_site_response_20260918/PLAN.md --output workspaces/mace_site_response_20260918/prepared_replay_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_site_response.py report --prepared workspaces/mace_site_response_20260918/prepared_v1/preparation.json --output workspaces/mace_site_response_20260918/report_replay_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_site_coordinates.py -v
```

[Compact result and exact source hashes](RESULT.json). Full per-mode results,
source mappings, gradient components and pinned implementations are retained
in `workspaces/mace_site_response_20260918/prepared_v1` and `report_v1`.
