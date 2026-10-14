# ExaA metal placement: eight-cell MACE test4130 running

Preparation4129 passed both sources and all eight geometries on actual Slurm.
The proper PQQ alignment RMSD is0.022439Å; each transferred metal moves0.364311Å.
All nonmetal coordinate lines remain unchanged, both metals use identical
coordinates, and existing bound/source/cap/overlap checks pass. The target's
original donor motion remains at its existing0.8Å bound; no limit was changed.

Scoring4130 now runs on seven shared standard CPUs with the pinned warm float64
MACE worker. Eight new vacuum evaluations; eight exact originals reused. No
optimization, solvent, DFT, new state or class/calibration change. WatchPID1068645
is armed; exact rearm command is scoring_v1/SUBMISSION.json in the workspace.

Read scoring_v1/execution/FINAL_COLLECTION.json and EXECUTION.json when terminal.
Check actual cell/runtime provenance and complete8/8 denominator. Interpret
per-host delta_D under the frozen reciprocal-placement rule in PLAN.md; neither
class labels nor old bands enter this diagnostic. Failures stay unavailable.
Root owns the next decision; no solvent/optimization follow-up is submitted.

Workspace prefix: workspaces/exaa_metal_position_20261014/.
ConfigSHA0874c1088e9eeb6ea29fb0408916b7c1926b5ca90396772392c2657ef4ba2154.
Preparation took3scheduler seconds on1CPU; model cost remains unmeasured until
actual receipts. EXECUTION.json links pinned commands/manifests/results.
Previous eight-cell DFT reference is complete; do not rerun its delayed wakes.
