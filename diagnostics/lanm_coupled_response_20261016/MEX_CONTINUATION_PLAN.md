# Finite Mex continuation from the admitted seven-coordinate pilot

Declared 16 October 2026, before continuation calls. Root owns submission.

## Question and unchanged physical model

Does further bounded own-energy motion in the same seven coordinates materially
change the paired accommodation work after pilot 12717? Its own La and Dy works
are −0.252489 and −0.354867 kcal/mol, while substantial target and hinge loads
remain. This is EF3 response with three fixed Nd spectators, not protein affinity
or a stationary minimum. No outcome sign selects the permitted continuation.

Keep the exact 2049-atom physical inventory, 209-atom local region and four caps,
1840 shifted external charges, ff19SB/TIP3P plus formal +3 global charge vector,
and all source mappings from the admitted pilot. The Hamiltonian remains
MACE vacuum(Q) + GFN2 embedded(Q;MM) − GFN2 isolated(Q) + existing MM ledger +
native fixed-charge global ALPB. Local physical La/Dy multiplicities are 1/6;
published f-in-core GFN2 uses effective UHF 0. No whole-protein SCC, global MACE
spin assumption, new solvent prescription, charge fit or force correction.

## Immutable starts and basis

Read `workspaces/lanm_coupled_response_20261016/mex_response_v1/execution/RESULT.json`,
SHA256 `2c5c50010f502ee0a09262cefa8934b5237716dc7290fe6a8cbc187a103eef32`.
Require its passed four-cell result and both own descents. Reuse exact accepted
La_proposal/La and Dy_proposal/Dy coordinates, native component energies and full
assembled gradients. No origin or warm-start molecular replay.

Pin the pilot MODES.json (`b2f227d576fc9b1fb2f855408f5591d228230449f8cc67e0a4e70cceed1a2d1d`)
and PROPOSALS.json (`b7dc055a2853c0170c15339fb3e6f0c0b0e1bab44c7a8c4c65ac59194dc76c84`).
Retain source-defined order and exact moving membership: Leu83 N–CA complete
hinge, Asp84 CA–CB, Asp86 CA–CB, Glu95 CA–CB, then target x/y/z. Reconstruct
coordinates from the immutable original physical source and total q at every
trial. Original source coordinates and severe-contact baselines are never rebased
to a pilot or continuation endpoint. Waters, spectators and fixed complement stay
fixed. Preserve all pilot protected groups, bond, angle and chirality checks.

Use original-source column normalization and the symmetric inverse square root
of its normalized Gram matrix (relative cutoff 1e−10). Require retained rank 7
and orthonormality error ≤1e−10. Optimize in these fixed Cartesian-length metric
coordinates; update the exact current torsion Jacobian for every force projection.
The start q and gradient are copied from 12717; whitening may introduce only
float64 coordinate-reconstruction roundoff ≤1e−12 Å.

## Frozen finite search

Use the original LBT own-energy L-BFGS/Armijo pattern: empty initial history,
history length 8, descent reset to negative current gradient if needed, Armijo
coefficient 1e−4 and trial multipliers 1, 1/2, 1/4, 1/8. No alternate direction,
constrained fallback or retry after the four trials. Curvature acceptance remains
s·y > 1e−10 ||s|| ||y|| and s·y > 1e−16. Stopping thresholds in the represented
Cartesian span are maximum atom gradient ≤0.1 and Cartesian RMS ≤0.03
kcal/mol/Å; a bound or exhausted trial budget does not imply convergence.

Each metal gets at most 12 accepted steps and 24 new own composite evaluations.
Each trial is bounded by 0.02 Å maximum physical atom displacement and 0.5°
maximum torsion increment. Use a conservative sum-of-arcs scale and then check
actual physical displacement. Relative to the original source, all atoms and
the target remain within 0.20 Å and each total torsion within 0.20 radian.
These envelopes transfer the existing LBT finite-search rule before inspecting
the continuation outcome. Keep the pilot's stricter bond tolerance 1e−10 Å,
angle cosine 1e−10, protected distance 1e−9 Å, fixed complement 1e−12 Å and
severe-contact thresholds (heavy 1.0 Å, H 0.55 Å, worsening allowance 1e−7 Å).
Geometry failures cost no energy call. An admitted exact geometry/state cache
hit costs no new molecular call. A backend failure stops that search as
unavailable and preserves its last accepted point and all attempted receipts.

## Execution and paired denominator

The adapter first performs zero-call real-artifact preflight for both starts,
source graph/membership, current cap/field maps, saved gradients and all states.
Run La scout first. Only an admitted scout with all guards and own descent may
release Dy under the identical source, plan and settings. Do not change bounds,
basis or search settings after the scout. Native MACE stays warm within a stage;
two one-thread local GFN calls can run concurrently. MM uses the pinned original
OpenMM Reference environment through an owned CPU subprocess, and unchanged
native ALPB uses its qualified binary. Root uses one shared GPU/12-CPU allocation
per molecular stage, no memory/time overrides. Scheduler receipts are authoritative.

After both searches, retain original common3 separately. Build common5 from the
original origin, both pilot proposals and both admitted terminal continuation
geometries; both metals see the same complete menu. Exact coordinate deduplication
retains all declared aliases. Reuse the six original cells, both terminal own
cells and exact matching admitted cache entries; allow at most two new terminal
opposite-metal composite cells. Rejected trial geometries are never extra pool
candidates. Keep both terminal slots in the denominator if either search fails;
then common5 is unavailable while original common3 remains valid.

Maximum new calls across both searches and pool: 50 MACE, 100 local GFN2,
50 MM configurations/200 MM group queries and 50 ALPB. Staging the unused
opposite-metal task is not a molecular call. Record actual attempted/completed
calls, model loads, application wall/CPU/GPU intervals and scheduler costs.
The pilot's roughly 30-second four-cell application suggests minutes for this
bounded continuation; serial MM startup and search depth make this only a rough
resource extrapolation, not a measured saving or promised duration.

## Report and limits

Report every accepted/rejected step, original-source and incremental own work,
component work, terminal residual loads, stopping reason, geometry bounds and
availability. Preserve original common3 and complete common5 selections and
their paired Dy−La accommodation change. No metal-specific raw energy offset is
interpreted as affinity. No new basis, larger envelope, further optimization,
parameter change or cross-protein conclusion follows automatically.

Prior implementation pins: pilot `f11777c7c49706c6798f4534ed7bf70b908f1fc1dbf8d815de0d39e6283b35d0`;
LBT search pattern `ee1db81bdf3db90c303f824d7a55b52b00c2bb64f2422903aed60c25931d8ca1`.
All numerical preparation, validation and evaluation run inside Slurm.
