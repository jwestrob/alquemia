# Possible throughput improvement: selected metal Hessian

Read-only source inspection, 2026-09-18; **no new inference or derivative test**.
This is a later engineering option if the current transfer result merits it.
It does not change HYBRID_GGR_TRANSFER_PLAN.md or its running grids.

The installed pinned MACE implementation exposes the coordinate graph with
`compute_force=False, compute_hessian=False`: ScaleShiftMACE.forward passes the
same positions tensor through prepare_graph, and get_outputs leaves forces
uncomputed in this configuration. With autograd enabled, a potential route is
one energy forward, grad(E, positions, create_graph=True), then three backward
calls for the selected metal's x/y/z gradient components, retaining only that
metal's three coordinates from each result. This computes a3x3 projected
learned Hessian, not a dense whole-protein or quantum Hessian.

Do not enable the installed generic `compute_hessian=True`: its
compute_hessians_vmap constructs an identity with3N rows and requests the whole
matrix. The loop fallback still iterates all3N forces. Neither matches the
small physical coordinate subspace needed here.

Source inspected in the existing environment:
`workspaces/mace_hybrid_20260916/software_v1/venv/lib/python3.11/site-packages/mace/modules/models.py`
(ScaleShiftMACE.forward) and `modules/utils.py` (prepare_graph, get_outputs,
compute_hessians_vmap/loop). Existing `scripts/mace_omol_gradients.py` uses
non-reentrant edge/product checkpointing. Its second-derivative support and
memory/runtime remain untested. A working first derivative does not establish
that the second derivative fits or is correct.

A future contained engineering check should compare the selected matrix and
sphere prediction against the already executed full/fine grids, with unchanged
weights, charges, mask, geometry and radius. Count all recomputation and memory.
It needs no new quantum calculation to establish numerical equivalence. Do not
claim speedup, physical curvature validation, or replace the current grids until
actual second-derivative outputs support those claims. No training is involved.
