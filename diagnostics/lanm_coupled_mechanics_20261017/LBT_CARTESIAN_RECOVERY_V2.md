# Cartesian saved-parent replay recovery

Preflight12954 failed before molecular calls in the historical matched-response
loader's exact ten-mode lineage comparison. Saved-only diagnostic12956 reproduced
that failure and isolated its sole false conjunct: exact float-list equality
between the archived whitening transform and its regenerated eigensystem result.
Every state/source/point/active-mode identity check passed.

The maximum transform difference was 1.1102230246251565e-16 for7CCO and
1.942890293094024e-16 for7CCN, for both metals. The archived/regenerated original
whitening residuals were respectively1.11e-15/6.66e-16 for7CCO and
1.33e-15/1.22e-15 for7CCN, inside the unchanged1e-10 parent criterion.
Exact immutable failure and diagnostic pins are in
`LBT_CARTESIAN_RECOVERY_V2_EVIDENCE.json`.

`LBT_CARTESIAN_V2.py` changes only saved-parent loading and recovery provenance.
`LBT_CARTESIAN_REPLAY_V2.py` contains explicit copies of the two required parent
loader bodies; historical source files and their functions remain unchanged.
The archived transform is authoritative and must exactly match the independently
pinned original basis-preflight row. Regenerated and archived transforms must
agree within a64-epsilon float64 representation bound. Both must still pass the
original1e-10 whitening criterion, and their original/current physical basis
maps must agree within componentwise dot-product roundoff bounds. Every original
source/point/coordinate/charge/state/geometry check remains in the replay path.

The original matched-energy evaluator and Cartesian optimizer are unchanged.
New Cartesian trials never use this historical ten-mode transform. No energy,
force, geometry, state, optimizer setting or scientific tolerance is changed.
The four starts,32 finite-difference cells and search budgets remain fixed.

Both new source files passed AST parsing only. Root owns the fresh Slurm
preflight and subsequent execution; no numerical work ran on login.
