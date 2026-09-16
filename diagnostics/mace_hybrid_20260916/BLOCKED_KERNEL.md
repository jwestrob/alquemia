# Bounded-memory realspace execution

Approved scope: [MEMORY_AGREEMENT.md](MEMORY_AGREEMENT.md). Implementation:
`scripts/mace_blocked.py`, opt-in via the existing manifested runner. The original
isolated software and reference outputs are retained. No production change.

For displaced-charge receiver i, source j and Gaussian width w, the existing
backend uses `h_w(r)=erf(r/(2w))/(r+epsilon)`, epsilon=1e-6 Angstrom, and
`F_i,w=C sum_j q_j h_w(r_ij)`, C=FIELD_CONSTANT/(4pi). Pairs from different
batch systems, and all displaced sites belonging to the same original atom,
are excluded. The original feature projection, dipole offsets, normalization,
self-interaction coefficients and inclusion flags are retained.

The replacement sums512-by512 site blocks. Its analytic backward recomputes
each block with

```
h'_w(r) = [ exp(-(r/(2w))^2)*(r+epsilon)/(sqrt(pi)*w)
            - erf(r/(2w)) ] / (r+epsilon)^2
dL/dq_j = C sum_i,w (dL/dF_i,w) h_w(r_ij)
dL/dr_i = C sum_j,w (dL/dF_i,w) q_j h'_w(r_ij) (r_i-r_j)/r_ij
```

The corresponding negative coordinate contribution accumulates on source j.
At coincident positions, the norm derivative is zero, matching torch's original
norm convention. Electrostatic energy is `0.5 sum_i q_i F_i`, plus the original
self terms where enabled. Autograd includes both the explicit q factor and its
occurrence inside F, and propagates both charge and coordinate derivatives
through all learned field updates. Forces retain the minus-gradient convention.

Saved custom-function tensors are charges, positions, batch IDs and widths;
pair matrices/edge lists are neither built globally nor retained. Peak pair
workspace scales with tile size squared, stored inputs with atom count. Pair
work remains quadratic. Local neural-network activations are a separate memory
cost; host offload remains an authorized fallback. Float64 remains mandatory.
Training, Hessians and higher derivatives fail explicitly; this implementation
provides energy and first derivatives for the approved inference experiment.

## Verification and provenance

CPU kernel comparisons use the four actual computed densities and original
geometries. Tiles17 and64 exercise block boundaries. Compare field/energy values,
charge derivatives and coordinate derivatives against the original autograd
implementation. A combined batch contains those same four real systems and
checks cross-system exclusions. Saved-tensor hooks check that no quadratic
pair arrays are retained. Numerical tolerances were recorded before results.

The full-model CUDA comparison is a separate four-core execution. The executor
refuses full-protein work until energies, forces and density coefficients agree
with the pinned successful reference collection at the declared tolerances.
New cache keys include the kernel ID, tile size and implementation hashes.
All original physical input bytes are reused; no core result is substituted for
a full-system result and no failed attempt can satisfy this gate.

Current implementation manifest: `workspaces/mace_hybrid_20260916/blocked_v6/manifest.json`.
Exact status, execution IDs and measurements belong in the companion result
record after collection. This design note alone does not assert a successful
full-protein evaluation or validation of the large checkpoint.

## Local neural-network workspace

The pair rewrite alone passed the core comparisons but exposed independent
local-network memory requirements in the full protein. `mace_local_memory.py`
now evaluates the unchanged nonlinear interaction in4096-edge blocks. All
neighbors contribute to message and density sums before node normalization
and nonlinearities. An indices-only accumulation adjoint avoids torch.index_add
retaining the large edge-message source values. Complete block recomputation
avoids TorchScript wrapping the checkpoint early-stop exception.

The original per-atom product module runs in256-atom blocks, also checkpointed.
These modules operate independently on atoms after their neighbor sums, so
batching changes no physical neighborhoods. Core validation deliberately uses
128-edge and17-atom blocks to exercise boundaries. All sizes are recorded in
the manifest and actual result. Weights and buffers remain bitwise unchanged.

The new implementation is intended for compatible pretrained checkpoints with
these interaction/product classes and monopole/dipole electrostatics. Larger
weights have not been evaluated. Unsupported architectures, training, LAMMPS,
periodic evaluation or higher derivative requests fail explicitly.

The per-atom sparse tensor products in field updates and the local-energy
readout use the same node batching. Global fields and global charge constraints
are unchanged; only each independent atom's tensor-product workspace is batched.
