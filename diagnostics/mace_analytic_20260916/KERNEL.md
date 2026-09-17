# Analytic monopole/dipole kernel

This opt-in model change is the zero-displacement derivative expression of the
previous realspace representation. It retains its radial regularizer and saved
normalization/self-interaction coefficients. No model weights are changed.

For receiver i, source j, d = r_i-r_j, r = |d|, width w, and epsilon = 1e-6 A:

```
f(r) = erf(r/(2w))/(r+epsilon)
A = f'(r)/r
B = (f''(r)-A)/r^2
phi_i = C sum_j [q_j f - A (p_j . d)]
grad(phi_i) = C sum_j [q_j A d - A p_j - B (p_j . d) d]
E_pairs = 1/2 sum_i [q_i phi_i + p_i . grad(phi_i)]
```

C is the unchanged backend FIELD_CONSTANT/(4*pi). Self pairs and cross-batch
pairs are excluded. Coincident distinct physical centers fail explicitly;
excluded distances are replaced only inside the masked arithmetic. Existing
analytic GTO self terms are then added under their original enable flags.

The scalar field projection uses the saved l0_factors. Vector projections use
l1_factors * old_offset times the analytic gradient, removing the 1/offset
already included in the saved finite-difference factor. The backend's spherical
[y,z,x] convention is preserved. The old offsets no longer displace coordinates;
using their product with the saved coefficient preserves normalization.

`scripts/mace_analytic.py` sums 256-by-256 physical-atom blocks. Its custom first
backward rebuilds one analytic block and uses autograd on that expression.
This supplies coordinate, charge and dipole derivatives through the learned
field-response layers. It is analytic differentiation, not a finite-difference
force or numerical DFT gradient. Local neural edge/node checkpointing remains.
Only linear-size inputs are saved between blocks; pair work remains quadratic.
Training, Hessians, higher derivatives and unsupported architectures are rejected.

The independent test reference differentiates the scalar radial function in
Cartesian coordinates with torch.func.jacrev and vmap on the four small real
cores. It does not reuse the optimized f'/f'' formulas. Field values and the
first derivatives of a field-dependent scalar objective are compared, together
with batch exclusion, rotation of energies/forces/features and bitwise weight
identity. Tests deliberately cross block boundaries at sizes 17 and 64.

The epsilon regularizer prevents this from being exactly the unregularized
Gaussian Coulomb functional. It is retained to isolate removal of finite
stencils. Self terms are inherited, not silently rederived under a new convention.
Actual model outputs and changed energies must be measured; improved numerical
symmetry alone does not validate La/Ca discrimination or the hybrid partition.
