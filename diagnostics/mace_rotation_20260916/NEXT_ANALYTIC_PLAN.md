# Proposed analytic-dipole correction pilot — not yet approved

## Question and engineering choice

Can analytic monopole/dipole evaluation remove the demonstrated fixed-axis
finite-displacement error while keeping the memory and runtime advantages?
Use the same medium checkpoint and consumed structures. Implement an opt-in
research version; keep current finite-displacement outputs and baseline intact.

Evaluate analytic derivatives of the current radial kernel
`erf(r/(2w))/(r + 1e-6 Angstrom)` for scalar and dipolar source/projection terms.
Retain widths, normalization, exclusions, self terms, global charge constraint,
float64 weights and boundary convention. The small denominator regularizer
remains explicit: this is not a claim of an exact unregularized Gaussian
Coulomb functional. Use bounded blocks and recomputation for first derivatives.
Do not build a whole-protein Hessian or use numerical DFT derivatives.

Removing finite displacement changes the model's field inputs. It therefore
requires a new implementation/method version and validation of learned-density
response; it is not automatically a validated accuracy improvement.

## Tests and run inventory to agree

1. On the same four archived core densities, verify analytic tensor contractions
   and first derivatives against an independent small-core autograd reference
   of the same radial kernel. Check self/excluded terms, charge closure,
   rotational covariance and energy/force signs. Kernel/autograd comparisons
   use absolute 1e-8 plus relative 1e-10, as in the prior kernel checks.
   No fitted offset or new label.
2. Eight complete core energy/force calls: four unrotated and the same four
   rotated by the existing 37-degree transform. Compare changed energies and
   fields with the preserved finite-displacement results; do not require those
   approximations to be identical.
3. If core numerical checks pass, four full-protein calls: La/Ca primary and
   the same rotated pair. Preserve all 9,141 atoms, states and water inventory.
   Reuse existing DFT endpoints to recompute the hybrid partition residual.

Total proposed new model calls: **12**; new DFT: **0**. Same A5000 allocation
share (16 CPUs, 64,474 MiB host RAM); roughly minutes based on the old backend,
with new-kernel timing measured rather than assumed. No project compute/time
budget. Keep the cluster QOS limit and all receipts. No relaxation, entropy,
new biological case, large checkpoint, classifier fit or automatic promotion.

Existing physical thresholds stay fixed: 0.01 kcal/mol endpoint/contrast
rotation error, 0.001 eV/Angstrom force error, 1e-5 e charge closure and
2 kcal/mol partition residual. Report the last separately; solving rotation
is not expected by itself to solve partition sensitivity. A failed physical
check is a result, not permission to adjust a threshold or model parameter.

Deliver actual changed energy components, learned-density shifts, derivative
checks, old/new costs and all unsuccessful outcomes. Only after this numerical
repair should a larger-checkpoint capability test or broader accuracy benchmark
be proposed. Current evidence does not establish broad La/Ca discrimination.
