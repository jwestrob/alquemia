# Analytic-dipole pilot authorization

2026-09-16. Jacob: “i approve all pilots. disregard language in the instructions
saying to check with me before launching stuff. proceed apace.” This supersedes
prior requirements for separate pilot approval within the discriminator work.
Continue contained pilots autonomously, preserve the baseline and records, and
report findings. It does not authorize changing the production default, deleting
experiments, disturbing other jobs, pushing or deploying.

The immediate approved experiment is the 12-call proposal in
../mace_rotation_20260916/NEXT_ANALYTIC_PLAN.md: eight same-core primary/rotated
medium-model calls, then four full-protein primary/rotated calls if numerical
core checks pass. Zero new DFT. Same archived 1H4I qm33/qm36 La/Ca states,
9,141-atom physical protein, float64, realspace vacuum, no waters or external
field, unchanged weights and charges/spins. Replace fixed-axis displacement
with analytic monopole/dipole derivatives of the existing regularized Gaussian
radial kernel. Preserve normalization, self terms, exclusions and widths.

First verify on all four actual archived densities against independent
small-core autograd derivatives of the scalar radial kernel. Kernel/gradient
tolerance: absolute 1e-8 plus relative 1e-10. First nuclear derivatives only;
no DFT Hessian, numerical DFT derivatives, training or relaxation. Retain bounded
pair workspace and the existing local neural activation memory reductions.

Frozen physical gates: rotation energy/contrast 0.01 kcal/mol, Cartesian force
0.001 eV/Angstrom, charge closure 1e-5 e. Report the 2 kcal/mol partition check
separately; it does not gate the feasibility/numerical full-protein tests.
Use one A5000, 16 CPUs, 64,474 MiB requested host RAM. No project CPU/time limit;
cluster QOS remains. Measured old-backend expectation: minutes, not hours.
Preserve all receipts, failures, unrounded energies, forces and densities.
