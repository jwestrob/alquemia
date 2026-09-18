# Frozen-response solvent transfer: qualify the energy before scoring

## Question and scope

Can a matched change of reaction-field model include the permanent protein
multipoles and both AMOEBA induced-dipole sets without dropping or double
counting polarization terms? This is an implementation and numerical test on
the already consumed GGR 2FW0/2FVY Ca/La states. It does not change a biological
label, geometry, source charge fit, water inventory, or production protocol.

Stage A executes four native operator replays, one per archived endpoint.
Each replay uses the existing GK cavity, multipoles, supplied driving fields,
and converged induced dipoles. It performs one Born-radius evaluation, one
permanent-field evaluation, one mutual-field evaluation, three static GK
electrostatic energy evaluations, and one nonpolar component evaluation.
It performs no iterative response solve, quantum calculation, MACE inference,
force calculation or continuum solve. One 64-CPU allocation; measured costs and
all failures retained. Existing native kernels and library remain unchanged.

The four inputs come from `trial_gk_expansion_native_v1/collection_job_1201162.json`.
Only primary geometry, Ca/La, tight supplied-density response states are used.
Standard-field controls are not substituted. This first stage is useful even
if the conductor discretization fails its separate numerical qualification.

## Energy derivation

Use the native Tinker conventions, with electric fields in e/angstrom^2,
dipoles in e*angstrom and the recorded conversion constant c=electric/dielec.
Let A be the common mutual-response operator and F_d,F_p the two actual
driving fields. The stationary bilinear functional is

    H(mu_d,mu_p) = c/2 [mu_p^T A mu_d - mu_p^T F_d - mu_d^T F_p].

At A mu_d=F_d and A mu_p=F_p it reduces to the native induction energy
`I=-c/2 mu_d^T F_p`. This is a stationary expression, not a claim of a minimum
over two independent dipole vectors. Verify both equations using native
`ufield0d`, not a newly implemented interaction tensor.

Let P contain the projected quantum-source monopoles plus the actual permanent
environment multipoles. For a symmetric reaction operator R, the part of the
total stationary functional that depends on R is

    C_R = 1/2 (P+mu_p)^T R (P+mu_d)
        = E_R(P+(mu_d+mu_p)/2) - E_R((mu_d-mu_p)/2),
    E_R(v) = 1/2 v^T R v.

The second equality is an algebraic way to evaluate the cross term. It does
not replace AMOEBA's two response states with one averaged physical state.
It includes permanent/permanent, permanent/induced, and induced/induced
reaction terms. Native GK fields independently give

    C_GK = E_GK(P) - c/2 [(mu_d+mu_p)^T F_RF(P)
                                  + mu_p^T F_RF(mu_d)].

Stage A checks these two evaluations against each other, and checks H against
the archived induction contraction. The nonpolar energy is recorded separately
to compare E_GK(P) with the archived static solvation component.

A later frozen-response transfer would be `Delta=C_new-C_GK`, evaluated with
the same P, mu_d and mu_p in both models. Add `Delta_Ca-Delta_La` to the old
hybrid raw contrast. Unchanged intrinsic QM, direct-density interaction,
vacuum mutual response, and MACE short-context terms cancel within this
transfer. Environment-only terms cancel between endpoints only after their
physical identity is demonstrated. No aquo reference or calibration is implied.
This is a perturbative evaluation at GK response, not self-consistent response
to the new boundary. No gradient or relaxation term follows from it.

For discretized ddCOSMO, reciprocity is approximate. Its quadratic scalar
defines the symmetric part of the numerical reaction operator; the above
energy-difference form must retain the existing reciprocity/convergence gates.
No exact variational property of the unsymmetrized discretization is asserted.

## Frozen Stage A acceptance

- Source coordinates, atom order, charges, radii and global multipoles must
  reproduce their archived inputs; maximum numeric discrepancy 1e-12 in native
  units, with original files and mappings pinned.
- Born radii agree within 1e-10 angstrom. Source induced dipoles remain zero.
- Native static GK energy agrees with archived solvation minus native nonpolar
  energy within 1e-7 kcal/mol.
- Direct d/p reaction fields agree within 1e-10 e/angstrom^2.
- The two C_GK expressions agree within 1e-7 kcal/mol.
- The bilinear stationary functional reproduces the archived induction
  contraction within 1e-5 kcal/mol. This is tighter than the 0.1 kcal continuum
  discretization screen and allows the declared finite 1e-7 Debye induction
  convergence tolerance; retain actual residuals as well as this energy test.
- Independently computed response residual RMS, after multiplication by each
  site's polarizability and conversion to Debye, must be <=1e-7 Debye. Only
  responsive atoms contribute, using the native normalization over n atoms.
- Mutual operator reciprocity energy discrepancy <=1e-7 kcal/mol.

Missing/nonconverged/incompatible inputs fail explicitly. No new solver-backed
full score is authorized by a passing parser/algebra test. Declare its actual
source representation and convergence inventory before a later transfer pilot.

## Source verification

The installed Tinker 26.2 source is pinned at
`87050685eff8840d312e2a332cc82c33f63c7c3d`. Relevant native routines are
`induce0c`, `dfield0d`, `ufield0d`, `egk`, and `enp`. The existing native field
accounting experiment already verified `I=-c/2 mu_d^T F_p` on real states.
The [Tinker polarization manual](https://tinkerdoc.readthedocs.io/en/latest/text/forcefield/polarize.html)
documents the two fields, common mutual operator and bilinear functional.
Its informal description of the bilinear stationary point as a minimum is not
needed here. The formula above follows directly from the two linear equations.

The [coupled QM/AMOEBA/ddCOSMO paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8444335/)
identifies separate direct/adjoint continuum and induced-dipole variables.
Its fully coupled method does not validate this frozen-response approximation.
Primary article full text was blocked by a browser challenge in this session;
the implemented energy must be justified by the native equations and explicit
Stage A tests, not by an unverified claim of equivalence to that paper.
