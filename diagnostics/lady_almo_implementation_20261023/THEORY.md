# Variational AO-block ALMO implementation

## Scope

`almo.py` implements real closed-shell absolutely localized occupied orbitals
over a full PySCF RHF/RKS backend. The application is the declared fixed-geometry
PBE0 metal + whole-shell effective-singlet partition. The caller supplies the
unchanged molecular Hamiltonian, basis/ECPs, integration grid, fragment AO
membership, and doubly occupied orbital counts. This implements an established
ALMO variational constraint with a new numerical optimizer; it does not establish
the physical validity of a lanthanide effective-singlet state or the completeness
of an energy-decomposition interpretation.

The occupied-support restriction follows the published ALMO/SCF-MI construction:
orbitals are expanded only in their own fragment's atomic basis, while orbitals
belonging to different fragments can overlap. See [Khaliullin, Head-Gordon and
Bell, JCP 124, 204105 (2006)](https://pubmed.ncbi.nlm.nih.gov/16774317/) and
[Khaliullin et al., JPCA 111, 8753–8765 (2007)](https://pubs.acs.org/doi/10.1021/jp073685z).
The optimizer here is not a code translation of either paper's SCF algorithm.

## Density and derivative

Let `S` be the full AO overlap and `C` the AO-by-occupied real coefficient
matrix. Columns are grouped by fragment; each column is zero outside its own
fragment's AO rows. Require `M = C.T @ S @ C` to be positive definite. Then

```text
D = 2 C M^-1 C.T
Tr(D S) = 2 nocc
D S D = 2 D.
```

The inverse occupied-overlap matrix incorporates interfragment
antisymmetrization. `D` itself generally has off-diagonal fragment blocks. A
block-diagonal density or a zeroed interfragment Fock is a different problem.

For the full backend energy `E[D]` with `dE = Tr(F dD)`, define `A = C M^-1`.
Direct differentiation gives

```text
dE/dC = 4 (F A - S A C.T F A).
```

Only the entries permitted by the coefficient support are variations of the
ALMO problem. `ALMOProblem.evaluate` returns the full derivative to make
independent checks possible; the optimizer selects the allowed blocks.
`M^-1` is applied with Cholesky solves. The code checks finite input, exact
fragment partition, total electrons, support, occupied rank/conditioning,
`Tr(D S)`, and the actual density idempotency residual on every evaluation.
It rejects near-singular overlap without silently deleting basis functions.

## Constraint-preserving parameterization and search

For each fragment, construct a square whitening matrix `X_f` such that
`X_f.T S_ff X_f = I`, and write `C_f = X_f U_f`, `U_f.T U_f = I`.
The within-fragment orthonormality is a choice of occupied frame, not an extra
physical restriction: any full-rank supported occupied block spans the same
subspace after this normalization. Each fragment therefore contributes a
Grassmann manifold, the set of its occupied subspaces.

The whitened gradient is `G_f = X_f.T (dE/dC)_f`. Its horizontal projection is
`g_f = G_f - U_f (U_f.T G_f)`. The stopping gradient is the combined Frobenius
norm of these projected blocks in Hartree per dimensionless orbital variation.
The polar retraction of `U_f + alpha eta_f` keeps the next frame orthonormal
without mixing AO support between fragments.

The implementation uses limited-memory BFGS with tangent-vector projection
transport, a positive local orbital-gap preconditioner, and a bounded Armijo
line search. Absolute local Fock energy differences with a fixed positive floor
define only the optimizer preconditioner; they do not change the energy or its
derivative. Invalid occupied-overlap trial points are rejected. A failed BFGS
line search permits one documented reset to preconditioned gradient descent.
All trial energies, rejections, steps, restart attempts and accepted states are
recorded. Floating-point energy slack is explicitly recorded per trial.

Success requires both the specified projected-gradient norm and accepted energy
change criteria. No energy-only convergence is allowed. A stationary initial
guess has zero preceding change and may qualify immediately. A finite budget,
backend failure, invariant failure or failed search yields a nonconverged result
containing the last accepted density/orbitals and history. Invalid initial inputs
or an initial evaluation failure raise before there is an accepted energy;
the calling run harness must archive that exception. A converged result certifies
first-order stationarity, not a minimum, uniqueness, or a physical spin state.

## Backend contract

The full potential is `mf.get_veff(mol, D, hermi=1)`, the Fock is `hcore + veff`,
and the energy is `mf.energy_tot(dm=D, h1e=hcore, vhf=veff)`. PySCF's tagged
potential is passed intact to its energy routine, preserving its Coulomb and
exchange-correlation bookkeeping. There is no Fock masking, DIIS, density
mixing, solvent addition, or implicit ordinary SCF in an ALMO evaluation.
Backend implementation references: [PySCF RKS energy/potential source](https://pyscf.org/_modules/pyscf/dft/rks.html)
and [PySCF restricted SCF source](https://pyscf.org/_modules/pyscf/scf/hf.html).

An existing numerical grid is retained. If none exists, the configured grid is
built once without density-dependent point deletion, then retained throughout
the solve. Matching ordinary/fragment comparison calculations must explicitly
preserve their declared grid conventions. The public backend should remain
unchanged while the problem is active. Complex, unrestricted and solvent-wrapped
backends are rejected; custom energy/Fock wrappers require separate validation.
The current implementation evaluates the full potential each time rather than
accumulating differences between successive trial densities.

`release_constraints` explicitly runs the ordinary SCF kernel from the final
ALMO density on the supplied backend. It changes only the requested convergence
settings and SCF state; it does not add a fallback. The caller should use an
independent same-Hamiltonian backend when preserving prior SCF state matters.

## Interpretation and admission

The independent test harness owns numerical validation and a separate published
benchmark. Required tests include independent density reconstruction, orbital
finite differences, occupied-block gauge invariance, rank/support rejection,
one-fragment recovery of ordinary SCF, constraint release, and the expected
variational energy ordering within numerical tolerances. All scientific
calculations run under root-owned Slurm `standard_expand` submissions.

For identical Hamiltonians and a sufficiently minimized unrestricted solution,
the variational relation is `E_full <= E_ALMO <= E_frozen`. Local SCF roots can
complicate the first comparison, so a violation must be reported and investigated,
not corrected by relabeling roots or shifting energies. Frozen-to-ALMO lowering
describes polarization within the specified fragment AO spaces. Further lowering
upon release describes the relaxation enabled by that specified support
constraint, with the chosen finite basis and fragment electron counts. It does
not alone supply a thermodynamic La/Dy exchange free energy.
