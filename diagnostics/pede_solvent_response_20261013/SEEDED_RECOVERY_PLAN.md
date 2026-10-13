# Conditional calcium vacuum recovery — 13 October 2026

## Question and current decision

Can an explicitly read, converged same-metal water wavefunction initialize the
unchanged sample1 calcium vacuum Hamiltonian if its primary atomic-guess attempt
fails? This tests numerical initialization; it does not change chemistry to
improve a class decision.

At16:14UTC job3372 was unconverged after50 DIIS iterations and8 TRAH macro rows.
The TRAH gradient increased167→224, but every macro step was accepted, energy
decreased substantially, trust radius grew and negative-gap warnings diminished.
An independent read-only review supports continuing the original finite attempt.
The health alert is evidence for review, not proof of stagnation or permission
to accept an unfinished energy. No limit extension or concurrent duplicate.

## Fixed recovery, prepared only

- One fresh sample1 Ca vacuum endpoint. Preserve the failed/unfinished original
  and its result in the primary eight-cell denominator.
- Copy the actual converged sample1 Ca water `endpoint.runtime.gbw`; SHA256
  `1d718711e1aebe6720d8b4e5bbf0fab90cba9a930cc2aeb096f843a1c71421a3`.
  Pin its accepted source collection and source state as well as the seed bytes.
- Preserve the primary vacuum XYZ byte for byte:192 atoms,14 caps, no explicit
  waters, charge−3, physical singlet,838 electrons. No solvent-state, donor,
  cofactor, proton, element or geometry changes.
- Preserve native ORCA6.1.1 r2SCAN-3c, DefGrid3, TightSCF, basis and composite
  corrections. The target remains vacuum. Explicitly retain MaxIter125, the
  existing effective native bound. No additional gCP/dispersion/xTB parameters.
- Change initialization only: explicit `Guess MORead` and `MOInp "initial.gbw"`
  with NoAutostart. Use the already executed native MORead staging pattern.
  Do not force another spin, solver tolerance or method to obtain convergence.

## Admission and reuse

Use a separately named minimal preflight/adapter; do not bypass or relax the
primary frozen-input checks. Require the actual runtime MORead/read-file lines
outside echoed input, reject a model-potential fallback, and verify the printed
vacuum Hamiltonian, state, native basis, normal termination and convergence.
Retain final solver residuals, populations, wavefunction, execution receipt and
all failures. The seed initializes orbitals; it is not the target solvent model.

Parser checks may use archived real MORead and primary outputs, with corrupted
copies explicitly restricted to negative parser tests. They cannot demonstrate
that the new Ca recovery converges. All tests and molecular work use Slurm.

Only an original terminal failure or a separately documented reviewed stop can
activate this numerical recovery. If the primary converges, do not launch this
contingency merely because it was prepared. Root owns that decision and all
submissions. Original collector3377 remains responsible for the primary matrix;
any later recovery comparison must remain separately named and reuse the seven
successful original cells exactly. One converged seed does not establish the
electronic ground state, a solvent-model fix or improved biological accuracy.

## Resources and ownership

Preparation has no molecular cost. A needed recovery is one endpoint on12 MPI
ranks under current shared Slurm QoS, with allocation-derived memory and the
existing native binding fixes. No memory/time override. Small parser checks use
one CPU. Record actual costs and arm completion/health wakes before yielding.

Agent `la_parameter_scope` owns preparation in
`reference_preparation/seeded_recovery/` and matching
`workspaces/pede_solvent_response_20261013/reference_sample1_Ca_seeded_v1/`.
Root owns the scientific decision, execution, comparison and current checkpoint.
