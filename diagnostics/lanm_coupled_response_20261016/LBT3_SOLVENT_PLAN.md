# Whole-peptide solvent and response scout

Declared16October2026 after the fixed eight-cell MACE result and real La/Dy
parameter qualification12488, before new whole-peptide calls. Root owns execution.

## Why this experiment

The whole-peptide vacuum model selects the La-conditioned7CCO structure for La
and Lu-conditioned7CCN for Dy, whereas the matched pocket selects7CCO for both.
However, fixed-source La/Dy contrasts still differ by72.05kcal/mol, almost as
much as the pocket's73.69. Scaffold representation affects basin selection but
does not itself remove starting-structure dependence. This is a mechanistic
lead, not validation of preference or an experimental population.

Before optimization, determine whether that selection survives a consistently
included solvent term and whether the corresponding forces can support a small
covalent-preserving accommodation. The chemically intact256-atom peptide is small
enough to evaluate directly; no additional classical protein energy is added on
top of MACE. The native peptide mechanics supplies topology and geometry checks,
not duplicated bonded/nonbonded energy.

## Fixed state, endpoints and model

Reuse the exact whole-peptide preparation12469, both7CCO/7CCN model1 sources,
both metals, charge−2. All atom identities, proton inventory, terminal assumption,
dry explicit-water inventory and initial coordinates stay fixed. The NMR water
absence and unresolved experimental terminal chemistry remain limitations.

For each source/metal define

E(q) = E_nativeOMOL,vac(q) + E_authorGFN2,ALPBwater(q) − E_authorGFN2,vac(q).

The force is the negative derivative of these same three terms, including the
GFN2 self-consistent charge response. This is a separately named exploratory
Hamiltonian, not the PQQ frozen-charge solvent policy. MACE uses physical La1/Dy6
multiplicities; originalGFN2 has spin-independent f-in-core/uhf0 bookkeeping.
That hybrid state approximation is explicit and not an open-shell validation.
No solvent term or old native-Dy component is silently reused.

Run exactly eight standaloneGFN2 analytic energy/gradient endpoints:
two sources × La/Dy × vacuum/ALPBwater. Reuse the already-qualified installed
xTB6.7.1 executable, parameters, state/temperature/accuracy conventions. No
geometry optimization, SCF restart, blanket retries or new environment. Retain
all missing/failed cells. These small local-peptide vacuum calls do not restart
the failed whole-LanM vacuum-subtraction experiment.

Obtain four nativeMACE analytic forces on the exact whole origins in a separate
warm GPU job, with the same pinned model/head/state/runtime as12471. Recomputed
energies must reproduce saved fixed evaluations within1e−7eV; these four calls
are counted as new derivative evaluations. No full-protein MACE potential or
Dy f-electron accuracy is inferred merely from executable support.

## Numerical admission and comparison

Keep the existing standalone endpoint gates: actual GFN2 method/version,
electron count/spin/net charge, accuracy0.0001,500maximumSCCiterations,300K
electronic smearing, normal convergence/termination, charge closure≤1e−5e and
JSON/printed energy consistency≤1e−7Eh. ALPB water uses gsolv,298.15K solvent,
normal230-point surface and no ionic screening. Analytic engrad atom order,
coordinates≤1e−6angstrom and energy≤1e−8Eh must match the same evaluated input;
reject a reported numerical gradient. Record full gradients in Eh/bohr and MACE
forces in eV/angstrom; convert once when assembling the composite force.

Report every source separately and the common two-source matrix. Row minima
are mathematical selection, not population. Report the solvent contribution to
the paired/source contrasts, any change in selected geometry and availability.
Do not compare with unsolved Mex solvent rows, sum EF-site energies, calculate
Kd, fit a label threshold or silently replace the existing vacuum result.

In parallel prepare a geometry-only catalogue from the actual whole-peptide
graph: peptide phi/psi, supported heavy-atom sidechain bridge torsions and metal
translations. Preserve rings, conjugated groups, covalent bonds and attached H;
no synthetic caps or disconnected fragments move. Reuse existing kinematics,
physical-displacement normalization and rcond1e−10 rank diagnostics. The same
catalogue serves both metals, selected later by individual and differential
forces. No score fitting or Asn-only prescription from a different protein.

Before using a new composite gradient to optimize, test a small shared physical
direction with central finite differences at0.001 and0.0005angstrom maximum
displacement, under a separately pinned real-source direction. Fixed acceptance
is residual≤0.2kcal/mol/angstrom +5%of the derivative magnitude, with both step
results and their difference reported. This is a later finite interface test;
the present twelve molecular calls do not establish it. No numerical DFT
gradients, full Hessian, long trajectory or unbounded optimization.

## Execution

Eight CPU endpoints use four independent one-thread workers on standardSlurm.
Four MACE force calls use one warm GPU with current granted resources. Derive
memory from allocation; no memory/time override. Geometry-only preparation uses
oneCPU. Pin inputs, executables, helper code and actual controls, count all calls,
failed attempts and allocated CPU/GPU time, and arm completion wakes. Keep
PQQ production, canonical midpoint and all old LanM failures unchanged.
