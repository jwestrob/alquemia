# Fixed-charge solvent forces: finite interface test

Declared 16 October 2026 before new native solvent calls. This is a small
qualification of an existing physical component, not a new classifier.
Root owns implementation and Slurm execution. Production remains unchanged.

## Purpose and energy

An affordable protein mechanics model needs solvent forces without whole-protein
SCC. The already-used xTB 6.7.1 `TBorn` implementation provides `addGradient`
alongside `getEnergyParts`. Bind both to the same fixed physical charge vector:

`G(q,R) = Gelec(q,R) + Ghb(q,R) + Gsasa(R) + Gshift`.

Here q is held fixed as coordinates R change. This includes cavity, Born-radius,
hydrogen-bond and ALPB shape derivatives; it is not a charge-response derivative.
No SCC, MACE, quantum calculation, optimization, fitted charge or new solvent
parameters occur. This is distinct from the current self-consistent
GFN2(ALPB) minus GFN2(vacuum) composite. It does not qualify a future full model.

## Exact inputs and implementation

Reuse all four intact LBT3 source/metal geometries and their actual vacuum and
water population artifacts from accepted job12510, each256 atoms and charge−2.
Preserve the archived eight-decimal populations without renormalization. La and
Dy keep their real atomic identities. These are consumed real fixtures.

Reuse the built xTB6.7.1 static library, module files and mctc-lib from the
successful PedE native replay v3. Compile only a separate thin wrapper with
the original pinned compiler and libraries, inside a fresh workspace on Slurm.
Never modify the old build or environment. Pin source, modules, static/shared
libraries, compiler, wrapper, binary and inputs.

Use the unchanged native GFN2-water factory: ALPB/P16, GBOBC, gsolv,
298.15K,230 surface points, zero salt. Input XYZ is angstrom; the native
coordinates are Bohr, gradient Hartree/Bohr. Initialize the inout gradient to
zero before `addGradient`. The exact source ignores the zero-length qsh argument.
Report force as minus the gradient with explicit units.

## Frozen finite test and gates

1. Four water-charge calls must reproduce the four archived native component
   energies. Use the existing eight-decimal charge-rounding bounds plus1e−9Eh;
   bounds above1e−6Eh fail as uninformative. Check charge closure within1e−5e.
2. Four vacuum-charge calls establish new fixed-charge energy/analytic gradients.
   No comparison with the old self-consistent energy is called an identity test.
3. On the7CCO La/Dy fixtures, test two real-atom directions: selected-metal x and
   Asn5 OD1 x. Centered coordinate displacements are±1e−4 and±5e−5angstrom.
   These16 native calls are derivative probes of real fixtures, not new observed
   structures or a numerical quantum gradient. For each direction/step require
   `abs(FD - analytic) <= 1e−5 Eh/angstrom + 1e−4*abs(analytic)`.
   Report residuals and both step sizes, not just acceptance.
4. Repeat the unchanged7CCO/La vacuum-charge call and jointly translate all its
   atoms by(1.25,−0.75,0.5)angstrom. Require energy error≤1e−10Eh and maximum
   Cartesian gradient error≤1e−9Eh/Bohr. These tolerances check numerical
   reproducibility/invariance, not biological accuracy.

Total26 component calls plus one isolated wrapper build. No automatic tolerance,
grid or charge changes after a failure. Retain all results and stop the dependent
force use if a required gate fails. A rotation/grid study and a full-protein
charge convention remain separate prerequisites for a full-model release.

This stage does not project capped QM populations onto Mex, choose whole-protein
spin states, add Coulomb interactions or combine incompatible Hamiltonians.
The separate local-field scout and whole-protein interaction ledger must resolve
those questions before a new preference comparison.
