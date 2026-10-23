# Paired La/Dy reference chemistry — initial panel

21 October 2026. Root owns submissions. This protocol is part of the activated
explicit-solvent exchange goal, not a continuation of the old whole-peptide
optimizer. All native calculations require Slurm. Local input generation and
saved-data checks use one Python process with all numerical thread controls set
to one, within this agent's four-thread allowance.

## Scientific question

Does the previous frozen-4f electronic approximation preserve La/Dy differences
in water-versus-donor interactions and forces? How does that depend on donor mode,
distance and first-shell coordination? These fixed-configuration results qualify
ingredients, not binding free energies or a protein selectivity prediction.

## Inputs and finite inventory

The exact four supplied structures in Wait et al., JPCB 2023,
DOI 10.1021/acs.jpcb.2c07237, supplement s001 are used: water, acetate bidentate,
acetate monodentate and acetamide. La is replaced by Dy without changing any other
coordinate. Each dimer has central, compressed and stretched configurations;
only the metal moves by -0.15/0/+0.15 Angstrom along the metal-to-donor-centroid
axis (one coordinating oxygen, or both acetate oxygens for bidentate acetate).
The supplied geometry is not assumed to be a minimum of either new method.

Two explicit aquo controls contain eight and nine water molecules. Eight oxygen
positions form a square antiprism at 2.50 Angstrom. The nine-water control adds
one axial water and redistributes the other eight into two staggered four-member
rings at z/r=+/-0.35; all oxygen radii remain 2.50 Angstrom. Water internal
geometry is 0.9572 Angstrom and 104.52 degrees, with the hydrogen bisector
pointing away from the ion. These deterministic stress fixtures are not measured
structures or equilibrium samples. La and Dy use identical coordinates.

There are 28 physical configurations: 4 dimers x 3 distances x 2 metals,
plus 2 aquo clusters x 2 metals. Two electronic recipes give 56 fixed energy and
gradient cells. The first pilot is exactly eight central cells: water and
bidentate acetate x La/Dy x both recipes. Remaining cells require actual pilot
admission and root resource review; no blind automatic launch or retry.

## Electronic recipes

1. `ecp_pbe0`: PBE0-D4, def2-TZVP for light atoms, the exact archived lcecp-1-TZVP
   metal basis/ECP and AutoAux-J auxiliary. La ECP46 and Dy ECP55 leave eight
   metal electrons for Ln3+; both effective singlets. Dy physical multiplicity6
   is separately recorded and its 4f9 shell is frozen in the core.
2. `ae_pbe0`: PBE0-D4, scalar ZORA, ZORA-def2-TZVP light atoms,
   SARC-ZORA-TZVP on La/Dy and SARC/J. No ECP. La singlet, Dy sextet with explicit
   4f electrons. No spin-orbit coupling; a converged determinant is not a
   ground-state or multiplet qualification.

Both use ORCA6.1.1, analytic gradients, DefGrid3, VeryTightSCF, PModel initial
guess, NoAutostart and MaxIter300. RIJCOSX is used for both. This changes basis,
core and scalar-relativistic representation together: differences cannot be
uniquely assigned to frozen-4f error. Agreement does not validate PBE0 itself.
The next reference tier requires correlated small-system and basis/BSSE tests,
and hydrated ligand competition; it is not silently supplied by these results.

## Resources and admission

Each pilot singleton uses four MPI ranks on one standard node,
`--account=banfield --partition=standard --qos=standard_expand --ntasks=4`.
Default memory only, no memory/time flags. MaxCore is floor(0.75 x granted
memory MiB / MPI ranks), with the actual grant and input recorded. One numerical
thread per rank. ORCA is invoked directly and spawns its own MPI workers.

Admission requires normal native termination, explicit SCF convergence, a full
finite gradient, matching engrad energy/coordinates/elements, expected electron
count and multiplicity, and correct ECP/ZORA evidence. Every printed convergence
residual, spin expectation and warning is retained; native admission is not
physical-state validation. Failed or interrupted results remain immutable. No
automatic SCF fallback, changed multiplicity, optimization or continuation.

## Balanced interpretation

For donor L, the reaction is Ln(H2O) + L -> Ln(L) + H2O. The **Dy-minus-La
exchange contrast** is E(DyL)-E(DyW)-E(LaL)+E(LaW); isolated ligand and water
energies cancel exactly. Counts and net charge balance in the four-cell cycle,
including acetate (complex +2) versus water/amide (+3). A negative contrast means
greater Dy stabilization relative to water at these fixed coordinates, not a
solution affinity. No absolute La-versus-Dy total-energy subtraction is reported.

The aquo coordination contrast similarly compares [E(DyW9)-E(DyW8)] minus
[E(LaW9)-E(LaW8)]; the released/added water cancels. Report full forces and metal
radial projections, matched metal-force differences and radial deformation
energies. The +/-0.15 Angstrom points are physical distortions, not a tight
finite-difference derivative certification.

Uncorrected dimer energies are not directly comparable to the paper's
counterpoise-corrected DLPNO-CCSD(T) interaction energies. The published
monodentate-acetate discrepancy (about49 kcal/mol for AMOEBA22) is preserved as a
known limitation, not repaired by relabeling the geometry or changing the gate.

## Later required chemistry qualification

After the pilot, add matched hydrated ligand replacement configurations from
the actual explicit-solvent model and qualify representative many-body response.
Use counterpoise/basis convergence and correlated references on tractable
closed-shell cases; test explicit-f/SOC/state sensitivity where it affects the
La/Dy contrast. Freeze distinct fitting and validation configurations. This
initial panel does not fulfill the complete Stage2 or authorize protein claims.
