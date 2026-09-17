# Next numerical step: coupled CPU polarization on real frameworks

Declared after the parameter-only result and before any Tinker energy output.
Active-goal autonomy applies. This is a solver implementation control toward
the MACE hybrid, not a biological comparison or a metal-affinity score.

## Question and input selection

Does the pinned CPU solver converge on the three existing real GGR1GLG,
alpha1F6S and alpha6IP9 protein/water frameworks, honor the restored source
freeze through its actual GK-coupled response, and run at a practical cost?

Reuse exactly `tinker_capability_v1` physical/framework mappings, source
coordinates, native AMOEBA2018 atom types and full parameter file. Keep source
hydrogens, disulfides, waters (0/2/3) and protonation unchanged. The selected
metal is explicitly not in these framework-only controls; no result may be
reported as a whole-system La/Ca energy or evidence that metal omission is
acceptable in a scorer. No generated molecular fixtures or biological labels.

## One model, twelve requested energy evaluations

For each of the three frameworks evaluate four predetermined states:

1. All framework sites polarizable, mutual convergence target1e−5.
2. The already exported QM-source support frozen with the post-initialization
   adapter, target1e−5. GGR uses its recorded extended-core support.
3. The same frozen support, target1e−7.
4. State3 after one recorded proper rotation and translation. Use the existing
   project rigid-transformation utility and preserve the exact transform in
   the manifest before execution; no geometry perturbation or optimization.

Use the native CPU GK-coupled solver, no periodic boundary, no spatial cutoff
(native isolated-system1e12Å sentinel), no added ionic screening. Set internal
dielectric1 explicitly; native GK source hardcodes solvent78.3 in both field
and energy routines, which will be verified/pinned instead of using an ignored
input keyword. GK uses the native AMOEBA2018 SOLUTE
radius inventory and pinned defaults (GKC2.455, descreen offset0.30Å, neck/tanh
and hydrogen-descreening flags from the native initializer). Export actual
per-atom radii/scales and global flags before evaluating; fail on unexpected
missing parameters. This is distinct from the earlier OpenMM Bondi-based
inventory and cannot inherit its solvent energies.

Activate native permanent multipoles, polarization and GK solvation only.
Disable bonded, vdW, charge-transfer and other potential terms through native
keywords/settings. Preserve and identify any nonpolar contribution included
inside native `es`; do not rename it pure reaction-field energy. Freeze the
explicit remaining solvent settings in the prepared manifest before execution.
They are engineering choices, not fitted to the GGR/alpha ordering.

Keep permanent charges/multipoles and damping identical across states1–3.
Only the allowed induced variables and stated SCF tolerance change. Apply the
source mask after `mechanic`, retain nonzero damping parameters, and record
both vacuum and solvent induced dipoles. The inspected `induce0c` zeroes these
arrays at the start; verify the resulting frozen values rather than assuming
the parser-mask test guarantees them. Do not use the GPU code for this model.

## Energy accounting and acceptance, frozen before outputs

The native `energy()` function returns the activated sum. Record unrounded
`em`, `ep`, `es`, total, all other component values, solver settings and status.
Tinker energy units are kcal/mol. Verify total equals `em+ep+es` within1e−8
kcal/mol and every disabled component is zero. Record repeated internal SCF
work if the native routines invoke induction more than once per requested
energy evaluation; execution cost includes all such work.

Require successful native convergence, finite energies/radii/dipoles, unchanged
state inventory and frozen source dipoles ≤1e−12 eÅ for both vacuum and solvent
solutions. Require the state2–3 energy change and state3–4 rigid-transform
change each ≤0.01kcal/mol, the existing numerical tolerance. This is far below
the2kcal/mol partition scale used by the hybrid; no biological classification
is used to choose the tolerance. Report the individual components as well as
the total. No monotonicity gate is assumed for AMOEBA's distinct d/p scaling
schemes. No partition-invariance or La/Ca accuracy claim from these controls.

Capture elapsed time, actual/allocated CPU use, memory, per-call timing and
all failed attempts. Use the existing64CPU standard/memory allocation pattern,
one node, existing excluded nodes, explicit finite task manifest, no project
time/CPU budget. No DFT, MACE, charge fitting, gradients, trajectories or
optimization. No further energies outside this12-call inventory without a
new recorded experiment definition; authorization remains autonomous.

## What follows a pass

If affordable and numerically consistent, implement the complete source-charge/
multipole/covalent-boundary model and the full-cavity subtraction in
AMOEBA_ACCOUNTING.md, with a separately declared electric-field quality check
and native La/Ca state accounting. Reuse the existing QM densities and MACE
short terms where compatible. A failed solver check instead identifies the
specific technical/physical defect to repair. This control cannot qualify a
new discriminator, generate an aquo reference or justify production promotion.
