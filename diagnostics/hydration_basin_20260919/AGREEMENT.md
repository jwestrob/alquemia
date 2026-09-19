# Coupled water response and occupancy continuation — approved 2026-09-19

Jacob: “go for all of it!”, responding to the successful radial check and the
proposed coupled translation/rotation and occupancy free-energy continuation.
Standing autonomy applies. Preserve the baseline, inputs and prior experiments.

## Stage A: coupled physical response

Use the 20 nonempty metal endpoints in the completed 12-pattern/24-endpoint
occupancy table. Reuse actual native DFT/CPCM energies and analytic gradients.
Every retained variable H2O has three center-of-mass translations (Å) and three
exponential-map rotations about its center of mass (radian). All protein, metal,
outer waters and water internal geometry remain fixed. Use actual H/O masses,
physical inertia and all cross-water curvature terms. No cap/fragment motions.

Define the anchored potential in Cartesian space, before changing coordinates:
V(q) = E_DFT(0) + E_MACE(x(q)) - E_MACE(x(0))
       + [g_DFT,x(0)-g_MACE,x(0)] dot [x(q)-x(0)].
This includes the rotation-map chain rule and its second derivative. It is an
approximation to local DFT/CPCM response, not self-consistent dielectric response.
Native unmasked OMOL0-100M float64 checkpoint; target r2SCAN-3c/NoAutostart/
CPCM(Water)/DefGrid3/TightSCF EnGrad; unchanged charge, spin, inventory, assembly.

Compute cheap projected Hessians from analytic MACE gradients at ±0.001 and
±0.0005 per physical coordinate; no numerical DFT derivatives. Rotations use
one Å per radian as a numerical scale; physical mass/inertia removes that choice
when comparing eigenmodes. Symmetry residual must be ≤max(0.05,0.001*matrix norm)
and coarse/fine residual ≤max(0.05,0.01*matrix norm), in the declared scaled units.
Keep raw unsymmetrized matrices, symmetric matrices, eigenvalues and all errors.
Never clamp negative eigenvalues or delete soft modes to create stability.

Also search the anchored cheap potential from the actual center, using SLSQP
with analytic gradients, at most 200 iterations, ftol1e-10 eV. Trust region per
water: COM translation≤0.20Å, rotation≤0.35rad. These bounds define a local model,
not a production compute budget. A boundary optimum cannot supply unbounded
harmonic entropy. Record projected residual gradients and minimum curvature;
flag changes of basin, overlap, numerical failure and unbounded local response.

One GPU/16CPUs/64474MiB, existing finite executor, 20 tasks. This is approximately
800 initial finite-difference cheap force calls plus local searches and minimum
curvatures. Measure actual wall time/memory/calls; no project time/compute cap.

## Stage B: native physical checks before using free energies

On the four already consumed pilot centers (1F6S11 and6IP9110, both metals), test
(a) the lowest mass-weighted curvature direction, retaining negative values, and
(b) the local proposed displacement direction. If the latter is absent/too small
or nearly parallel (>0.95 absolute cosine) use the largest-curvature eigenvector.
Normalize each direction to max atomic movement≤0.05Å and rotation≤0.10rad.
Run both signs: 16 native analytic energy/gradient endpoints, plus up to four
native checks of interior, stable proposed minima. Use the same existing runner.
Freeze the exact finite manifest before submission; no selection by class labels.

Local energy/even-curvature/gradient criteria retain the radial pilot scales:
energy error≤max(0.02kcal/mol,25% of |DFT even energy|); even-energy error≤
max(0.005kcal/mol,25% of |DFT even energy|); projected gradient error times path
amplitude≤max(0.02kcal/mol,25% of the actual gradient change times amplitude).
Evaluate native energy at proposed minima and actual physical residual forces;
force consistency is necessary, not proof of a complete basin ensemble.

## Stage C: conditional free energies and occupancy

Only numerically stable, physically checked local models can contribute basin
terms. Use physical masses, moments, identical coordinate measures and coupled
curvature. Include nonzero linear terms and check trust-region/basin extent;
no unvalidated log determinant. Retain documented alternative water orientations,
site exclusions, and missing states. No favorable-state cherry-picking.

Reconcile internal-water vibration and solvent terms with the existing liquid
reference before an absolute occupancy claim. If a term/backend is absent,
report a conditional component or missing status rather than zero. A constrained
model estimate must identify its assumptions and is not full binding free energy.
Further native endpoints for validated arrangements must have a fresh finite
manifest; do not run a production rescore or change the baseline/default.

Write actual results, limitations, costs, executable commands and a vault note.
The old native orientation comparators are separate; collect when completed.
