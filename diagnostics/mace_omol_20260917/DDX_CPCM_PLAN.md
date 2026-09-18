# A conductor-like protein-boundary component, chosen on engineering evidence

## Why this physical approximation now

The finite-dielectric PCM stage is the cost/convergence bottleneck. The first
GMRES Ca endpoint spends490.249seconds in that system and4.351seconds in the
subsequent single-layer system. Its energy agrees with the native value to
1.309e-6kcal, narrowly failing the predeclared1e-6 equivalence test; the native
solution's true residual also fails the new1e-9 screen. Both failures stand.
No biological classification has been computed with this solver.

Test the conductor-like approximation as a distinct model. ORCA's C-PCM uses
the conductor solution scaled by`f(epsilon)=(epsilon-1)/epsilon`. The pinned
Psi4/ddX host applies exactly this factor to native ddCOSMO energies; ddX
itself returns the unscaled conductor energy. Thus the expensive first PCM
system can be avoided by an established physical approximation, not by
misreporting a missing term as zero.
[ORCA C-PCM definition](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html),
[ddX host scaling](https://raw.githubusercontent.com/psi4/psi4/master/psi4/driver/procrouting/solvent/ddx.py).

This does not reproduce the complete ORCA baseline: the present sphere radii,
surface treatment and frozen projected source charges differ from its quantum
density and Gaussian cavity. It inherits no reference, threshold or prediction.
Do not call it exact finite-dielectric PCM or full quantum CPCM. Keep ongoing
finite-dielectric jobs, their outputs and the failed numerical screens intact.

## Frozen state, energy and execution inventory

Reuse the already pinned full GGR2FW0/2FVY physical source inputs and the exact
original eight groups/twenty roles (three grids, rigid transforms, zeros and
repeats). All twenty are new conductor solves; no PCM success can satisfy this
model's cache. Same default source charges, actual coordinates/radii, common
metal sphere, explicit atoms, eta.1/shift0,eps78.3,salt0,FMM12/12,direct Coulomb
source potential, matching monopole psi,64threads, matrix-free storage and
maxiter1200. Fresh native Model per state, no automatic retry. Use the original
qualified stock ddX0.9.0 environment, not the kernel-adapter build.

Only the model changes to native`cosmo`, with the documented energy factor:

`U_CPCM = ((78.3-1)/78.3) * U_conductor`

`U_conductor = 0.5 * dot(psi,x_conductor)` in Hartree.

Retain the raw conductor energy and coefficients, factor, scaled energy and
single Hartree-to-kcal conversion. Report Ca-minus-La source-self contributions
and their between-structure differences. No full environment or hybrid energy,
affinity direction, electronic response, DFT, MACE, gradient or relaxation.

One64CPU/128GiB allocation, all eight groups serial. Measure setup and source
costs as well as solves, CPU/memory and allocation. No project time/compute
budget; numerical iteration limits remain explicit and failures visible.

## Acceptance and interpretation fixed before conductor outputs

Keep original numerical tolerances: source potential1e-8au, psi1e-12,
contraction/zero/repeat1e-8kcal, reciprocity0.05kcal, primary-to-refined0.1kcal,
rigid0.05kcal, including endpoint, Ca-minus-La and between-structure changes.
Scale energy/cross-contraction checks consistently. Passive energies remain
nonpositive. All physical-state and paired-coordinate checks are unchanged.
The coarse-to-primary change remains separately reported; a coarse failed
reciprocity check cannot disappear if finer grids improve.

Compare against actually available PCM values as a physical-model difference,
not an implementation-equivalence test or a target to tune. Unavailable PCM
values stay missing. This experiment can establish computational feasibility
and numerical behavior of the conductor component; predictive value requires
separate coherent full-energy accounting and benchmark work.
