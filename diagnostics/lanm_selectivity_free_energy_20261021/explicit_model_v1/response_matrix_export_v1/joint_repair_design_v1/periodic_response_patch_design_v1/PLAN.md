# Periodic h/H implementation patch plan

**Source/design only; no engine change or execution release.** Extend an isolated copy of the qualified vacuum CPU refactor (`source_manifest 9ccee436…`), retaining the exact response functional and native AMOEBA terms. The running fit keeps its original binary. Exact file paths, hashes and unresolved dependencies are in `DEPENDENCIES.json`.

## 1. Initial scope and the established tail fact

First target: one metal in the complete neutral, orthorhombic BIO water/Na/Cl system, native CPU MUTUAL/Thole PME with conducting boundary. Preserve all BIO intramolecular terms, salt parameters, polarization-group scaling, intrinsic alpha and fixed POLPAIR. No donor-response, GPU, replica/polymer, solvent-continuum, predictor, alternative polarization solver or nonorthogonal-cell support is added by this plan.

For each primary water oxygen retain the **27 labelled images** `n ∈ {−1,0,1}³`, with `r_I = R_O − R_M + L n`. This is not minimum image. The saved geometry-only evaluation bounds the omitted radial sums in both exact prepared boxes by T0 ≤ 1.1115596485…×10⁻⁸⁰ and T1 ≤ 2.3826936607…×10⁻⁷⁸ Å, using inherited D = 43.492549421087 Å. Exact outward endpoints are authoritative; the displayed values are abbreviated. This is a finite saved-box fact, not a uniform trajectory certificate or a selected physical cutoff.

The implementation computes a conservative finite-image functional and its exact analytic derivatives for a fixed retained image set. It approximates the declared infinite-image model. A numerical error budget for that approximation, retained-sum roundoff, kernel evaluation and incomplete induction remains to be declared before execution. The saved T0/T1 alone do not admit h/H errors or self-consistent total forces.

## 2. Concrete native changes

Line references below are to the pinned refactor source.

| File / routine | Narrow change |
| --- | --- |
| `polar.f`, response storage | Replace 28-entry work arrays and one-water indexing with dynamic image-pair storage. Store primary O index, integer image label, r, f, grad f, phi, phi′/r, other-image sum and **per-image** beta derivative. Keep h(3,n) and beta(n) on primary atoms. No image dipole variables and no global coordinate Jacobian. |
| `induce.f:7923`, `localresponsesetup` | Separate immutable parameter/topology setup from per-coordinate geometry. Add a separately guarded periodic branch with translated roles; leave the qualified vacuum branch and absent/all-zero bypass intact. Replace the n≤28/all-atom-water restrictions only in the new branch. |
| `induce.f:8109`, `localresponsegeometry` | Enumerate all 27 images of every eligible O. Compute the same regularized radial/contact functions and positive prefix/suffix environment sums over **image entries**. Add each image field into its primary O rather than assigning/overwriting it. Evaluate H once at the primary metal. |
| `induce.f:7878/7898`, `dfield0` / `ufield0` | Retain the existing common hooks exactly once: add the same accumulated h to both permanent-field channels; subtract beta times each channel's own trial dipole after native real/PME/self field evaluation. No additions inside lower kernels. |
| `induce.f:8245`, `localresponsespectrum` | Keep dense certification for the existing small vacuum route. The periodic route needs the separately reviewed scalable certificate on the exact active operator; never call this dense routine on the full box or bypass it merely because PCG converges. |
| `induce.f:8381`, `localresponsegradient` | Generalize its two-pass contraction to image entries. Reuse primary u/v for every image of that O. Scatter each image gradient to the primary O and its opposite to the primary metal, and accumulate image-displacement virial. |
| `epolar1.f:19/2488`, `epolar1` / `epolar1a` | Move the existing vacuum-only correction call to one common post-dispatch hook after native frame/reciprocal derivatives; remove the old inner call. Periodic c/d branches currently lack this correction. Ensure it executes once for every supported energy-gradient call. |
| `epolar3.f:19/570`, analysis dispatch / `epolar3a` | Move the pairwise analysis correction to one common hook for all supported pairwise branches, removing the old a-only addition. Add only `−K/2 sum_i u_i·h_i` to `ep` and matching `aep` entries. Guard it by the pairwise route. |

The default `epolar.f:1797` / `epolar0e` and gradient energy `epolar1e` already contract `u*udirp/alpha`, where `udirp=alpha(Ep+h)`. They require **no extra H energy or additional h energy**. `epolar3e` similarly already includes h if that alternative is used. The initial implementation must keep these route distinctions explicit to prevent double counting. The active `moments.f` half-sum readback remains appropriate; periodic lab dipoles need a declared representative convention and are not a sampling observable certified here.

## 3. Dual response and full-system roles

With `K=electric/dielec`, preserve

```
D = A_native_periodic + H;  D u = Ed+h;  D v = Ep+h
U_pol = −K/2 uᵀ(Ep+h)
delta_gradient = K/2 [uᵀ H′ v − (u+v)ᵀ h′].
```

Never assume Ed=Ep or u=v. Native alpha, direct guesses, preconditioning, PCG/peek and damping remain unchanged initially; H changes the operator through the common mutual-field hook. Native permanent frames, their torques, PME/self terms and salt-induced response use the newly solved u/v. Adding the correction to baseline dipoles would be wrong.

Bind the independently prepared namespace: metal 500/501→900/901, O36→402, H37→403. Validate reciprocal intact-water graphs and actual role/parameter records; do not identify water solely by a numeric type. Na405 and Cl414 receive no h/H but retain native permanent and induced interactions. Nonpolar/zero-alpha sites remain outside the induced unknown set and still contribute their native permanent fields. The stability certificate operates only on active positive-alpha DOFs. Keep one metal and all selected response sites active. Type403 here means water H: the separate donor prototype's O403 API must not be imported or enabled.

The candidate-specific BIO overlay remains unbound until terminal selection/admission. It must supply the same a,b,c,d and fixed damping as its immutable candidate; this plan chooses no coefficients. For conducting PME, assert the actual native boundary is `TINFOIL`. In this source, merely providing `EWALD-BOUNDARY` sets `VACUUM` regardless of its text; do not render a misleading `EWALD-BOUNDARY TINFOIL` line.

## 4. Image derivatives, reactions and cell work

Use the reviewed definition's O(number of image pairs) contraction, with `w_I=u_O+v_O`, `z_I=3 phi_I r_I·w_I`, and `C_J=(a+b eta_J)z_J + sum_(I≠J) b f_I z_I/(1+q_I)²`. Exclude exactly image J, not all images of its primary water. Prefix/suffix sums preserve this distinction.

For each image J,

```
G_J = K/2 [(u_M·v_M) beta_,J grad f_J
            − C_J grad f_J
            − 3 f_J(a+b eta_J) (phi_J w_J
                       +(phi′_J/r_J) r_J(r_J·w_J))].
```

Store `beta_,J` or its per-image vector. The current `lrdb(:,O)` representation cannot be aggregated and then reused for each image: that would multiply/cross-count stiffness derivatives. Preserve the analytic contact limit `phi′/r → 0` and the admitted series, with no clipping.

For coordinate changes, `dr_J=dR_O−dR_M+(dL)n_J`. Add G_J to O, −G_J to M. For homogeneous strain at fixed fractional coordinates, `dr_J=E r_J`; add `sym(r_J ⊗ G_J)` to native virial exactly once. This includes the lattice-image term. Using primary O−M displacement alone gives the wrong cell derivative. Retain all native virial/frame terms. The native `lattice.f:132` stores lattice vectors as **rows**, so the Cartesian shift is `transpose(lvec)*n`; the initial orthogonal case hides this orientation issue but the implementation must state it.

## 5. Wrapping and rebuild timing

Native `bounds.f` wraps whole molecules using their centers of mass. The infinite response is invariant under this reindexing, but a newly chosen literal −1…1 cube is generally a different finite sum. Preserve the initially retained **physical** image labels across a pure wrapping operation: if `R_O′=R_O+L k_O` and `R_M′=R_M+L k_M`, update each label to `n′=n−k_O+k_M`, leaving r_J unchanged. Track whole-water/metal integer image counters; do not independently wrap O versus H or apply `image()` to each response displacement.

This convention preserves exact finite-sum energy/force continuity under coordinate bookkeeping. It does not solve unlimited diffusion: retained physical images may eventually cease to approximate the infinite environment at the old bound. Recentring the 27-image set is an approximation change, not a harmless relabeling. Until a uniform coordinate/cell domain and error-controlled update rule are declared, such a boundary must stop admission rather than silently replacing the image set. Initial fixed-box checks and small coordinate/strain perturbations can use the fixed set; unrestricted sampling remains blocked.

Recompute all image geometry, h, H and derivative coefficients at every new coordinate/cell energy or gradient evaluation, including displaced testgrad calls, integrator states and cell proposals. Keep them fixed within one mutual solve. Rebuild topology/allocation only when topology or parameters change. A rejected cell move must restore/rebuild both native PME state and response state. The first implementation should favor explicit rebuilds over an unproven cache.

`udirect1` initializes spline/qfac state during `dfield0c`; current setup calls the dense spectrum routine **before** that initialization. Periodic ordering must be: current lattice/native list state → response geometry → current native direct/PME preparation and both RHS → scalable spectral admission → normal dual solve → native plus response derivatives. A certificate probe must preserve RHS and trial arrays. Invalidate it on coordinates, cell, cutoff/grid/damping or response changes; no stale PME certificate or state reused across a cell update.

## 6. Dependencies before any release

The patch is implementable without dense global A or an O(n²) derivative tensor: added geometry/gradient work is linear in 27 times the primary-water population, and H application is primary-site diagonal work. This is an operation-count statement, not a benchmark or stability proof.

Execution still needs: exact candidate/BIO composition; an explicit periodic numerical-error budget and supported finite domain; Planck/root's admitted scalable stability method under the frozen relative 1e−6 spectral criterion; the PME adjoint/even-qfac and numerical-application conditions; and a finite declared parity/derivative panel for both channels, energy routes, full Cartesian/frame forces, symmetric cell strain, image reindexing, zero bypass and salt coupling. No new thresholds or test inventory are selected here. A failed/inconclusive certificate must remain a failure, not trigger dense full-box construction, solver fallback, coefficient changes or silent cutoff growth.
