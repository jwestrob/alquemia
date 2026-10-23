# Exact local response in explicit-water sampling

## Recommendation and present status

**Retain native AMOEBA/PME and its two-channel induced-dipole solver; add the
reviewed local field, metal diagonal block, and their conservative derivatives.**
The Hamiltonian does not require a dense global response matrix. First extend
the qualified CPU route to a neutral, one-metal periodic water box; port that
same route to Tinker9 afterward. This is conditional implementation design,
not evidence that the fitted model transfers or that sampling is released.
No fit outcome was read, no coefficients changed, and no engine, build or
molecular program was executed for this report.

The current implementation is numerically qualified only for the declared
1–28-site, one-metal/intact-water, nonperiodic CPU MUTUAL route. The installed
Tinker9 AMOEBA route has the necessary field/solver/PME infrastructure but no
`POLAR-LOCAL-RESPONSE` implementation. Existing HIPPO exchange-polarization
kernels are not a shortcut: their charge-penetration Hamiltonian differs.
The older EXCHPOL portability report is source evidence, not qualification of
this new regularized correction. Exact inspected files are pinned in
`SOURCE_EVIDENCE.json`.

## Same energy, without a dense matrix

Preserve intrinsic alpha, pdamp, POLPAIR, all native permanent multipoles,
scalings and mutual tensor T. Let

```
A0 = diag(1/alpha) - T;  D = A0 + H
d = Ed + h;  p = Ep + h
D u = d;  D v = p
L = K/2 [uᵀ D v - uᵀ p - vᵀ d],  K = electric/dielec
U_pol = -K/2 uᵀ p
delta_gradient = K/2 [uᵀ H' v - (u+v)ᵀ h'].
```

This dual functional is stationary, not jointly minimized. Use `uind=u`,
`uinp=v`, `udir=alpha*d`, `udirp=alpha*p`. Add h to both direct fields after
the native field calculation; return `T*x-H*x` from every mutual-field action
in both channels. This includes initial residual and PCG actions. Preserve
the final peek and check its actual D residual. Original positive diagonal
preconditioning can remain an approximation; changing alpha would also change
the model. Build h/H once per geometry and rebuild them after any coordinate
or cell change. Predictor history, if later enabled, supplies guesses only.

H is three diagonal entries per metal, so applying it costs O(number of
metals). Field construction costs O(Q), with Q the retained metal–water image
pairs. Global PME and self-consistent coupling to every other polarizable atom
remain native. A correction evaluated afterward on the old dipoles is not
equivalent. Dot-product energy routes already include h through udirp;
pair-analysis routes require only `-K/2*u·h` once, not another H energy.
Electrical observables must use permanent moment plus `(u+v)/2`, with native
unit conversion and a consistent periodic origin/branch. Do not assume u=v
in an arbitrary host. Matrix export stays optional on small diagnostic frames.

### Remove the dense coordinate Jacobian exactly

The present `lrdh(3,3,n,n)` and dense eigensystem are diagnostic implementations,
not requirements of the functional. The force can be contracted in two passes
per metal. With the unchanged SPEC shapes, set

```
r_i = R_Oi-R_M; q_i = sum_(j!=i) f_j; eta_i = q_i/(1+q_i)
w_i = u_Oi+v_Oi; z_i = 3*phi_i*(r_i dot w_i)
Rsum = sum_i b*f_i*z_i/(1+q_i)^2
C_j = (a+b*eta_j)*z_j + Rsum - b*f_j*z_j/(1+q_j)^2

G_j = K/2 { (u_M dot v_M)*(d beta/d f_j)*grad(f_j)
            -3*f_j*(a+b*eta_j)*[phi_j*w_j
                         +(phi'_j/r_j)*r_j*(r_j dot w_j)]
            -C_j*grad(f_j) }.
```

Here beta and its derivative are exactly those in the reviewed integration
checklist. G is an energy gradient: add it to O_j, its opposite to M, and
`sym(r_j outer G_j)` to the native virial. This includes every other-water
cross derivative without storing an n-by-n Jacobian. Use stable positive
prefix/suffix sums for other-water terms and the existing contact series for
phi, phi'/r and f'/r; never normalize a zero displacement. Native local-frame
torques still act through the new u/v. The added global-vector field introduces
no independent permanent-frame torque. Multiple metals can add fields/forces
to a shared water, although their actual qualification is outside the minimal
one-metal extension below.

## Source implementation map

Paths in this table are relative to the pinned CPU `source/` or Tinker9 tree.

| Route | Minimal change and retained work |
|---|---|
| CPU `polar.f`, appended helpers in `induce.f` | Replace fixed 28-site/single-metal storage by selected-site maps and O(Q) buffers. Preserve exact absent/zero bypass. Permit nonselected host atoms without applying the correction to them. Respect native nonpolarizable-site masks rather than requiring every atom to have positive alpha. |
| CPU `induce0a`, `dfield0a/b/c`, `ufield0a/b/c` | Move the current vacuum-only additions to common post-dispatch wrappers, so list and complete Ewald fields receive h and -H*x once. Remove the old inner hook to prevent double counting. Keep MUTUAL; do not enable OPT/TCG/IELSCF or mixed CHGPEN/EXCHPOL as a side effect. |
| CPU `epolar0e`, `epolar1e`; `epolar3a/b/c/d/e` | Reduced-energy paths already use corrected udirp. Extend the single h correction only to pair-analysis routes and retain component accounting. Preserve native Ewald self/reciprocal terms; audit the selected surface boundary explicitly. |
| CPU `epolar1a/b/c/d` | Call the O(Q) gradient/virial contraction once after the selected native route. Preserve native reciprocal, self and permanent-frame torque derivatives evaluated with new u/v. |
| CPU `moments.f` | Retain dual-average conjugacy; add periodic molecule-unwrapping/origin conventions to diagnostics rather than subtracting charged COM moments. |
| Tinker9 `src/amoeba/field.cpp` | Common `dfield`/`ufield` wrappers already surround Ewald/non-Ewald work. Add h after native/external direct fields and -H*x after complete mutual fields. No additional PME mesh is needed for these short-range functions. |
| Tinker9 `src/acc/amoeba/induce.cpp`, `src/cu/amoeba/pcg.cu`, `src/cu/induce.cu` | Keep two-channel PCG1 and native residual/peek conventions. Upload fixed coefficients/site maps, build per-geometry h/H on device, and use the common field hooks. New O(Q) reductions/scatters need CUDA and OpenACC implementations; qualify only the backend actually built. |
| Tinker9 `src/amoeba/epolar.cpp`, `src/cu/amoeba/epolar.cu`, `src/acc/amoeba/epolarnonewald.cpp`, `epolarewald.cpp` | Add pair-analysis h energy only where needed and the explicit local gradient/virial once. Preserve native dot-product energy, reciprocal/self derivatives and torques. |
| Tinker9 `src/amoeba/emplar.cpp`, `src/cu/amoeba/emplar.cu` | The fused CUDA route independently calls induce and dot-product energy. It needs the same one-time local force/virial dispatch; otherwise a standalone epolar patch silently misses the sampling route. |

Type 500/501–36 selection, intact water H connectivity, metal +3 monopole,
fixed kappa and original water parameters remain explicit. Allowing unrelated
counterion/protein atoms through the guard preserves their original Hamiltonian
terms but changes their induced response self-consistently. It does not qualify
donor/protein energetics or release held validation. Other water types need an
explicit verified mapping, not automatic class substitution.

## Two consequential production gaps

**Periodic extension and tails.** The frozen f is exponentially decaying but
nonzero at every finite separation. A minimum-image hard cutoff changes the
Hamiltonian and can make the added vector field discontinuous when an image
changes. Define the periodic model explicitly as converged lattice-image sums
of the same functions: each primary metal's environment sums water images,
and each primary water receives contributions from metal images. Use the same
image displacement in energy, force and cell-strain derivatives. These sums
converge absolutely for a nondegenerate cell. A finite image/neighbor list needs
a declared residual-tail bound for fields and derivatives, checked against an
expanded-image reference. A smooth compact-support replacement would be a new
model and is not proposed as an unnoticed optimization. For a few metals,
visiting all real waters is already O(M*Nwater), not O(Nwater²). Begin with
matched conducting-boundary neutral boxes; other Ewald surface conventions and
charged-cell backgrounds require their own consistent audit.

**Stability certification.** The local inverse-response factor in [1/4,4]
does not prove D is positive definite. Current per-solve dense probing/Jacobi
must be replaced explicitly before scalable sampling. PCG convergence,
positive curvature on visited directions, or unconstrained Lanczos Ritz values
alone do not certify the full spectrum. A useful exact reduction, conditional
on baseline A0 being SPD, is

```
D = A0 + U B Uᵀ; G = Uᵀ A0^(-1) U
D is SPD iff I + G^(1/2) B G^(1/2) is SPD.
```

U selects the 3M metal coordinates; at most 3M baseline matrix-free solves
and a small matrix replace the dense update test. This works with signed B,
including zero entries. It proves positivity, not the stronger frozen
`lambda_min(Dbar)>1e-6*lambda_max(Dbar)>0`, where
`Dbar=sqrt(alpha)*D*sqrt(alpha)`. Preserving that criterion requires certified
baseline/extremal bounds or another reviewed matrix-free certificate; ordinary
estimated eigenvalues cannot silently replace it. Retain the dense criterion
on small parity fixtures, stop on unproved/failed stability, and introduce no
clipping or fallback. This remains a real implementation prerequisite.

## Finite conditional qualification, before sampling claims

These are bounded development stages, not submitted jobs or bound inputs.
Use frozen diagnostic coefficients, not fit-outcome-selected favorable frames.

1. **Vacuum equivalence:** reuse the existing 57-program/636-internal-FD panel
   and its admitted raw reference for each new backend (CPU refactor first,
   GPU second). Keep its absent/zero/ion/water identities, both metals,
   separate and combined response, Cartesian and rigid-water derivatives,
   dual-field conjugacy and matrix checks. No new QM or new accuracy threshold.
2. **Periodic static and virial:** use the two existing neutral explicit-water
   endpoint boxes, with deliberately matched parameters, PME/cutoff/surface
   settings and verified water mapping. Off/on CPU/GPU E/F/virial is 8 calls;
   whole-water reimaging of each enabled box on both backends is 4 more.
   CPU central differences along three fixed coordinate directions per metal
   (metal translation, rigid-water rotation, other-water motion) need 12
   energies; six symmetric cell strains per box need 24 energies: **48 calls**
   total. For a finite image-tail gate, use three prospectively fixed image
   sets `|n_x|,|n_y|,|n_z|<=1,2,3`; the first is the enabled static call and
   the other two add 8 calls across both metals/backends, for **56 calls**.
   Require a residual-tail bound as well as observed convergence; failure is
   inconclusive, with no automatic image expansion. Freeze steps and
   numerical/strain budgets before release. Do not substitute all-Cartesian
   FD on the 11,958-atom boxes. Native E/F parity uses
   inherited 1e-4 kcal/mol and 1e-3 kcal/mol/Å tolerances; virial budgets must
   be declared from strain/numerical precision, not invented after outcomes.
3. **Integrator gate:** after static qualification, off/on, La/Dy, CPU/GPU,
   and dt/dt2 gives 16 finite NVE segments. A proposed bound is 100 steps at dt
   and 200 at dt2 for the same duration: 2,400 force updates. Bind the existing
   sampler timestep, identical velocities/seed/constraints and residual rules
   prospectively. Examine energy drift/timestep scaling and induction failure;
   a thermostat must not conceal force inconsistency. No long trajectory or
   thermodynamic-accuracy claim follows from this gate.

The smallest implementation advance is therefore a CPU common-hook/O(Q)
refactor, followed by a declared periodic image convention and certified
stability strategy. GPU kernels follow the same functional once that reference
exists. This preserves the current model rather than using a new polarization
model to bypass missing implementation work.
