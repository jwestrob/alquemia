# Finite local energy and response kernel v1

Authority: ../PHASE_1.txt and the model-form decision. This phase implements a
small differentiable model and manufactured qualification. No coefficients are
fitted; no native/quantum executable, optimizer, scheduler or periodic engine
is called. Root owns evaluation of admitted cached fixtures in integration_v1.

## Chosen finite basis

One metal and intact water, acetate and amide ligands are explicit metadata.
Acetate's two oxygens belong to one ligand: they never generate a three-body
term with each other. Terms sum over its equivalent donor oxygens. Water frame
atoms are both equivalent H atoms; acetate O frames use the shared carboxyl C;
amide O uses carbonyl C and optionally an explicitly declared plane-defining N.
No nearest-atom inference or molecular fragmentation is performed by the kernel.

At each donor O, R=O-M, r²=R·R, n=R/sqrt(r²+δ²). Unit-like covalent bond
vectors use the same δ; their sum, smoothly normalized with dimensionless
softness ε, gives b. Both H order and equivalent O order are immaterial.
The default prospective scales are δ=0.1 Å, ε=0.1, r0=2.5 Å, Rc=6 Å.
These are prototype model choices, not physical distance/accuracy thresholds.

w=(max(0,1-r²/Rc²))³ is C2 at Rc and regular at contact. rho=(r²-r0²)/r0².
Each metal/intact-ligand pair sums the following donor contributions:

* ΔV2=w[v0+vr rho+va(n·b)]. The frame is directed; H/O permutation symmetry
  does not impose an artificial head-tail symmetry on the mechanical energy.
* Δh2 on O = w[hr n+hf b].
* ΔH2 on O = w[krr nnᵀ+kff bbᵀ].

Each pair of DISTINCT intact ligands sums cross-donor contributions i,j:

* ΔV3=wi wj[v0+va(ni·nj)+vf((ni·bi)+(nj·bj))/2]. The last primitive
  permits independent mechanical changes upon ligand rotation at fixed O sites.
* Δhi=hn wi wj nj and Δhj=hn wi wj ni.
* ΔHij=kc wi wj ni njᵀ and ΔHji=ΔHijᵀ.

Seven named pair coefficients and five named triple coefficients separate
energy value/slope/orientation, drive direction and response modes. Pair
coefficients are keyed by ligand kind; triples by the unordered pair of kinds.
A caller selects the metal's parameter object. This is a fixed minimal basis,
not a fitted amplitude expansion. In particular it does not claim a complete
donor-frame tensor basis or full short-range energy surface. Defaults are all
zero. Signed stiffness coefficients are allowed only while the FULL coupled
operator passes the numerical SPD/solve checks; local block positivity is not
a substitute. No universal useful coefficient box is asserted.

The mechanical ΔV is independent of induction. h/H supply directional response
within the existing quadratic electronic class. There is no new charge-flow,
induced-quadrupole or electronic saturation variable. Pair/triple labels apply
before relaxation; observed E3 requires complete fragment cycles after solving.

## Units and dual stationary energy

All tensors are CPU float64. Coordinates Å; A0,H Å⁻³; Ed,Ep,h and the uniform
field F are e/Å²; u,v and permanent moment eÅ; ΔV and energies kcal/mol.
K is supplied literally from the native electric/dielec authority; it is never
replaced by a modern constant. Flattening is site-major xyz, with explicit atom
indices. The kernel accepts at most 200 response sites (600-dimensional dense
operator); it is not the periodic implementation.

G=A0+H; G u=Ed+h+B F; G v=Ep+h+B F, where B repeats the 3x3 identity.
U=U_other0-K μ_perm·F+ΔV-K/2 uᵀ(Ep+h+B F).
Both native fields are preserved, including unequal Ed/Ep. Induced moment is
Bᵀ(u+v)/2; susceptibility is BᵀG⁻¹B (Å³). Explicit correction gradients follow
K/2[uᵀ(dH)v-(u+v)ᵀdh]+dV, alongside ALL baseline derivatives when supplied.

`field_au_to_native` retains the historical field conversion Ha/(K a0) from
the supplied constants. `response_to_source_au` returns source-conjugate units
and separately χ/a0³; rounded native constants can make the two slightly
different. No silent conversion or constant update occurs.

A0 used by the solver must be exactly symmetric. An explicit upstream helper
`symmetrize_native_operator` accepts raw arrays only with a caller-supplied
inherited asymmetry bound and its admission authority; it returns the symmetric
part and measured raw asymmetry/change. Raw arrays/pins must remain preserved.
There is no implicit symmetrization and no default newly invented bound.
The kernel computes all eigenvalues and requires min>0 plus successful full
Cholesky. Solve residual is bounded by 100*d*eps times the declared magnitude
scale. These are ordinary-double numerical checks, not interval proofs or a
trajectory stability certificate. The eigenvalues and raw residual are returned.

## Derivative and provenance interface

`Baseline` contains A0, Ed0, Ep0, site order, K, optional nonpolar/permanent
baseline energy U_other0, optional permanent moment and authority string.
Missing U_other0 means polarization+local energy only; missing permanent moment
means induced moment only. Charged-system permanent moment needs the adapter's
explicit origin. The kernel does not establish baseline compatibility; data
adapters/root bind geometry, parameters and exact native/source provenance.

Ligand-only fragments use metal_atom=None and receive exactly zero local
correction; ion-only fragments use empty ligands. Their baseline response still
solves normally. Spectator atoms are allowed and never inferred to be ligands.

`evaluate(...,baseline_mode='frozen')` explicitly labels coordinate derivatives
partial_frozen_baseline_not_total_force even if an energy value is available.
`evaluate_provider` requires provider(xyz,cell) to supply differentiable A0,
fields, other energy and permanent moment, with an explicit completeness
assertion. This is a provider contract, not something inferred from requires_grad.
Total forces/cell derivatives are valid only for a truthful complete provider.
No cached missing baseline derivatives are fabricated or set to zero as physics.
The tests demonstrate the difference on a complete manufactured provider.

## Prospective manufactured qualification

One deterministic test invocation, one CPU numerical thread, no real fixture
replay here. Small manufactured intact ligand systems (<=17 atoms) test:
zero correction/exact baseline algebra; translation/rotation and equivalent H/O
and ligand permutations; no same-acetate triple; C2 cutoff/separation/contact;
unequal-channel energy stationarity; field autodiff and finite differences;
full coordinate and affine-cell derivatives with a complete provider; rejection
of non-SPD, asymmetric and incomplete-provider inputs; explicit partial-force
scope; and literal source unit conversions. Nonzero coefficients are synthetic
test values, not a candidate and not a calibration.

Analytic covariance/zero/dual/autodiff comparisons use 2e-10 scaled arithmetic
tolerance. Central coordinate/field/strain finite differences at 1e-5 and
5e-6 use 1e-6 scaled tolerance. C2 boundary derivatives are exactly zero at/outside
Rc and decay from inside. No accuracy threshold or data fitting follows.
TEST_KERNEL.py writes immutable RESULT.json with failures, source hashes,
thread/version settings and wall/CPU costs. Preserve failed attempts; do not
relax tolerances to obtain a pass.
