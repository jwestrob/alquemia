# Finite screening with qualified baseline response matrices

**Recommendation:** first compare a short-range water driving-field correction
with an independent, environment-dependent **metal** stiffness correction. Use
the exact saved operators to measure leverage, stability and coupled energy/
dipole costs before implementing native derivatives. This is a four-amplitude
diagnostic family, not a fitted or demonstrated repair. Execution waits for the
La-dimer scout and the complete 54-state operator qualification.

**Known obstruction:** one signed stiffness amplitude moves Dy dimer and trimer
susceptibilities in the same direction, whereas their observed offsets require
opposite directions for exact agreement. Drive cannot change susceptibility.
Do not run a fit to rediscover this; the four amplitudes diagnose leverage and
tradeoffs, with a six-amplitude collective extension only a conditional option.

## Evidence and scope

The consumed inventory retains ten energy and fifteen force checks, 34 literal
pair dipoles, four C/P/B cycles, eight E3/mu3 cycles, and four independent state
susceptibilities. The latter's two steps, two estimators and derived increments
are consistency/dependency information, not 28 independent targets.

The oxygen-radial saved Jacobian is locally weak: La-dimer susceptibility needs
a change near -3.347 a0^3, versus the declared first-order envelope about 0.113.
Its metal-site contribution is 10.853 of the model's 18.820 a0^3. That is a
partition of this model's response, not a quantum atomic attribution. Zero-field
pair and collective dipoles remain too small, while added-water susceptibility
is too large. A global polarizability increase therefore has the wrong combination
of effects. The failed EXCHPOL and D/M probes establish their tested tradeoffs;
they do not exclude a different response functional.

Use all 54 frozen metal-containing/ion geometries and the 46 cached metal-free
terms. Reuse the four field-state aliases. Preserve original PBE0-D4 energy/force
targets separately from electronic common-parent-basis E3 and moment targets;
do not substitute CP or correlated-method values. Reserved outcomes remain excluded.

## Minimal family and prospective domain

There are two dimensionless amplitudes per metal, `(a_La,l_La,a_Dy,l_Dy)`.
All native charges, multipoles, water polarizabilities, Thole lengths, mutual
tensor T, geometries and nonpolarization energies stay fixed. For each water
oxygen, let r and n point from the metal to that O. Retain the previously
declared radial shape

```
f(r) = [(1+t+t^2/3) exp(-t)]^2 / [(1+5+25/3) exp(-5)]^2,
t = 2 r/Angstrom;  f(2.50 Angstrom)=1.
h_O = a_M f(r) (3/r^2) n;  h_other=0.
s = sum_water_O f(r);  g(s)=2s/(1+s).
H_MM = [exp(l_M g(s))-1] I_3 / alpha_M;  H_other=0.
```

h has native field units e/Angstrom^2; H has Angstrom^-3. Alpha_M is the actual
baseline 1.6597/0.7280 Angstrom^3. The new g is a **prospective phenomenological
shape**, not a physical derivation: it is zero for a bare ion, one for the
2.50-Angstrom dimer, and smoothly bounded below two. It preserves the separated
ligand limit and metal-free water. It changes the metal's environment-dependent
effective local alpha to `alpha_M exp(-l_M g)`; unchanged parameter-file bytes
must not be described as unchanged effective response. Unlike the old native
POLARIZE-alpha probes, this leaves pdamp and T fixed.

A transparent proposed screening domain is `|a_M|<=0.1`, `|l_M|<=log(2)`, with
`|l_M g(s)|<=log(2)` on **every** inventory state: a local metal-response factor
between one-half and two. These are finite diagnostic bounds, not fitted
coefficients or physical priors. Root will declare the actual finite evaluation
points and numerical conditioning criterion after operator admission. No loss,
weights, physical equivalence threshold or parameter selection is defined here.

## Exact fixed-geometry predictions

Let K=electric/dielec, B=0.529177210544 Angstrom, and H_a=627.509474063 kcal/mol.
Use actual native constants/pins. Let W repeat a three-component laboratory
field at every induced site, with `c=H_a/(K B)`. At applied field F in atomic units,

```
A_t=A_0+H;   d_t(F)=Ed_0+h+c W F;   p_t(F)=Ep_0+h+c W F;
A_t u=d_t;  A_t v=p_t;
L=K/2 [u^T A_t v-u^T p_t-v^T d_t];
U_pol=-K/2 d_t^T A_t^(-1) p_t.
```

L is stationary, not a joint minimum in independent u/v. Report the exact
conjugate total laboratory moment and susceptibility:

```
mu_au = mu_perm,au + W^T(u+v)/(2 B);
chi_au = H_a/(K B^2) W^T A_t^(-1) W.
E_t(F) = E_0(0)+U_pol,t(F)-U_pol,0(0)-H_a mu_perm,au dot F.
```

Use the native dual solves even where Ed=Ep is proven for intact waters. These
equations retain charged-state origins and their linear field energy; neutral
cycles use the literal recorded coordinates and coefficients. Compute all C/P/B
and seven-term E3/mu3 values from constituent predictions, never from a sum of the
eight sampled E3 values. Retain both field-step/estimator reference values and
their limitations; the candidate is exactly linear in field at fixed geometry.

For any parameter q, exact sensitivities are
`u_q=A_t^-1(h_q-H_q u)` and `v_q=A_t^-1(h_q-H_q v)`. Thus drive changes zero-field
moments/energies but has zero susceptibility derivative; stiffness changes both,
with `chi_q=-H_a/(K B^2) W^T A_t^-1 H_q A_t^-1 W`. Positive H decreases each
same-state susceptibility, although increments need not decrease. This separates
the two mechanisms without assigning the quantum discrepancy uniquely to induction.

## Decision, identifiability and limits

First require zero-parameter replay of native u/v, energy, moments and saved field
columns, both export solves, exact site order, and unmodified matrix symmetry.
For every finite point retain the spectrum and condition number of
`diag(sqrt(alpha)) A_t diag(sqrt(alpha))`, successful positive-definite factorization,
and residual-based solution-error bounds. Positive local inverse alpha alone
does not establish stability; never repair a failed matrix by symmetrization,
eigenvalue clipping or a solver fallback.

The four-amplitude susceptibility Jacobian has zero drive columns and rank at
most two. Independent zero-field dipoles and energies must identify drive; report
column scaling/rank rather than assume identifiability from observation counts.
There is already one exact directional limitation: Dy dimer is about 0.232 a0^3
below reference while Dy trimer is about 0.905 above. A single signed stiffness
amplitude cannot move both toward their exact values, and drive cannot alter chi.
Those offsets are not established physical-error tolerances. If the joint evidence
requires a separate collective direction, the smallest nested extension adds one
metal-stiffness coefficient per metal multiplied by a fixed smooth other-water
factor that vanishes for a single water. Its shape/domain must be separately
declared; it is not part of this screen and is not presumed necessary physics.

A useful result would show stable, appreciable leverage with compatible changes
across pair, trimer, full-shell energies/moments and all four susceptibilities.
Failure means this finite family/domain cannot reconcile those consumed constraints,
not global impossibility. Large moment gains accompanied by energy deterioration
are a rejected tradeoff, not a selected candidate. New response observables remain
descriptive until a prospective development objective is declared.

Fixed-geometry matrices provide no dA/dR, dEd/dR or dEp/dR. The future force is the
stationary derivative `K/2[u^T A'_t v-u^T p'_t-v^T d'_t]`, including derivatives of
f, g, radial directions and native permanent-frame torques. All fifteen original
force checks remain required and **unevaluated for this family**; saved baseline
forces cannot certify them. No native correction, fitting, molecular call,
GPU/periodic/protein qualification or production release follows from this plan.

## Pinned evidence

- Consumed TARGET_INVENTORY.json: `61d73a16` (full hashes and source paths in its artifacts table).
- Operator FIXTURE_INVENTORY.json: `0348cc547ca3cc37ae17341aa38fc345946744d9e83f07cb77bc2f731e8fe5ff`.
- Saved oxygen-Jacobian actual RESULT: `8642b3b22c4d57f88c301df965ee947992c1418c6a2e7f39f94bc4ac46eaf5fd`;
  independent actual review `694a9143185f4335cf9807b1182e85af31c29d05d0c81b0a2a3496edb67610bd`.
- Susceptibility comparison RESULT: `41dd4eb34c99d38295ac2efadb0ef63b2607e917c6efbe9a4c64d0bae36c5d4d`.

Preparation only; zero matrix screening, fitting, compiler or molecular calls.
