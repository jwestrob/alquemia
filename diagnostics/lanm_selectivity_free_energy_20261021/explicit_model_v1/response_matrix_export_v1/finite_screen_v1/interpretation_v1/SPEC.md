# Proposed regularized pair/collective response and short-range energy repair

**Scope:** one isolated CPU implementation of a conservative 12-parameter water
development model: `(a,b,c,d,R,epsilon)` for each metal. Eight dimensionless
response amplitudes and four existing metal–O buffered-14-7 parameters are
adjustable in a subsequently declared joint fit. Root has declared the domain
below; this specification chooses no fitted coefficients, optimizer or candidate.
It releases no execution.
It is a phenomenological hypothesis, not a unique charge-transfer or polarization
attribution. The regularized drive makes it a distinct family from the completed
finite screen; those grid predictions are not predictions of this new family.

## Fixed geometry functions and selector

Use only the explicit metal types 500/501 and intact water O type36 in the
qualified original vacuum MUTUAL route. Metal-free water, H interactions,
permanent multipoles, water alpha, original POLPAIR gamma, pdamp and mutual
tensor T stay unchanged. For each metal M define r_i=R_Oi−R_M, r=|r_i| and

```
t=2r/Angstrom
f_i=[(1+t+t²/3) exp(−t) / ((1+5+25/3) exp(−5))]²
s=sum_i f_i; Q=sum_i f_i²; g=2s/(1+s)
q_i=s−f_i; eta_i=q_i/(1+q_i)
eta=(s²−Q)/(1+s²)
```

The normalized radial shape is the already-declared f(2.50 Å)=1. The rational
environment factors are fixed prospective shapes, not fitted coordination
counts: `0<=eta_i,eta<1`, `0<=g<2`; eta_i and eta vanish for a single water.
Use sums per metal, with no global water count or outcome-selected neighbors.

Freeze each original positive metal–O gamma and pdamp before changing R/epsilon:

```
kappa=gamma0/(pdamp_M*pdamp_O)^3             [Angstrom^-3]
s3(r)=−expm1(−kappa*r³)
phi(r)=s3(r)/r³, with phi(0)=kappa
h_Oi=3 f_i (a_M+b_M eta_i) phi(r_i) r_i      [e/Angstrom²]
sigma_M=g(s) [c_M+d_M eta]
H_MM=[exp(sigma_M)−1] I3/alpha_M            [Angstrom^-3]
```

Add h to **both** native fields. All other H blocks vanish. Multiple metals
contribute additively to a water's h; each H/environment is metal-local.
Near contact, h tends `3*kappa*f(0)*(a+b eta_i)*r_i`, avoiding the diagnostic's
singular added field. Stable series evaluate phi and its derivative near zero;
this is evaluation of the same function, not clamping. All functions are
translation/rotation invariant or covariant as appropriate. Bare ions have
h=H=0; separated ligands give f→0 and cross-environment effects disappear.
Separated clusters recover additive energies asymptotically. No finite-distance
truncation or neighbor-list force discontinuity is introduced here.

H changes effective metal alpha to `alpha_M*exp(−sigma_M)` while leaving T and
native damping lengths fixed. Root declares `a,b in [−1,+1]`, `S=ln(4)`, and
box coordinates xi,zeta in [−1,+1] with `c=S(xi+zeta)/4`, `d=S(xi−zeta)/4`.
Thus `2(|c_M|+|d_M|)<=S`, bounding the local inverse-response multiplier between
1/4 and4 for **any geometry**. These are development bounds, not physical priors.
All eight response amplitudes start at zero. R/epsilon start at the exact
selected17676 values and retain their original declared intervals; POLPAIR is
frozen. Xi/zeta are fit coordinates, distinct from the induced vectors u/v below.
This bound does not guarantee positivity of the total coupled response matrix.
Require the existing normalized SPD criterion `lambda_min>1e−6*lambda_max>0`
on every evaluated configuration, with no clipping or fallback. A finite drive
and positive local diagonal do not themselves establish dynamics or global
potential stability.

## One energy and its observables

Let `At=A0+H`, `d=Ed0+h`, `p=Ep0+h`, with `A0=diag(1/alpha)−T` and
`At*u=d`, `At*v=p`. The dual functional is stationary, not a joint minimum:

```
L=K/2 [uᵀ At v−uᵀ p−vᵀ d],  K=electric/dielec
U_pol*=−K/2 uᵀp
E=E_nonpolarization,baseline−E_MO_vdw,baseline
  +E_MO_vdw(native buffered-14-7; R_MO,epsilon_MO)+U_pol*
```

The existing native buffered-14-7 term and its analytic derivative are reused,
with the same metal–O pair override convention. Other nonpolarization terms are
unchanged. The R/epsilon correction is field independent: it can compensate
pair E/F costs but cannot by itself change dipoles/χ or generate a nonzero
seven-term E3. Collective energy must come from the response functional here.

For a laboratory field F in au, add `Hartree/(K*Bohr)*W F` to both d and p,
and include `−Hartree*mu_perm,au·F`. Therefore

```
mu_au=mu_perm,au+Wᵀ(u+v)/(2*Bohr)
chi=Hartree/(K*Bohr²) Wᵀ At^−1 W.
```

Use this actual conjugate total moment, common origins and complete literal
fragment coefficients. The eight response amplitudes offer distinct algebraic
roles: a controls pair drive; b adds environment-dependent drive without changing
an isolated pair; c sets pair metal stiffness; d adds environment-dependent
stiffness. Drive leaves fixed-geometry χ unchanged. A single pair/trimer χ pair
can constrain c/d only if its sensitivity columns are independent; pair and C
moments constrain a/b after coupling to stiffness. R/epsilon act through E/F.
Observation counts do not prove joint identifiability.

## Required derivatives

For any coordinate differential, with n_i=r_i/r_i and k=2 Å^-1,

```
df_i=f_i' n_i·dr_i
f_i'=2k f_i [(1+2t_i/3)/(1+t_i+t_i²/3)−1]
ds=sum df_i; dQ=2 sum f_i df_i
deta_i=(ds−df_i)/(1+q_i)²
deta=[2s(1+Q) ds−(1+s²)dQ]/(1+s²)²
dg=2 ds/(1+s)²
dphi/dr=3[kappa*r³ exp(−kappa*r³)−s3]/r⁴
         → −(3/2) kappa² r² as r→0
w_i=3 f_i phi_i; w_i'=3(f_i' phi_i+f_i phi_i')
dh_i=(a+b eta_i)[w_i dr_i+w_i' r_i(n_i·dr_i)]
     +b w_i r_i deta_i
dsigma=(c+d eta) dg+g d deta
dH_MM=exp(sigma)/alpha_M * dsigma * I3
```

Here the plain d in `(c+d eta)` is the parameter; `dsigma`, `deta`, etc. denote
differentials. At exact contact evaluate the finite vector/Jacobian limits
without normalizing a zero vector. Every other water contributing to eta or s
receives its derivative; metal forces include the corresponding reaction.

The total polarization gradient is

```
dU=K/2 [uᵀ(dA0+dH)v−uᵀ(dEp0+dh)−vᵀ(dEd0+dh)].
```

Relative to the existing native gradient evaluated with the new u/v, the added
term is `K/2[uᵀ dH v−(u+v)ᵀ dh]`. Retain the dual cross product u*v, native
dT/dr, permanent-frame torque derivatives, and existing short-range gradients.
Parameter sensitivities use `u_theta=At^-1(h_theta−H_theta*u)` and the analogous
v equation. In particular, `chi_theta=−Hartree/(K*Bohr²) Wᵀ At^-1 H_theta At^-1 W`;
drive derivatives are zero. These provide a rank check in the eventual joint
fit without treating duplicated steps/estimators as independent observations.

## Smallest implementation and qualification step

Implement a shared h/H plus coordinate-derivative helper in a versioned copy of
the original CPU baseline, with an explicit disabled path and a scoped water
selector. Parse/read back all eight coefficients and frozen gamma/pdamp; do not
activate a global alternative water model. In `induce0a`, add h before native
udir/udirp construction and include H in every mutual solve/operator and residual
path. Keep `udirp=alpha_original*(Ep0+h)` for the energy contraction; changing
native polarity/pdamp to represent H would conflate the frozen T and energy units.
Preconditioning may remain the original positive alpha; the solved operator
must be At. Incompatible solvers/routes must fail explicitly when enabled.

`epolar0e/epolar1e` then obtain U* through the native u·udirp/alpha contraction;
do not add H energy a second time. Pair-analysis `epolar3a` needs the added
`−K/2*u·h` contribution to its native `−K/2*u·Ep0`. Analytic `epolar1a` requires
the derivative correction above. Audit all energy/analysis entry points and
ensure `moments`/field response report the conjugate dual moment. Planck's
separate native integration checklist binds precise source paths.

Before any fit, a finite source/numerical qualification must establish disabled
and zero identity, pure-water/ion invariance, actual coefficient application,
new energy/moment agreement with the same regularized saved-matrix functional,
SPD/solve residuals, contact/separation limits, and analytic versus numerical
coordinate gradients including additional-water motion and local-frame torque.
The qualification must distinguish pair/collective drive and stiffness; existing
unregularized matrix outputs or D/M gradient tests do not qualify this new code.
Root declares the exact bounded call inventory. No additional quantum campaign
is needed to construct the consumed development objective.

Use the frozen joint DATA_CONTRACT and the separately owned objective. All
original25 E/F checks remain; eight consumed E3/μ3 cycles are regression outside
the fit, not fresh validation. A finite development fit must report coupled
tradeoffs, Jacobian rank/conditioning, bounds and all failed checks. No optimizer
extension, radial-shape search, held-out tuning, periodic/GPU/protein qualification
or production release is included in this specification.
