# Conservative native integration review

Source-only review of the original CPU vacuum MUTUAL path and the proposed
regularized eight-response family. No candidate code, build or molecular call
has been qualified here. The old unscreened finite-grid result stays distinct.

## One operator and two fields

Keep intrinsic `polarity`, `pdamp`, original mutual T and POLPAIR coefficients
unchanged. Write `D=A+H`, `d=Ed+h`, `p=Ep+h`, where A=diag(1/alpha)-T. Solve
`D u=d` and `D v=p`. Native arrays remain `udir=alpha*d`, `udirp=alpha*p`,
`uind=u`, `uinp=v`; never replace alpha by its effective response value in the
parameter reader, which would also change damping lengths.

- Add h to both direct fields before forming udir/udirp in induce0a.
- The guessed initial residual needs `-H*u` and `-H*v`; every PCG action needs
  `+H*conj` and `+H*conjp`. Equivalently one common mutual-field wrapper can
  return `T*x-H*x` for both channels. Implement exactly one route, not both.
- Original positive preconditioning can remain an approximation. It must not
  silently redefine alpha/T. Preserve both channels, residual stop and final
  peek. Check the actual post-peek residual for D, not the old A; its bound uses
  `I-peek*D*diag(alpha)`. Signed H requires actual total-D stability separately.
- Recompute geometry functions at every displaced geometry and every relevant
  native entry. An analysis EM run induces twice. Do not accumulate h/H across
  solves, cache stale coordinates or patch only the first call.
- Operator export must describe the actual D and d/p, including the explicit
  self block and added fields. The old baseline export's direct-field replay
  and mutual-only probe cannot be relabeled as a candidate operator unchanged.

## Energy routes and the factor-of-two trap

The stationary functional is
`L=K/2[u^T D v-u^T p-v^T d]`, K=electric/dielec. It is stationary in two
independent channels, not a joint minimum. At convergence,
`U_pol=-K/2*u^T p`; no separate `+K/2*u^T H v` is added to this reduced energy.

Native energy uses epolar0e; gradient energy calls epolar1e from epolar1a.
Both already contract `-K/2*uind*udirp/polarity`. With the corrected udirp they
need no extra h or H energy. Native analysis instead selects pairwise epolar3a,
which sums `-K/2*u^T Ep`; it needs **exactly** `-K/2*u^T h` once, with consistent
atom/component accounting. Adding the full stationary-functional coupling
`-K/2*(u+v)^T h` to this reduced energy would double count part of the response.
Cross-check all three routes. Do not hide a difference by adjusting the cached
nonpolarization energy. Other pairwise/list/periodic branches are unqualified.

## Explicit gradients and native torques

Keep the existing native derivative contractions of A, Ed and Ep evaluated at
the new u/v, including dscale/pscale/uscale and permanent-frame torque transfer.
Add once to the **energy gradient**

`delta_grad = K/2 [u^T H' v - (u+v)^T h']`.

The force is its negative. For an isotropic metal block H=beta I, the stiffness
contraction is `(u_M dot v_M)*grad(beta)`, not `|u_M|^2` in a general dual system.
The added radial field depends on global geometry, not permanent local axes;
it creates no independent permanent-multipole torque term. Existing frame
torques still change through u/v and must remain intact, including hydrogen
forces. Other-water screening creates cross-water gradients even though h is
stored only on O sites.

Let `v_j=R_Oj-R_M`. If `G_j` is the explicit correction differentiated with
respect to v_j, add G_j to O_j's energy gradient and `-sum_j G_j` to the metal.
Accumulate the native virial consistently with `sum_j sym(v_j outer G_j)`;
the source's pair convention uses `-r_ik outer grad_i`. Do not substitute force
for gradient in this update. Translation and rigid-rotation balance are useful
checks, not substitutes for Cartesian/local-frame finite differences.

## Smooth screened-field derivatives

Use fixed `kappa=gamma0/(pdamp_M*pdamp_O)^3` from the original resolved native
charge-field damping, `s3=-expm1(-kappa*r^3)`, and `phi=s3/r^3`. Define
`h_i=3 B_i phi_i v_i`, `B_i=f_i(a+b eta_i)`. The regularized field tends linearly
to zero at contact: phi tends to kappa and
`phi'=3[kappa*r^3*exp(-kappa*r^3)-s3]/r^4`, tending to
`-3*kappa^2*r^2/2`. Use stable small-r series and vector formulas, not an
undefined unit direction or cancellation of two tiny numbers.

For `t=2r`, `P=1+t+t^2/3`, the normalized S2U derivative is
`f'=-(4/3) f*t*(1+t)/P`; `f'/r` tends to `-8*f(0)/3`.
For `s=sum f`, `Q=sum f^2`, `s_i=sum_(j!=i) f_j`, use

```
eta_i = s_i/(1+s_i)                 d eta_i/d f_j = (1-delta_ij)/(1+s_i)^2
eta   = (s^2-Q)/(1+s^2)             d eta/d f_j = 2[(s-f_j)(1+s^2)-s(s^2-Q)]/(1+s^2)^2
g     = 2s/(1+s)                    d g/d f_j = 2/(1+s)^2
sigma = g(c+d eta)                  d sigma/d f_j = g_j(c+d eta)+g*d*eta_j
beta  = expm1(sigma)/alpha_M        d beta/d f_j = exp(sigma)*sigma_j/alpha_M
B_i,j = delta_ij(a+b eta_i)+f_i*b*(1-delta_ij)/(1+s_i)^2
```

Form all coordinate derivatives with `grad_vj f_j=f'_j v_j/r_j`. The exact
field Jacobian is
`J_ij=3 phi_i v_i outer grad_vj(B_i) + delta_ij*3 B_i[phi_i I+(phi'_i/r_i)v_i v_i^T]`.
Then `G_j=K/2[(u_M dot v_M) grad_vj(beta)-sum_i J_ij^T(u_i+v_i)]`.
Compute other-water sums directly (and s^2-Q as twice the distinct pair products
if useful), preserving exact zero collective terms for a single water.

All shapes depend on distances and symmetric sums: translation invariance,
rotation covariance and water-permutation invariance follow. Pure water and
bare ions have no correction; separated waters give h,H tending to zero.
For nonnegative f, `0<=eta_i,eta<1`, `0<=g<2`. Root's
`c=S(u+v)/4,d=S(u-v)/4`, S=ln4 and u,v in [-1,1], gives
`|c|+|d|<=S/2`, hence `|sigma|<=S` and inverse-response multiplier [1/4,4].
This bounds the local block, not the full D spectrum. No clipping/fallback is
authorized. The later R/epsilon changes do not alter fixed kappa or these
response derivatives; their existing native derivatives must still be checked.

## Supported path and minimal actual evidence

Use an isolated original-baseline source copy. Initially restrict to the
reviewed single-metal/water, active positive-alpha, vacuum/no-list MUTUAL path.
Reject unsupported Ewald/replica/polymer/solvent, alternate-Thole/CHGPEN/EXCHPOL,
prediction/IELSCF/OPT/TCG, mixed old pair-response controls and finite cutoffs.
Static external field support requires an explicit consistent energy/moment
test; the time-dependent field route remains unsupported. Nonzero parameter
readback must prove the feature is applied, while absent/zero controls preserve
the original path. Native moments use uind; for general Ed!=Ep the conjugate
moment uses `(u+v)/2`, so intact-water equality must be verified explicitly.

The bounded qualification adapter should reuse the existing pinned inputs,
parsers and native tolerances: (1) absent/zero plus pure-water/bare-ion identity;
(2) separate drive/stiffness and mixed nonzero operator/energy-route comparisons
for both metals; (3) Cartesian finite differences on asymmetric water fixtures
with other-water cross terms, metal/O/H motion and rigid-water torque directions;
(4) one small matched static-field energy/conjugate-moment test per metal;
(5) finite dual induction and actual D stability on every evaluated control.
Freeze exact coefficients, call counts, FD/field steps and old tolerances before
execution. These are numerical checks only, with no physical pass/fit/selection.
