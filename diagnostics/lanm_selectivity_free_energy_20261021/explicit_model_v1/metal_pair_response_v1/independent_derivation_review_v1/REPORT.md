# Direct/mutual pair damping: independent derivative review

**The proposed dual-field stationary energy gives the required derivative.**
This is a source-only review of the existing nonperiodic, double-loop MUTUAL
Thole route, with fixed scalar polarizabilities and reciprocal pair kernels.
It does not qualify a new executable or a parameter value.

Let `u = uind`, `v = uinp`, `Ed = field`, `Ep = fieldp`,
`A = alpha_inverse − T_mutual`, and `K = electric/dielec`. For symmetric A,

```
U(u,v;R) = K/2 [uᵀ A v − uᵀ Ep − vᵀ Ed]
∂U/∂u    = K/2 (A v − Ep)
∂U/∂v    = K/2 (A u − Ed)
```

At the two converged native solves, `A u = Ed`, `A v = Ep`, reciprocity gives
`uᵀ Ep = vᵀ Ed`; therefore `U = −K/2 uᵀ Ep`. This matches `epolar1e`:
`udirp = alpha*fieldp`, so `−K/2 Σ uind*udirp/alpha` is the same energy.

For any Cartesian coordinate x, stationarity cancels the implicit dipole
derivatives:

```
dU/dx = K/2 [uᵀ (∂A/∂x) v − uᵀ (∂Ep/∂x) − vᵀ (∂Ed/∂x)]
      = −K/2 [uᵀ (∂T_mutual/∂x) v
               + uᵀ (∂Ep/∂x) + vᵀ (∂Ed/∂x)]  # fixed alpha
```

The physical force is minus this derivative. Native `dep` accumulates the
energy gradient; do not introduce an extra sign into its pair updates.
Incomplete induction leaves the residual terms
`K/2 [(Av−Ep)ᵀ du/dx + (Au−Ed)ᵀ dv/dx]`, so both native solves must converge.

## Exact native contraction map

Locations refer to the pinned original `epolar1.f` / `induce.f`.

| Native location | Required coefficient family |
| --- | --- |
| `induce.f:866–925`, `dfield0a` | Direct D kernel for both permanent fields; retain their distinct `dscale` and `pscale` exclusions. |
| `induce.f:1394–1430`, `ufield0a` | Mutual M kernel for both induced fields, retaining `uscale`. |
| `epolar1.f:450–522` | D factors `sc3/sc5/sc7`, their coordinate derivatives `rc3/rc5/rc7`, and derived `dsr*`/`psr*`. |
| `epolar1.f:540–619` | Direct D factors for charge-flux potentials and permanent dipole/quadrupole torque fields. A guard may exclude charge flux from the first implementation. |
| `epolar1.f:624–716` | D factors and derivatives in the permanent-field spatial tensors. |
| `epolar1.f:848–857` | `dEd/dR` contracts with **uinp** and `dscale`; this uses D, not M. |
| `epolar1.f:861–869` | `dEp/dR` contracts with **uind** and `pscale`; this also uses D. |
| `epolar1.f:919–951` | Reset `sc3/sc5/rc3/rc5` to M before the `uscale`-weighted cross contractions of uind and uinp. |
| `epolar1.f:1205–1229` | Retain balanced site-gradient and virial accumulation unchanged. |

The prefactor `f = K/2` is already present at `epolar1.f:237`. Both dual-field
direct contractions are necessary; dropping either or adding another half
would change the energy derivative.

The stock alternate-direct reset at `epolar1.f:898` uses
`min(pti,thole(k))`. That is not the authoritative POLPAIR lookup. A new active
pair must use its exact `thlval(jpolar(i),jpolar(k))` baseline for both
`gammaD = gamma0*(1+deltaD)` and `gammaM = gamma0*(1+deltaM)`, including the case
`deltaM=0` with nonzero `deltaD`. Never carry D factors into the mutual block.
Absent/zero and unmatched pairs must retain the original path.

For cubic damping `t=gamma*(r/(pdamp_i*pdamp_k))^3`, the existing factors are
`s3=1−exp(−t)`, `s5=1−(1+t)exp(−t)`, and
`s7=1−(1+t+0.6t²)exp(−t)`. The native recurrences satisfy
`rc3=rr3*grad(s3)`, `rc5=rr5*grad(s5)`, and `rc7=rr7*grad(s7)`.
Use the selected D or M gamma consistently in both factors and derivatives.

## Limits of this derivation

The two-variable cross functional is stationary, **not a joint minimum in
independent u and v**; its Hessian has positive and negative branches. Positive
definiteness is a requirement on the physical response operator A, not on this
dual-field functional. Positive damping coefficients alone do not prove A is
positive definite for every geometry.

The derivation supports the proposed narrow implementation architecture. It
does not establish force correctness for a new build, arbitrary modes, periodic
routes, OPT/TCG, charge penetration, EXCHPOL, or charge flux. The initial source
must guard unsupported branches and retain zero-edit behavior. No engine,
build, parameter probe, molecular calculation, or finite-difference test was
performed in this review.
