# Compact-weight error amplification contract

**Scope:** one water-response metal environment, the frozen universal parameter box, and the admitted compact scalar errors. No selected coefficients, actual geometry, induced solution or spectral result is used here. [BOUNDS.py](BOUNDS.py) implements the formulas as pure optional-input arithmetic. Manufactured controls check the inequalities; no physical evaluation is authorized.

## Inputs and norms

Use image-labelled weights f=F and f_c=FS with Δf=f−f_c≥0; zero compact weights remain conceptually in the infinite index set. Let T=T0 bound ΣΔf. For one declared directional differential, let DΔ bound Σ|dΔf| and WΔ bound ΣΔf|dr|:

| Direction | DΔ | WΔ |
| --- | --- | --- |
| One Cartesian coordinate, or relative motion with every |dr|≤1 | T1 | T0 |
| Homogeneous strain with operator norm≤1 | V1 | V0 |

Rescale the directional bounds for a larger relative-motion/strain norm. In particular, simultaneously moving metal and oxygen in opposite unit directions can have |dr|=2. This convention cannot be silently reused for an arbitrary collective mode.

Provide compact retained sums s_c=Σf_c, D_c=Σ|df_c| and W_c=Σf_c|dr|, or upper bounds. The optional universal fallback uses N≤3980, F0=C, B=47/8 Å⁻¹:

```
s_c ≤ N F0
D_c ≤ B N F0,       W_c ≤ N F0            (unit relative motion)
D_c ≤ 21 B N F0,    W_c ≤ 21 N F0         (unit strain).
```

Do not substitute D_c≤B s_c: S′/S diverges near the outer endpoint. A direct O(Q) sum can be much tighter. The fallback needs no actual coordinates. Alpha_M>0 and kappa>0 are the unchanged native values (or suitable positive lower/upper bounds); kappa=gamma/(pdamp_M pdamp_O)^3. Let V=phi(r)r with phi=(1−exp(−kappa r³))/r³ and its contact limit. Uniformly, |V|≤K0=kappa^(2/3), ||JV||≤K1=5 kappa. Use a kappa upper bound for both.

All h errors below bound the sum of primary-water vector norms after image aggregation, hence also the full Euclidean field norm. H errors are operator norms for the one metal's scalar block. Do not multiply an H error once per image.

## Field and explicit derivative

Set A=|a|+|b|≤2 and B_b=|b|≤1. The implementation uses A=2, B_b=1 universally. Define

```
B0Δ=K0 T                         B1Δ=K0 DΔ+K1 WΔ
B0c=K0 s_c                       B1c=K0 D_c+K1 W_c
E_h = 3 [A B0Δ+B_b T B0c]
E_dh= 3 [A B1Δ+B_b T B1c
          +B_b B0Δ(D_c+DΔ)+B_b B0c(2T D_c+DΔ)].
```

The terms involving B0c/B1c are essential: modified/removed weights change the collective factor of retained waters. Proof: with e(q)=q/(1+q), |e′|≤1 and |e″|≤2, write q=q_c+Δq (excluding exactly one image), |Δq|≤T and |dΔq|≤DΔ. Expand the differences of `df (a+b e)V`, `f(a+b e)dV` and `b f e′ dq V`; apply these inequalities termwise. This proof allows attenuation of a retained weight as well as omission. The old omitted-only shortcut B1Δ=(4K0+K1)T is not valid for the switched nearest weights.

An upper bound for the compact field derivative itself is
`B_dh,c=3[A B1c+B_b D_c B0c]`.

## Stiffness

The pinned box is |a|,|b|≤1, |c|+|d|≤S/2 with S=ln4. For all nonnegative weight arrays the existing geometry bounds are |eta_i|≤2, |eta_ij|≤8, |g_i|≤2 and |g_ij|≤4. Thus

```
L1=2|c|+6|d|≤3S
L2=4|c|+28|d|≤14S
L_beta =12S/alpha_M
H_beta =4(14S+9S²)/alpha_M
E_beta =L_beta T
E_dbeta=H_beta T D_c+L_beta DΔ.
```

These follow from exp(sigma)≤4 and interpolation along f_c+theta Δf in the nonnegative orthant. Consequently ||ΔH||≤E_beta, ||dΔH||≤E_dbeta; compact |d beta_c|≤L_beta D_c. The normalized operator perturbation is at most eta_H=alpha_M E_beta. The four VDW parameters do not enter these differences: their field-independent native term is identical between these two kernels.

## Explicit correction gradient at fixed dipoles

Let U_M≥||u_M||, V_M≥||v_M||, W≥max_water||u_i+v_i|| and K_abs≥|electric/dielec|. The native stationary correction is K/2 `[u_M·v_M d beta−Σ_i(u_i+v_i)·dh_i]`, so

```
E_explicit = K_abs/2 [U_M V_M E_dbeta+W E_dh].
```

This applies separately to force coordinates and the specified strain direction. Use both independent mutual channels: replacing u·v by |u|² is unjustified. Metal reaction and image virial conventions remain those of the kernel. This comparison holds u/v fixed and does not yet bound a change in the self-consistent gradient.

## Optional conditional dipole amplification

Let P=diag(sqrt(alpha)), D_c=A0+H_c and Dbar_c=P D_c P. A separately qualified numerical lower bound ell_c for lambda_min(Dbar_c) can be supplied; this module supplies no such qualification. If ell_c>eta_H, Weyl gives Dbar_infinity≥(ell_c−eta_H)I. With alpha_max the largest active intrinsic alpha and alpha_O,max the largest corrected water-O alpha,

```
e_u ≤ sqrt(alpha_max) [sqrt(alpha_O,max) E_h
                       +sqrt(alpha_M) E_beta U_M] /(ell_c−eta_H)
e_v ≤ sqrt(alpha_max) [sqrt(alpha_O,max) E_h
                       +sqrt(alpha_M) E_beta V_M] /(ell_c−eta_H).
```

These bound full Euclidean ||u_infinity−u_c|| and ||v_infinity−v_c|| when U_M/V_M are compact-solution bounds. Missing spectral input is `not_assessed`; a nonpositive denominator is `inconclusive`, not a proof of instability. A supplied qualified lower bound uses the existing numerical qualification conventions; this task adds no compulsory interval-factorization rule. Algebraic perturbation and floating solver/application errors remain distinct, and any latter error must be budgeted by the caller.

Optionally the **explicit correction part**, now allowing changed self-consistent dipoles, is bounded by

```
E_corr_sc ≤ E_explicit + K_abs/2 [
 (L_beta D_c+E_dbeta)(U_M e_v+V_M e_u+e_u e_v)
 +(B_dh,c+E_dh)(e_u+e_v)].
```

This still omits changes to the native A0/permanent-field/frame derivative contractions caused by changed u/v. Total self-consistent force/stress admission requires those native derivative bounds and cannot be inferred from a small E_corr_sc. Multi-metal extension would need summed field bounds and a correctly assembled maximum block perturbation; the initial implementation contract remains one metal.

## Implementation and provenance

`amplify(tails, direction, kappa_upper, alpha_m, ..., dual_norms=None, spectral=None)` accepts exact decimal strings denoting upper bounds (alpha_m is its fixed positive native value). It returns outward scalar expression enclosures using the inherited Decimal helper. Supplied retained/norm/spectral inputs require caller provenance; the pure function does not certify them. `INPUTS.json` pins the original SPEC/domain, fixed regularization helper, qualified refactor derivative source, compact actual result/root admission and inherited arithmetic. Only universal/manufactured inputs are exercised here; no real alpha/kappa/dipoles or actual bound amplification is evaluated.
