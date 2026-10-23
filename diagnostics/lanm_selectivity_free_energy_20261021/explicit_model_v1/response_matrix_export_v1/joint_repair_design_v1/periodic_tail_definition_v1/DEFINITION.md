# Periodic image definition and certified tails of the frozen water response

This is a symbolic, source-only definition. It selects no fitted coefficients, finite image radius, tolerance, cell, trajectory or new kernel. No numerical experiment, native call, optimization or fit outcome was read. The referenced SPEC defines the analytic functions; the native source implements their vacuum evaluation, including a small-argument series. The infinite sums below extend those functions rather than replacing them by a minimum-image cutoff.

## 1. Exactly which objects are summed

Let L be a nonsingular 3×3 cell matrix in Å, with finitely many primary metals m and intact primary water oxygens i. Current implementation scope is one metal; formulas are additive in primary metals without asserting multi-metal qualification. Use fixed molecular labels and consistent whole-water image representatives. Define

```
I=(i,n), n in Z³
r_mI = R_Oi - R_Mm + L n
f_mI = F(|r_mI|)
s_m = sum_I f_mI; Q_m = sum_I f_mI²
q_mI = sum_(J != I) f_mJ = s_m - f_mI
eta_mI = q_mI/(1+q_mI)
eta_m = (s_m²-Q_m)/(1+s_m²)
g_m = 2s_m/(1+s_m)
```

**The exclusion is exactly one image-labelled oxygen I.** Other images of the same primary water remain in q_mI. Excluding every image with primary label i would define a different functional. Consequently, a periodic cell containing one primary water generally has nonzero collective factors: its infinitely repeated waters are distinct neighbors. This does not contradict the isolated one-water limit.

For λ=2 Å⁻¹, P(t)=1+t+t²/3 and C=[P(5)exp(−5)]⁻², retain

```
F(r)=C P(λr)² exp(−2λr)
phi_m(r)=[1-exp(−kappa_m r³)]/r³; phi_m(0)=kappa_m
V_m(r)=phi_m(|r|) r
h_Oi = 3 sum_m sum_n f_mI [a_m+b_m eta_mI] V_m(r_mI)
sigma_m = g_m [c_m+d_m eta_m]
beta_m = [exp(sigma_m)-1]/alpha_m; H_Mm,Mm=beta_m I3
```

All other H blocks vanish. The original positive kappa, alpha, POLPAIR and mutual tensor remain fixed. The sum over n for M_m in the primary cell is equivalently the contribution from metal images (m,−n) to primary O_i. Each labelled metal–water-image relation is included once, with **no half-weight**. H is evaluated once per primary metal, not once per metal image.

The native periodic A0, permanent fields and boundary condition are separate, unchanged ingredients. The correction enters D=A0+H and both d=Ed+h, p=Ep+h; D u=d, D v=p. The energy is per primary cell, U_pol=−K uᵀp/2. This definition does not supply a charged-cell background or alter PME surface conventions. A future first qualification should bind the already proposed neutral conducting-boundary setting explicitly.

Whole-water reimaging R_Oi→R_Oi+L k simply reindexes n→n−k. Metal reimaging similarly reindexes its environment. The infinite functions are invariant. A numerical finite image list is only an approximation to these sums; it must not silently become a different cutoff Hamiltonian.

## 2. Conservative coordinate and strain derivatives

Differentiate every image using its same displacement. For primary coordinate variations,

```
dr_mI = dR_Oi - dR_Mm + (dL)n.
```

Every image of a moved primary water moves; q_mI excludes only df_mI, not all derivatives of that primary water. With fixed fractional atom coordinates and homogeneous infinitesimal strain E, dL=E L and dr_mI=E r_mI. Derivatives of the native periodic terms, permanent local frames and induced-solver energy remain native and use the new u/v.

One may apply the frozen SPEC differential formulas term by term. An O(number of retained image pairs) contraction per primary metal is also exact in the infinite limit. Let w_I=u_Oi+v_Oi (periodically repeated), z_I=3 V_m(r_mI)·w_I, and use the metal's own a,b,c,d:

```
Rsum = sum_I b f_I z_I/(1+q_I)²
C_J = (a+b eta_J) z_J + Rsum - b f_J z_J/(1+q_J)²
eta_,J = 2(s-f_J)/(1+s²) - 2s(s²-Q)/(1+s²)²
beta_,J = exp(sigma)/alpha *
          [2(c+d eta)/(1+s)² + g d eta_,J]
G_J = K/2 { (u_M·v_M) beta_,J grad F(r_J)
            -3 f_J (a+b eta_J) [phi_J w_J
                    +(phi'_J/r_J) r_J(r_J·w_J)]
            -C_J grad F(r_J) }.
```

G_J is an energy gradient with respect to that image displacement. Add sum_(m,n) G_m,(i,n) to primary O_i and the opposite sum over all its water images to each primary metal. This includes the reactions and all cross-environment derivatives, including the other images of the same primary water. No independent image atoms or dipole variables are introduced.

For strain, dU_correction=sum_J G_J·E r_J. Thus its work-conjugate tensor is sum_J G_J⊗r_J; the symmetric strain contribution is sym(sum_J r_J⊗G_J), consistent with the native symmetric virial convention. Keep native frame-torque, reciprocal/self and all other virial contributions. Do not infer the total force or stress by adding this correction to baseline dipoles.

The displayed kernel derivatives have finite contact limits. Use the same analytic limits/series as the frozen source, not normalized zero displacement or clipping. Absolute convergence and the bounds below justify exchanging sums and these first derivatives on a compact nondegenerate cell domain.

## 3. A computable lattice-tail certificate

Use an image cube ||n||∞≤N solely to describe numerical approximation, not a chosen physical cutoff. For each metal, let Nw be the finite number of primary water oxygens and choose certified geometric bounds

```
0<a≤sigma_min(L),  b≥||L||₂,
D≥max_i ||R_Oi-R_M||.
```

For example a=1/||L⁻¹||_F and b=||L||_F are valid mathematical bounds. Actual floating-point certificates need outward rounding/error bounds. A shell k=||n||∞ contains exactly 24k²+2 lattice vectors, and

```
max(0,ak-D) ≤ |r_mI| ≤ sqrt(3) b k+D.
```

Since F is nonnegative and monotone decreasing, define for p=0,1

```
T_p(N) = Nw sum_(k=N+1..infinity)
         (24k²+2) [sqrt(3)b k+D]^p F(max(0,ak-D)).
```

Then sum_omitted f_I |r_I|^p ≤ T_p. These bounds require no minimum metal–water distance or assumed radial density, and do not depend on cancellation between opposite images. They may be conservative for skew cells.

The remaining infinite scalar sum has an explicit finite geometric-series bound. Choose integer K> D/a and K≥N+1; evaluate the finitely many earlier shells directly. For k≥K, with rho=exp(−2λa),

```
F(max(0,ak-D)) ≤ C exp(2λD) P(λak)² rho^k.
```

The resulting shell polynomial has nonnegative coefficients and degree 6+p. For every needed integer d≥0,

```
sum_(k=K..infinity) k^d rho^k
 = rho^K sum_(j=0..d) choose(d,j) K^(d-j)
       sum_(l=0..j) Stirling2(j,l) l! rho^l/(1-rho)^(l+1).
```

This is a finite, positive-term expression with a rigorous upper evaluation available by directed rounding. Numerical overflow/underflow must not be interpreted as a zero tail; scaled/logarithmic interval evaluation is appropriate. No shell count or value has been tested here.

For Q, 0≤Q∞−Q_N≤F(0)T_0 (or use F² for a tighter degree-eight exponential bound). Stable positive sums/prefix–suffix sums should be used for q and s²−Q when required; this is arithmetic evaluation, not a changed definition.

## 4. Field and stiffness error bounds, including retained-neighbor changes

Omitting images changes the environment factors of *retained* pairs. A bound on the omitted direct field alone is insufficient.

For one metal, let S_T=sum_omitted f≤T_0 and retain any certified upper bound called T0. Let D_R=sum_retained |df_I| for the directional variation of interest, and T1≥sum_omitted |df_I|. From the frozen radial derivative,

```
|F'(r)|≤2λ F(r).
```

For one atomic Cartesian derivative |dr_I|≤1, so T1=2λT_0 suffices. For a collective coordinate use a declared uniform |dr_I|≤Lq, giving T1=2λLq T_0. For a unit-operator-norm homogeneous strain, |dr_I|≤r_I and T1=2λT_1(N). These directional choices distinguish force units from strain/virial units.

Define V0=kappa^(2/3), V1=5 kappa. Uniformly including contact,

```
|V(r)|≤min(kappa r,r^-2)≤V0
||dV/dr||₂≤min(5 kappa,4/r³)≤V1.
```

For retained pairs define B0_R=sum f|V| and B1_R=sum [|df||V|+f V1|dr|]. For omitted pairs certified bounds are

```
B0_T≤V0 T0
B1_T≤(2λ V0+V1) sum_omitted f |dr|.
```

One may use all retained/omitted waters to bound any primary-water component, or restrict the B sums to that primary water's images for sharper h bounds. D_R and T1 still include the complete metal environment. Put A=|a|+|b|. The full vector field and its directional derivative satisfy

```
||h∞−h_N|| ≤ 3[A B0_T + |b| T0 B0_R]
||dh∞−dh_N|| ≤ 3[A B1_T + |b| B0_T(D_R+T1)
                 +|b| T0 B1_R + |b| B0_R(2 T0 D_R+T1)].
```

These are triangle-norm bounds and can be summed over primary metals. They follow from |eta'(q)|≤1, |eta''(q)|≤2 and q∞−q_N=T0_actual for each retained image. They include omitted-pair derivatives and the change to retained-pair environment derivatives.

For metal stiffness, view eta=(s²−Q)/(1+s²) as a function of all nonnegative f coordinates. Along this domain, conservative global bounds are

```
|eta_,i|≤2; |eta_,ij|≤8;
|g_,i|≤2; |g_,ij|≤4.
Lsigma1=2|c|+6|d|
Lsigma2=4|c|+28|d|.
```

For clarity, the Hessian used in that bound is

```
eta_,ij = 2(1+Q)/(1+s²)² +4s(f_i+f_j)/(1+s²)²
          -8s²(1+Q)/(1+s²)³ -2 delta_ij/(1+s²).
```

Each of its four absolute terms is ≤2 because Q≤s² and f_i,f_j≤s. The frozen coefficient box gives exp(sigma)≤4. Therefore set

```
Bbeta1=4 Lsigma1/alpha
Bbeta2=4 (Lsigma2+Lsigma1²)/alpha
|beta∞−beta_N|≤Bbeta1 T0
|d beta∞−d beta_N|≤Bbeta2 T0 D_R+Bbeta1 T1.
```

The derivative bound follows by adding the omitted weights continuously from zero, bounding the change of each retained partial derivative with the Hessian. It is valid for signed c,d. H is the corresponding scalar times I3, so these are also block operator-norm bounds. Zero response amplitudes yield exact zero bounds and the original bypass.

## 5. What these bounds do and do not certify about force

At fixed finite periodic dipoles u,v, the correction's directional energy-gradient error is bounded by

```
|delta(dU_corr)| ≤ |K|/2 [sum_m |u_M||v_M| E_dBeta,m
                         + max_i |u_Oi+v_Oi| E_dh,total].
```

Here E_dBeta and E_dh are the bounds above, with aggregated primary-water triangle norm. The same formula applies to a Cartesian force component or a unit strain component using its appropriate dr bounds. It certifies the missing-image contribution and retained-environment changes in the explicit h/H force terms.

It does **not** by itself bound the difference of fully self-consistent total forces. The image approximation also changes u,v and therefore the native mutual and permanent-frame force contractions. A bound on those changes needs an independently certified inverse-response bound and native derivative bounds. Conditionally, if ||D∞⁻¹||≤G is known,

```
||u∞−u_N||≤G [||h∞−h_N||+||H∞−H_N|| ||u_N||]
```

and analogously for v; further total-force propagation must use actual native derivative bounds. The present vacuum-only dense SPD admission supplies no periodic G or uniform dipole bound. This note deliberately supplies no global stability scheme and makes no total-force guarantee from PCG convergence.

## 6. Finite-domain assumptions and numerical prerequisites

A fixed nonsingular cell, finite primary population, finite original positive kappa/alpha and finite admitted coefficients suffice for absolute convergence at a snapshot. A **uniform** certificate over NPT/strain/coordinate sampling additionally requires declared a≥a_min>0, b≤b_max, D≤D_max (or consistently wrapped representatives), bounded primary population and fixed positive damping/alpha bounds. If the cell degenerates or populations/density grow without bound, there is no useful uniform N or tail bound from the available assumptions. The existing source has no such sampling domain.

For the exact frozen box coefficients the geometric h/H and first-derivative bounds above are available. A uniform self-consistent force error still needs finite dipole/inverse-response/native-derivative bounds; those are currently unavailable. No conclusion about global potential stability or physical transfer follows.

Machine arithmetic and kernel-series errors are separate from the lattice tail. The native phi series through zeta^5, zeta=kappa r³<10^-4, has analytic remainder at most kappa*zeta^6/7!. Its phi' series through zeta^4 has remainder at most kappa² r²*zeta^5/280. Roundoff in f, exp, rational environments and shell certificates also needs a prospective numerical budget. The image tail must not absorb or conceal these errors.

The prior proposal of image cubes 1,2,3 supplies possible finite test sets, not a certificate that any one is sufficient. Root must bind an actual cell domain, requested error budget and numerical certification before any expanded-image experiment or sampling. This definition changes no running/native source and releases no test.
