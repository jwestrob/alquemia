# Independent review of the symbolic periodic tail definition

**Result: source and mathematics pass, with the stated conditional scope.** No correction to `DEFINITION.md` is requested. This review supplies no cell, image count, error budget, periodic implementation, or physical qualification. The exact definition, analytic SPEC, coefficient domain and native helper are pinned in `REVIEW.json`.

## Checked identities and bounds

1. **Image counting.** The environment index is the pair `(primary water, lattice integer)`. Excluding one such pair from q leaves the other images of the same water present. The sum into a primary oxygen is equivalently a sum over metal images with the opposite lattice integer. Stiffness is evaluated once per primary metal. There is no pair-energy half factor to insert into the driving field; the stationary polarization energy already has its specified K/2 factor. Reimaging permutes the infinite sum. None of these statements establish multi-metal native qualification.

2. **Conservative derivatives.** For each image, `dr=dR_O−dR_M+dL*n`; at fixed fractional coordinates and homogeneous strain this is `dr=E*r`. Summing image gradients onto the primary oxygen and their negative onto the metal preserves translation invariance. Differentiating `sum_I f_I(a+b*eta_I) z_I` with respect to f_J gives precisely C_J: its own prefactor term plus all other environment terms. The displayed G_J therefore equals the frozen dual-functional correction `K/2[(u_M·v_M)d beta−(u+v)·dh]`. It uses u·v and the dual sum, not a single-dipole substitution. The strain work tensor is G⊗r; symmetric strain gives the stated native symmetric virial convention. Native permanent-frame torques and periodic derivatives remain necessary.

3. **Radial tail.** For t=λr, `F'/F=−2λ*t*(1+t)/(3P(t))`. Thus F is nonnegative and decreasing, and `|F'|≤2λF`, since `t+t²≤3+3t+t²`. The exact cube-shell count is `(2k+1)^3−(2k−1)^3=24k²+2`. The singular-value and operator-norm inequalities give the stated lower and upper distance bounds. These justify both T_0 and the strain-weighted T_1 without an assumed density or minimum contact distance.

4. **Positive series evaluation.** For k>D/a, replace P(λ(ak−D)) by the larger P(λak), while retaining the exponential factor `exp(2λD)*rho^k`. The shell polynomial has degree 6+p and nonnegative coefficients. Expanding `(K+m)^d` and then `m^j=sum_l S(j,l)*(m)_l` gives the displayed finite Stirling expression exactly. Its terms are positive for `0<rho<1`. This proves a symbolic bound, not a directed-rounding implementation; inverse-cell bounds, exponential evaluation, overflow and finite shell evaluation still require the separately prepared calculator's admission.

5. **Contact-safe vector kernel.** For V=phi*r, `|V|≤min(kappa*r,r^-2)≤kappa^(2/3)`. Its Jacobian eigenvalues are phi (twice) and `3*kappa*exp(−zeta)−2*phi`. Triangle bounds give at most 5*kappa; alternatively `(3*zeta*exp(−zeta)+2*(1−exp(−zeta)))/r³≤(3/e+2)/r³<4/r³`. Hence the stated V0 and V1 are conservative and include contact.

6. **Field error propagation.** For retained pairs the omitted environment changes q by the total omitted f mass. The bounds `|eta'|≤1`, `|eta''|≤2` give the retained prefactor and retained derivative changes in the definition. Omitted pairs add both their direct radial derivative and the derivative of their environment. Combining those terms reproduces every term in the stated dh bound. A direct omitted-field sum alone would miss the retained-neighbor contributions; this definition includes them. The B sums can bound the sum of primary-water vector norms, so the later maximum-dual-dipole contraction is valid.

7. **Global stiffness bounds.** Write `eta=1−(1+Q)/(1+s²)`. Its first partial is `2s(1+Q)/(1+s²)²−2f_i/(1+s²)`, which satisfies the declared bound 2. Differentiation gives exactly the four-term Hessian in the definition. Each term has magnitude at most 2 using `Q≤s²`, `f_i+f_j≤2s`, and `s²/(1+s²)²≤1/4`; therefore 8 is valid. Combining `g≤2`, `|g_i|≤2`, `|g_ij|≤4`, `eta≤1` yields `Lsigma1=2|c|+6|d|` and `Lsigma2=4|c|+28|d|`. The pinned domain gives `exp(sigma)≤4`. The beta Lipschitz/Hessian constants and the retained/omitted directional derivative bound follow by the mean-value theorem on nonnegative f weights. This remains valid for signed c,d and extends to the convergent infinite sum.

8. **Series remainder.** The native phi polynomial ends with `−kappa*zeta^5/720`; the first omitted alternating term is `kappa*zeta^6/7!`. Its derivative polynomial ends with `−kappa²*r²*zeta^4/48`; the first omitted term is `kappa²*r²*zeta^5/280`. On the native `0≤zeta<1e−4` branch, successive term magnitudes decrease, so these are valid analytic remainder bounds. They do not include machine rounding. The derivative polynomial is the derivative of the retained phi polynomial, preserving its internal consistency.

## Limits retained by this review

The explicit missing-image gradient bound holds at common fixed u,v. It does not certify the force difference between two self-consistent solutions. The stated conditional identity

`u_infinity−u_N = D_infinity^-1 [(h_infinity−h_N)−(H_infinity−H_N)u_N]`

supports the given inverse-response bound. Native mutual, permanent-frame and cell derivative contractions still need their own propagation bounds. The current vacuum SPD check supplies none of those periodic bounds. A uniform sampling guarantee also requires the declared nondegenerate cell/population/damping domain; a fixed-snapshot proof alone is insufficient.

No numerical molecular work, coefficient selection, calculator execution, control-suite rerun, native call, build, source mutation, or fit/outcome inspection was performed. This is an independent symbolic/source review of the frozen definition only.
