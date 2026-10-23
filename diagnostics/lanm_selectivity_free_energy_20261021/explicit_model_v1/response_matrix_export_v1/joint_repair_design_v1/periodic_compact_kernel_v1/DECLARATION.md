# Compact periodic response kernel prototype

## Frozen question and scope

Can a C² compact numerical approximation to the existing noncompact environment weight have a conservative minimum-image implementation and a uniform, outward scalar approximation bound? This is a mathematical prototype, not a native-engine patch or a fitted physical cutoff. The frozen 27-image plan and all fit/native sources remain unchanged.

Use r in Å, P(x)=1+x+x²/3, C=9 exp(10)/1849, and F(r)=C P(2r)² exp(−4r). Define t=(r−20 Å)/(1 Å) and S=1 for r≤20 Å, S=1−10t³+15t⁴−6t⁵ for 20<r<21 Å, S=0 for r≥21 Å. Replace F consistently by Ftilde=FS in every environment sum and h/H derivative. Original contact-regularized V, alpha and damping definitions are unchanged. No actual coefficients, dipoles, selected vectors, geometry files or response solves enter this task.

The uniform domain is orthorhombic cells with every side at least 44 Å and at most 3980 primary waters per metal environment. There is no upper cell-size bound. Positions are arbitrary, including contact. Minimum-image representatives have each component within its half cell side. Support 21 Å is strictly below half the minimum side, so at most one image per primary water contributes. Minimum image therefore equals the full compact lattice sum. Image-choice discontinuities lie outside support, where the kernel and its first two derivatives vanish.

## Bounds fixed before calculation

F is nonnegative and decreasing; |F′|≤4 F /Å. The product rF(r) decreases for r≥20 Å. For S, S′=−30t²(1−t)²/Å, hence |S′|≤15/(8 Å). S,S′,S″ match both constant branches.

Put a=44 Å, N=3980, R=20 Å, B=4+15/8=47/8 Å⁻¹. A nonnearest image shell k=||n||∞≥1 contains 24k²+2 image indices and has distance at least l_k=a(k−1/2). This follows componentwise from nearest representatives; it does not substitute D=a/2 into the old vector-radius assumption. Define

```
M0 = N sum_{k≥1} (24k²+2) F(l_k)
M1 = N sum_{k≥1} (24k²+2) l_k F(l_k)
T0 = N F(R) + M0
T1 = B N F(R) + 4 M0
V0 = N R F(R) + M1
V1 = B N R F(R) + 4 M1.
```

These uniformly bound, respectively, sum|F−Ftilde|, sum|F′−Ftilde′|, sum r|F−Ftilde| and sum r|F′−Ftilde′| over all image-labelled waters. All nonnearest compact terms vanish. Nearest-image loss starts only at R. For the weighted shell, monotonic rF gives the sharper l_k F(l_k), not an assumed upper radius. The proposed a k F(l_k) would also be conservative because l_k≤a k. V1 bounds weight error under unit-operator-norm homogeneous strain; translation/coordinate directional factors must be applied separately.

Evaluate the infinite sums as finite positive geometric-moment expressions: j=k−1≥0, l_k=22+44j, rho=exp(−176), and factor exp(−88). The polynomial (26+48j+24j²) P(44+88j)² (22+44j)^p has nonnegative coefficients and degree 6+p for p=0,1. Reuse the pinned Decimal outward arithmetic and positive Stirling geometric moments from TAIL.py. Use precision80, normal exponent domain ±999999, and fail on subnormal/underflow/overflow/clamping. No truncated lattice series or physical coordinate sampling is used.

## Finite implementation and controls

Implement a pure radial/environment prototype and one scalar-bound calculation. Tiny controls cover exact rational switch endpoint/positivity/slope identities, radial monotonic inequalities, finite-difference derivative and strain product rules including S′F, compact environment h/H chain derivatives on artificial inputs only, and nearest-image/full-compact-sum identity on artificial orthogonal cells. Finite-difference controls use fixed step 1e−5 Å and relative allowance 2e−7 (absolute floor 1e−12 for ordinary-sized artificial quantities); transition-kernel tests use relative scaling without that floor. These are numerical implementation checks, not physical tolerances. No optimizer, dependency installation, molecular call, native program, real-box evaluation or fit-outcome read. At most one local thread; retain source/config hashes, all results and measured wall/CPU time.

## Interpretation limits

The output intervals enclose the analytic upper-bound expressions; the physical errors lie between zero and their upper endpoints. They are scalar weight/first-derivative/virial-weight bounds only. Changes in retained-water collective factors must still be propagated through h/H. Full self-consistent energy/force/stress error additionally requires response stability, dipole and native derivative bounds. This prototype neither supplies that admission nor qualifies periodic PME, sampling, or physical transfer.
