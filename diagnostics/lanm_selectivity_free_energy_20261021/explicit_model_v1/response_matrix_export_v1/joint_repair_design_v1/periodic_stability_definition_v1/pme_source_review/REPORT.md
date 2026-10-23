# Reciprocal PME: conditional source proof and an upper bound

## Conclusion and scope

The inspected CPU implementation has the required reciprocal **minus sign**, transpose scatter/gather weights, and unnormalized forward/backward FFT pair. In exact arithmetic its reciprocal contribution is positive semidefinite after changing the sign from field to stiffness, provided the Fourier multiplier is finite, nonnegative, and conjugate-even, and the scatter visits precisely the same wrapped stencil as the gather.

The evenness condition is an essential refinement of the hypothesis in the preserved parent `ROOT_DESIGN.json`: native PME packs the two induced-dipole channels into one complex grid. Nonnegative multipliers alone do not prove that these channels remain independent. The initial proof domain should be **orthorhombic periodic boxes with conducting (`TINFOIL`) boundary**, exact current PME initialization, and checked stencil coverage and multipliers. No periodic implementation, configuration, molecular calculation, or numerical stability certificate is released by this note.

`REVIEW.json` pins the actual source files and parent declaration. Source is the already qualified CPU refactor, manifest `9ccee436f5611dde819c4b59bc9513b6651c4e00af3a40a0dcaaca1f047d503c`. Its enabled local-response route still rejects periodic execution. Existing vacuum qualification and PCG convergence provide no periodic proof.

## 1. Actual operator and adjoint

Let `p=bsorder`, `M=nfft1*nfft2*nfft3`, and let `u` contain the Cartesian dipoles at active sites. Define the fractional conversion

\[
 A_{ba}=nfft_b\,recip(a,b).
\]

At a wrapped grid point in site i's p-cubed stencil, define

\[
 d_{gi}=(\theta'_1\theta_2\theta_3,
          \theta_1\theta'_2\theta_3,
          \theta_1\theta_2\theta'_3),\qquad
 C_{g,ia}=\sum_b d_{gi,b}A_{ba}.
\]

Repeated indices, if admitted, must be summed at the same grid point. `umutual1` applies A to each dipole (`induce.f:3930–3945`); `grid_uind` adds exactly these three derivative products (`pmestuf.f:957–978`). Thus the real grid receives Cu and the imaginary grid receives Cv. No charge factor, electrostatic energy conversion constant, or extra polarizability occurs here.

`fphi_uind` accumulates the same derivative products in its first three derivative components (`pmestuf.f:1549–1609`, then `1635–1650`). `umutual1` converts those components with A-transpose (`induce.f:3986–4003`) and **subtracts** the result from each field (`4008–4013`). Therefore reciprocal field is −Ru, −Rv when the grid convolution is real, with

\[
 R=C^T F^*QFC,\qquad Q=\operatorname{diag}(qfac).
\]

The stencil alignment is literal: `getchunk` gives `grdoff−nlpts=2`, so the unwrapped scatter range is `igrid+2 ... igrid+p+1`, identical to the gather loop. `adjust` clips this range to a chunk and updates the spline offset; the gather wraps the same global indices. Chunks partition the grid, while `pmetable` marks the central and neighboring chunks. This proves the correspondence when each support reaches at most the immediate neighboring chunk and each wrapped support point is covered once. A sufficient structural check is p≤each grid dimension and each chunk width≥max(nlpts,nrpts); the native default p=5 and default allowed grids meet this check. Arbitrary spline orders/chunk settings must not be assumed covered. A future certificate should use the actual covered stencil, rather than a separately approximated deposition kernel.

## 2. FFT normalization and two-channel condition

`fft3d.f:40–47` creates FFTW plans with signs −1 and +1, and its forward/backward routines contain no division by M. The native FFTPACK route also has no normalization: `cfftf/cfftb` use the paired radix 2/3/4/5 butterflies with opposite imaginary signs (`fftpack.f:141–342,501–760,919–1178`). `kewald.f:52–60` restricts the selected dimensions to these factors. Accordingly, for the ideal DFT,

\[
 F^*F=MI,\qquad \|F\|_2=\sqrt M.
\]

The FFTW convention and legacy-Fortran dimension reversal are also explicitly documented by the upstream [DFT definition](https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html), [multidimensional definition](https://www.fftw.org/fftw3_doc/Multi_002ddimensional-Transforms.html), and [Fortran interface](https://www.fftw.org/fftw3_doc/Fortran_002dinterface-routines.html). These describe the API convention, not a measurement of a prospective compiled PME execution. `kewald.f:69–70` defaults to FFTPACK with one thread and FFTW with multiple threads; linking FFTW does not establish which route executes.

Let H=F*QF=X+iY. A real Q makes H complex Hermitian, so X is real symmetric and Y is real skew-symmetric. Native packing gives

\[
 H(Cu+iCv)=(XCu-YCv)+i(YCu+XCv).
\]

Hence **Q(k)=Q(−k)** is needed for Y=0 and independent applications of one real R to u and v. Neither u=v nor convergence of the separate PCG recurrences proves this condition. Under the condition,

\[
 x^TRx=\sum_k qfac_k\,|(FCx)_k|^2\ge0.
\]

## 3. Multiplier sign, masks, and box restriction

`udirect1`, `induce.f:2421–2472`, computes for each nonzero represented frequency

\[
 q_k=\frac{\exp[-(\pi/a_{Ewald})^2|h_k|^2]}
 {\pi V|h_k|^2\,bsmod_1(k_1)bsmod_2(k_2)bsmod_3(k_3)},
\]

or zero when its exponential cutoff is reached. For positive volume, nondegenerate reciprocal cell, positive finite spline moduli, and finite positive Ewald coefficient this is nonnegative. `dftmod` forms squared Fourier moduli, repairs small entries by neighboring nonnegative averages, and multiplies by a squared zeta factor (`pmestuf.f:187–234`). Those constructions motivate positivity, but division and finite arithmetic still require actual finite/positive checks; the source has no general proof that every user setting avoids zero denominators.

For periodic cells q(0)=0. The nonprismatic parity mask only replaces entries by zero, but is outside the proposed initial box domain. The separate nonperiodic multiplier and its nonzero origin are also outside this domain. The exponential mask is even in an orthorhombic box in exact arithmetic. The actual stored moduli and qfac should be checked, not inferred exactly symmetric from the formula.

The restriction to orthorhombic boxes matters. With even FFT dimensions the code represents a Nyquist component by −n/2. Its modular negative has the same represented Nyquist component. For a mixed frequency (−n/2,k,0), conjugation yields (−n/2,−k,0); in a nonorthogonal reciprocal metric the cross term changes sign, so |h|² and qfac need not agree. No qfac symmetrization appears in the inspected source. This is a source-level obstruction to the unrestricted real two-channel proof, even if a particular affected mode is small or masked to zero. No engine change is proposed here.

`ufield0c` adds positive self field `4*aewald³/(3*sqrt(pi))` and applies the separate cell-dipole subtraction only for `boundary='VACUUM'` (`induce.f:3843–3876`). Thus in the stated conducting domain the split is exactly

\[
 T=T_{real}+T_{self}-R.
\]

## 4. Computable upper bound without a dense response matrix

Let S=diag(sqrt(alpha)) on active Cartesian unknowns and B=CS. Zero-polarizability sites are excluded. The same active alpha convention used by the actual solver is required; native PCG's `polmin` must not be silently replaced by smaller nominal alpha.

Accumulate absolute row and column sums of B from the actual spline weights, including A and sqrt(alpha), and its squared Frobenius norm. Define

\[
 c_1=\max_j\sum_g|B_{gj}|,\quad
 c_\infty=\max_g\sum_j|B_{gj}|,\quad
 c_F=\sum_{g,j}|B_{gj}|^2.
\]

Then the normalized reciprocal stiffness obeys

\[
 0\preceq\bar R=B^TF^*QFB,\qquad
 \|\bar R\|_2\le U_{rec}
 =M\,\max_kq_k\,\min(c_1c_\infty,c_F).
\]

This takes O(N p³+M) arithmetic and O(M) temporary row storage, or streamed chunk storage; there is no global 3N-by-3N matrix or response solve. Grid aliases must be combined before cF is formed; for row/column sums, summing absolute individual contributions instead is a conservative bound. All bound arithmetic needs outward rounding or an explicitly justified error allowance. The M factor is mandatory because native inverse FFT is unnormalized.

This bound can be very loose: taking max(q) before accounting for cancellation in the derivative stencil loses the low-frequency gradient factor. It is a sufficient bound, not an expectation of practical admission. In particular it cannot repair a nonpositive real-space lower bound. The separate saved vacuum block-bound result already finds such inconclusive W9 bounds; this source proof does not overturn that result.

## 5. What remains before any certificate can be used

1. **Initialization identity.** `udirect1` builds qfac, moduli, splines and chunk tables; `umutual1` reuses them. Current `localresponsesetup` calls `localresponsespectrum` before `induce0a` calls `dfield0`. Replacing that spectrum routine in place would use uninitialized/stale PME data. A future implementation must initialize the exact current geometry/cell/grid first and invalidate it on changes, while retaining the current disabled periodic guard until separately qualified.
2. **Arithmetic error.** The proof above is exact arithmetic for the stated discrete PME operator. Stored conjugate pairs may differ from rounded moduli/exponentials, and FFT/deposition/interpolation are floating-point operations. For analysis only, define q_even(k)=(q(k)+q(−k))/2 and delta_q=max|q−q_even|. The multiplier-only perturbation on the packed real two-channel space is bounded by `M*delta_q*min(c1*cinfinity,cF)`. This is not an instruction to average qfac in the engine. FFT and scatter/gather errors remain additional terms; no rigorous numerical allowance has been established here. A loose equality tolerance without its propagated bound is not a proof.
3. **Integration.** With independently justified total error epsilon, root's sufficient bounds become `L_safe=L_real−epsilon`, `U_safe=U_real+U_rec+epsilon`, with the unchanged requirement `L_safe>1e−6*U_safe` and `U_safe>0`. Failing or unavailable bounds are inconclusive and must fail closed. No solver-convergence fallback establishes stability.
4. **Hamiltonian derivatives.** Reciprocal positivity alone does not qualify periodic energy, coordinate gradients, cell derivatives/virial, local h/H lattice tails, or a GPU implementation. Their exact consistency and finite-precision qualification remain separate requirements.

Source reads and mathematical derivation only: zero native calls, builds, fits, molecular numerical tests, or source edits. No fresh/held reference outcomes were read.
