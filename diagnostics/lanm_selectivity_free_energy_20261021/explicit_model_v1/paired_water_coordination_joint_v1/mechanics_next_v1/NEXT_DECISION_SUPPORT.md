# After the completed alpha diagnostic

**Decision: defer further H-Thole probing.** Weak H screening and the small fixed H-vdW contrast make it low priority. Root will weigh reference uncertainty before another mechanics parameter experiment. Do not start an alpha fit, add H parameters to a fitter, implement an adapter, or launch a new probe. The four-vector/40-program inventory below records the smallest plausible H-screening discrimination if this hypothesis is later reopened; it is not an active recommendation to execute.

## What the actual alpha result changes

17695 completed all 40 physical cells/80 programs, with independent raw review passing. None of the eight probes improves the aggregate objective or passes all controls. The W9-minus-W8 paired energy error spans only −4.371455 to −4.305389 kcal/mol around −4.338284. The 1% and 5% alpha columns differ in normed response by only 0.0715% (La) and 0.0630% (Dy), so small-step inconsistency does not explain this weak energy response.

The rank-eight unconstrained linear projection lowers loss to 0.141490 but requires alpha La −0.176908 and Dy −4.361643 Å³, an invalid negative Dy Thole coefficient, and changes outside the old trust region. Even this extrapolated projection passes only 23/25 controls. It supports neither a physical candidate nor another unchanged optimizer continuation.

## Source evidence: H screening is distinct but may be largely redundant

Native `kpolar.f` assigns `pdamp=alpha^(1/6)`; `damping.f` uses

`t_Mj = a_Mj r_Mj³ / sqrt(alpha_M alpha_j)`.

The fitted O coefficient is `POLPAIR 36 500/501`. H type 37 instead retains native minimum coefficients La 0.25 and Dy 0.139. At fixed alpha, an H-specific `POLPAIR 37 500/501` changes neither O screening nor intrinsic susceptibility. It is therefore a different physical edit from the alpha probe.

It is **not guaranteed to provide an independent residual column**. Locally, in fractional/log coordinates,

`D_alpha = D_susceptibility − ½ D_aO − ½ D_aH`.

Here `D_susceptibility` means varying the intrinsic response with all screening exponents held fixed. The O term is already in J6; alpha is measured. If the explicit susceptibility response is negligible or already represented, the H column repeats their span. Symmetry can suppress the metal induced dipole in W8; it must not be assumed suppressed in the asymmetric W9/dimer fixtures. The saved native outputs do not contain the induced dipoles needed to isolate that contribution directly.

Quantitative geometry/parameter inspection also predicts a weak H-screening perturbation. Both W8 and W9 have M–H distances 3.177364 Å. With water-H alpha 0.496 Å³, their exponents are La 8.838638 and Dy 7.420092. Native `dmp5=1−(1+t)exp(−t)` is already 0.998573 and 0.994956. A ±5% H coefficient change changes this factor by at most 0.000693 and 0.001944, respectively. These are coefficient calculations, not new energy/force evaluations, and do not bound the self-consistent total response.

Nor does the saved vdW decomposition strongly motivate an H-vdW detour. Since the O-pair epsilon contribution is exactly linear, the existing selected-origin epsilon probes give `E_MO=epsilon·dE_vdW/depsilon`. Subtracting this from saved vdW totals and taking the matched W9−W8 Dy−La contrast yields O −5.178226 and fixed H −0.111866 kcal/mol. Water–water terms cancel in the paired contrast. Printed-energy precision limits the last digits; this decomposition does not identify a uniquely correct QM term.

## Deferred finite test, only if subsequently declared

At the unchanged 17676 six-vector and original metal alphas, append only the explicit H pair coefficient: La `POLPAIR 37 500` at 0.2375/0.2625 and Dy `POLPAIR 37 501` at 0.13205/0.14595. Four fixed probes, each evaluating its five consumed fixtures and reusing the opposite five selected cells: **20 new physical cells, 40 analyze plus analytic-testgrad programs; no Cartesian FD**. Preserve original keys, O coefficients, charges, multipoles, vdW, water model, geometries and references. No new donor or reserved outcomes.

Report every unchanged 25-control score and the two finite ±5% H-response columns. Compare their residualized span against J6 and against J6 plus the measured 1% alpha columns, retaining the measured 5% alpha matrix as a saved scale check. Use the declared relative singular cutoff 1e−6 and report the absolute cutoffs. Report both raw column norms and residualized norms, including the W9 energy row and rigid-force fixture blocks. A new rank alone is insufficient; assess its magnitude and alignment with the remaining error. Center-midpoint curvature must stay explicit because H damping is nonlinear; a single ±5% secant does not qualify derivative precision.

If the H response is redundant or only materially useful outside the sampled range, close this local H-screening hypothesis. If it adds a sizeable, locally reliable response, that warrants a separately declared representability/physical-parameter question, not an automatic fit or release. No new programs are authorized or implemented by this note.

## Evidence pins

- Alpha analysis: `metal_alpha_response_v1/analysis_v1/RESULT.json`, SHA `0c0df0b2880d264b69aeab9f9f3162688045c489a31c05667765b43956e00a61`.
- Alpha execution: SHA `57db6e96a1c13a18ca3b80510d1d603be1123c2c1a7bce0f876078b710eedc41`; independent raw review: `b7ad46ea8d126c2799bf99e958b00447b50276a82a06495b6673ee6b4fcb3b5f`.
- Selected/J6/epsilon derivative authority remains terminal 17676 result `5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832`; compact consumed data `e2e58963297984305856210f93cd664f08903485d9413f5f7a3d10ba1a616250`.
- Native source pins are in the frozen alpha configuration: kpolar `01c5cba4…`, damping `8e7b99d2…`, induce `759e5447…`; buffered14-7 epsilon linearity is explicit in `tinker_source/source/ehal.f` lines 205–210. Original proposal is preserved.

All work for this note was saved arithmetic/source inspection at one CPU thread, with zero molecular calls or source changes. Provisional ECP-PBE0-D4 reference uncertainty, prior W9 geometry interpretation and the 17670 review-loss notice remain. Cartesian, fresh validation and sampling gates stay held.
