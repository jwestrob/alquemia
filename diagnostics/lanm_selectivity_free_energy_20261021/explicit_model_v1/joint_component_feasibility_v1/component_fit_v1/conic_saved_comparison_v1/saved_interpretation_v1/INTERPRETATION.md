# Saved first-proposal conic interpretation

The comparison strongly disfavors a material **first-proposal optimization error** as the explanation for limited improvement under these constraints. It does not establish Hamiltonian inadequacy or qualify the solver at the declared full tolerance.

## Evidence and numerical status

The exact saved result is `actual_v1/RESULT.json`, SHA256 `341472cdf50e3a1a8c9479a38cb6e79acf9fde6a32189a0e551de368df6c4f6e`. [VALUES.json](VALUES.json) pins that result, its stored conic data/row layout/linearized problem, original SLSQP proposal, configuration and source, and records the scalar summaries below. This interpretation reads only those saved artifacts; it invokes no solver, native program, reference or candidate evaluation.

At immutable origin 24, Clarabel returns q = **0.8427637922817901**, versus SLSQP **0.8427637921792542**: a difference of **1.02536e−10**. The largest difference between their twelve normalized increments is **3.78510e−7**. Both represent about 15.724% reduction in worst normalized squared component error, or 8.198% in its RMS, in this linearized model. The original campaign's first SLSQP call reported success after 49 iterations. This comparison does not examine later failed SLSQP proposals or certify their directions.

Clarabel stops after 20 iterations as **AlmostSolved**. Its scaled primal residual is **2.26834e−9**, dual residual **6.45134e−14**, and primal–dual objective gap **1.63133e−8 in t = sqrt(q)**; primal feasibility and gap miss the frozen 1e−10 requirements. Saved unscaled equation residual maxima are **5.35502e−8** and **9.57234e−13**, respectively. These use different normalization from the solver's reported residuals. The positive nominal dual objective, squared, is 0.8427637623298927—about 3.0e−8 below the primal q—but is not a rigorously certified lower bound.

All saved dual-cone margins and solver slack-cone margins are positive in ordinary-double arithmetic. The affine slack computed from the returned parameters has a **−4.37094e−11** P-energy cone margin; the independently evaluated squared epigraph margin is **−8.02521e−11**. Consequently the full fraction remains ineligible. The three predeclared smaller fractions remain strictly eligible in the saved linearized certificate; none is newly evaluated here. The independent actual review `c59a4c5a256abc079380f320c502b3cc518fffce28b1e17070408323962c2285` admits arithmetic reproduction only. It does not turn AlmostSolved into Solved.

## What restricts this proposal

The four fixed score ceilings are development safeguards, not separate experimental accuracy tolerances. The table uses each original score divided by its fixed ceiling; the actual weighted residual slice squared norm is one quarter of that score, a factor that cancels in the ratio.

| Block | Fixed original-score ceiling | Predicted score / ceiling | Local dual sensitivity of q to a unit relative squared-ceiling relaxation |
|---|---:|---:|---:|
| Old E/F | 0.15340308137508737 | 0.9999999999974 | −2.96109 |
| Pair dipole | 0.37816202481323913 | 0.9999999999986 | −7.75253 |
| Collective dipole | 0.0052286591834414045 | 0.9864180744792 | approximately zero |
| Susceptibility | 0.0011670073007337461 | 0.9999999970726 | −0.0238063 |

E/F, pair dipole and susceptibility therefore bind numerically, while collective dipole retains 1.35819% squared-score slack. The tabulated sensitivity is −t times the SOC dual head at a unit ceiling; it agrees closely with the original SLSQP squared-constraint multipliers. These are local shadow prices, not predictions for a finite relaxation or attribution of physical error. Pair-dipole preservation has the largest such price.

Both P and B epigraphs bind. Their dual heads are 0.8114467 and 0.1885533 and sum to approximately one. The La epsilon lower trust face and Dy epsilon upper trust face also bind. No individual physical gate is close to equality: maximum energy error is 0.14224 versus 1 kcal/mol, maximum force RMS 1.90856 versus 2, and maximum force component 4.69114 versus 5 kcal/mol/Å. Their dual multipliers are correspondingly negligible. The saved raw dual objective contributions are retained in VALUES; they are algebraic contributions to −b·z, not physical component energies.

The matching objective, direction, active constraints and multipliers provide substantive numerical evidence that SLSQP already found essentially the same constrained first-step solution. The remaining numerical defect matters for exact admission, but does not plausibly explain a large missing improvement in this particular linearized subproblem. This conclusion covers neither nonlinear response away from the origin nor model capacity elsewhere in the parameter domain.

## Single next scientific decision

The discriminating contrast is one **separately preregistered constraint-only comparison** at this same origin and exact saved Jacobian: remove the four development score ceilings together, retain all 25 physical gates, the same parameter domain and normalized ±0.025 trust region, and minimize the same P/B minimax objective. Compare objective evidence and report the four resulting scores without choosing a new native candidate. This tests whether the preservation policy obstructs useful local improvement; removing all four avoids selecting a favorable subset after seeing multipliers.

A materially lower q with independently supported numerical feasibility would identify a local preservation tradeoff. Little change under a sufficiently supported optimum would instead identify a limitation of this local response span, trust region and physical-gate combination; it still would not prove missing physics globally. An unresolved solver result would leave the question unresolved. No new physical tolerance is supplied by this test.

**No such solve is released here.** The existing Solved-gated ceiling-removal adapter remains unrun, with its gate unchanged. Any future contrast requires a separate prospective root decision and numerical admission specification, not retrospective relabeling of this result or an automatic retry.

Original comparison cost: 3.14355 s wall / 3.04954 CPU-s, one thread, one solver call, zero native calls. This saved-summary calculation: 0.06506 s wall / 0.05967 CPU-s, one thread, zero solver/native/reference calls. Earlier source reads are not included in that summary timer.

### Separately declared follow-up preparation

Root subsequently authorized the source-only `ceiling_tradeoff_v1` and its saved baseline certificate. The exact returned physical vector is unchanged; recomputed rho is 0.8427637923620422 and q_aux=rho+1e−8 is 0.8427638023620423. All851 original scalar inequalities pass strictly (minimum 1.41642e−12); the raw failed epigraph certificate remains stored. Baseline RESULT SHA256 `3b3ad7bdebf4a0df7c5a7339d0bd23e553b6cf312f09bc86bbbcf00a2b244e48` costs 1.27590 s wall / 1.24334 CPU-s, zero solver/native calls.

The new prospective contract also retains each relaxed return's raw certificate before applying the same auxiliary-only construction at each unchanged physical fraction. It requires strict retained physical constraints independently of the solver status and compares actual rho, excluding both numerical cushions. Thus a status string cannot veto a genuinely feasible diagnostic witness, nor can epigraph bookkeeping manufacture gain. Original solver statuses and gates are preserved. No new solve has run; the original Solved-gated ablation remains untouched and unrun.
