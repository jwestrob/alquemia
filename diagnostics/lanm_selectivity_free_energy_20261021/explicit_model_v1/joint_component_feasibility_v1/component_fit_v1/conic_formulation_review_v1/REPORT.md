# Exact conic form of the frozen linearized component proposal

## Finding

**Every objective and proposal constraint in frozen CORE.py has an equivalent convex second-order-cone form.** It needs 13 decision variables: the same 12 normalized parameter increments and one nonnegative auxiliary variable t=sqrt(q). There is no nonconvex physical constraint hidden in this particular linearized subproblem. This is a source/mathematical remapping by the numerical-core author, not an independently executed solver comparison.

The current SLSQP campaign, its status rules, coefficients, budgets and native execution remain frozen. No solver was run or installed, no new proposal was formed, and no further campaign outcomes were read. `REVIEW.json` pins the inspected source and declaration.

## 1. Identical linear model

Let x_current be the physical 12-vector, l/u its original domain bounds, W=diag(u−l), and z=W^-1(x_current−l). The proposed physical parameters are

\[
 x=x_{current}+Wd.
\]

The domain and trust constraints are exactly

\[
 \max(-0.025,-z_j)\le d_j\le\min(0.025,1-z_j),\quad j=1,\ldots,12.
\]

Use the existing `linearize()` output without recomputing a different model. Each residual vector is affine, `v_b(d)=r_b+J_b d`, with J formed from the same ordered 24 probes and their **actual normalized spans**, including asymmetric spans near the box. The original 531 residual ordering, physical energy/force rows, and P8/B8/P9/B9 cycles stay unchanged. No Jacobian at the selected candidate is inferred from a different origin.

Define the Lorentz cone `Q_(n+1)={(s,y): ||y||_2<=s}`. The following is a complete standard SOCP with an affine objective and affine maps into Lorentz cones or the nonnegative orthant.

## 2. Objective and all constraint families

| Frozen family | Identical conic form | Count |
|---|---|---:|
| P fractional squared-error epigraph | `(t, (r_P+J_P d)/sqrt(D_P)) in Q_3` | 1 |
| B fractional squared-error epigraph | `(t, (r_B+J_B d)/sqrt(D_B)) in Q_3` | 1 |
| Original old_EF/pair_mu/C_mu/chi ceilings | `(1, (r_b+J_b d)/sqrt(C_b/4)) in Q_(n_b+1)` | 4 |
| Category-balanced body-force RMS | `(1, diag(sqrt(w_g))*(r_g+J_g d)/2) in Q_(n_g+1)` | 15 |
| Physical energy error | `1−e_j(d)>=0`, `1+e_j(d)>=0` | 20 sides |
| Individual body-force component | `1−v_gk(d)/5>=0`, `1+v_gk(d)/5>=0` | 810 sides |
| Original domain plus trust | `d_j−lo_j>=0`, `hi_j−d_j>=0` | 24 sides |
| Auxiliary interval | `t>=0`, `sqrt(q_current)−t>=0` | 2 sides |

Minimize **t**. Here D_P and D_B are the fixed origin24 J_P/J_B denominators, while C_b are the four fixed original block-score ceilings. They do not become current-iterate ceilings.

The component residuals passed by `linearize()` are `(e_P8,e_P9)/sqrt(2)` and `(e_B8,e_B9)/sqrt(2)`, so their squared norms are precisely J_P and J_B. No additional factor two or error weighting is needed. The four original slices have sizes 413,102,12,4. Their squared norms are one quarter of the stored block scores; hence the denominator is **C_b/4**, or equivalently the radius sqrt(C_b)/2. Using sqrt(C_b) would relax those constraints by a factor four in score.

The physical energy gates comprise four balanced exchanges, four individual radial works and two paired radial works. The latter two are separate acceptance gates even though they are derived from the individual works and absent as separate rows from the eight energy residuals in the old objective. They must remain in the conic model.

For each of La, Dy and paired, the force geometries contain n=1,1,1,8,9 waters. The body vector has `3+6n` entries. Its first three weights are 1/6 and the other 6n weights are 1/(12n), summing to one. Consequently `||diag(sqrt(w))*v||<=2` is the existing RMS limit, while **every raw component** must also satisfy `|v_k|<=5`. The vector includes the original net-force and torque/1Å convention; it is not an unweighted Cartesian RMS. Retain the literal weights and rows from `linearize()`.

There are 405 raw body components: `3*(3*9+51+57)`. Thus the existing 851 scalar inequalities decompose as 6 block-norm inequalities +15 force-RMS inequalities +20 energy sides +810 component sides. The SOCP replaces the 21 quadratic inequalities with 21 Lorentz cones and retains the 830 linear physical sides. Adding variable bounds gives **21 Lorentz cones and 856 linear sides**. The t lower bound is redundant with nondegenerate component cones but may remain explicitly present. No auxiliary residual variables are necessary.

## 3. Exact equivalence and existence

On q>=0, t=sqrt(q) is a bijection with inverse q=t². For each positive component denominator,

\[
 \|r+Jd\|^2/D\le q
 \iff \|(r+Jd)/\sqrt D\|\le t.
\]

All other constraints are independent of q; `0<=q<=q_current` becomes the stated t interval. Since squaring is strictly increasing on nonnegative t, minimizing t and minimizing q have exactly the same optimal d set, with optimal objectives related by q*=t*². The possible nonuniqueness of d remains; this introduces no minimum-step or secondary objective.

The original squared-norm formulation is itself convex: its nonlinear inequalities are convex quadratic sublevel/epigraph constraints. The SOCP is a different numerical representation of the same subproblem, not a repair to its physics or objective. With consistent exact residual/score arithmetic, d=0 and q=q_current are feasible because the current incumbent is qualified. The bounded feasible set is compact, so an optimum exists. A solver failure therefore cannot establish that this closed proposal feasible set is empty. Lack of a *strictly improving* point is a different question.

Stored scores and stored residuals can have tiny arithmetic discrepancies; CORE checks their agreement and permits a small discrepancy in its zero-step startup check. A future numerical representation must preserve those literal authorities and the final exact-threshold certificate, rather than clipping a residual or silently expanding a ceiling to manufacture feasibility. The exact equivalence is between the same mathematical constraints; separate floating-point evaluation orders need independent checking.

## 4. Preserve fraction, status, and current/tie semantics

The optimizer variable change does **not** authorize changing the existing trial policy. Given a returned conic point `(d_ret,t_ret)`, convert its auxiliary return to `q_ret=t_ret²`. For the four original fractions f=1,.5,.25,.125 use

\[
 d_f=f d_{ret},\qquad
 q_f=(1-f)q_{current}+f t_{ret}^2,\qquad
 t_f=\sqrt{q_f}.
\]

Linear interpolation of t would be a different, generally tighter epigraph path. The frozen q interpolation must be retained. It is feasible in exact arithmetic when both endpoints are feasible because the original squared-norm epigraph is convex. No new fraction, clipping, or line search follows from the conic representation.

Reapply the existing scalar `math.fsum` certificate to the returned **q** representation, all 851 residual inequalities, literal physical parameter/trust bounds, all25 predicted gates and four ceilings. A reported successful optimizer status is insufficient. Failed optimizer statuses remain failed; the current policy accepts their fractions only when all residual margins are strictly positive. A point on a legitimate cone boundary can be feasible yet fail that stricter failed-status witness rule. This distinction is policy, not physical infeasibility.

Only complete, native-admitted eligible trial slots participate in `select()`. Its actual objective remains `q=max(J_P/D_P,J_B/D_B)`, and its tolerance remains `1e−10*max(1,abs(q_current))`. Do not substitute that same epsilon in t units. Scan the original fraction order; current wins ties, then the earlier fraction. Derivative probes remain excluded. Actual25 and fixed-four-ceiling checks are unchanged. Conic objective equivalence does not imply identical finite solver iterates or identical native candidates.

## 5. Numerical obstacles and interpretation

- Affine residual maps may be rank deficient, and paired energy/force rows create dependence. Some fixed ceilings may be active at the origin. There is no source guarantee of strict feasibility or a unique solution. These are numerical/duality considerations, not nonconvexity.
- Keep the normalized cone forms above to avoid unnecessary unit and scale disparity, but independently certify with the original scalar quantities. Rounded square roots, solver feasibility tolerances and objective transformations must not relax physical thresholds.
- A future conic solver would require a separately declared, pinned numerical-method qualification; the current SciPy/SLSQP source supplies no conic backend or conic dual certificate. This review installs nothing and authorizes no alternative solve or campaign continuation.
- A primal feasibility certificate proves a local secant witness. An independently checked conic dual lower bound and primal/dual gap could establish near-optimality of this **fixed linearized trust-box problem**. A solver status alone cannot. A certified lack of local predicted improvement still would not establish global inadequacy of the nonlinear Hamiltonian, parameter domain, or physical mechanism.
- Numerical method reliability must be separated from response-family adequacy before failed SLSQP statuses are used in scientific interpretation. The source establishes a concrete equivalent convex formulation that can be qualified prospectively; it does not claim such a repair has been performed or is guaranteed to succeed.
