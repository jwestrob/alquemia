# Exact CPU dual-field map and isolated prototype

The completed trimer comparison supports testing separate metal-pair driving and mutual-response kernels. All eight electronic E3 values are underestimated by 4.617–8.184 kcal/mol; collective dipole magnitudes are 58.8–85.6% too small, with all four La/Dy magnitude orderings reversed. These are consumed, provisional-reference diagnostics, not unique attribution to polarization or charge transfer. See the unchanged `TRIMER_RESULT_REPORT.md`.

The concrete prototype is in `../metal_pair_response_v1/`, with a separate source copy and five changed Fortran files. No original engine source, parameter or result is changed. This map supersedes only the schematic single-field functional in the earlier options report; that report remains preserved.

## One stationary energy and its derivatives

Use native `uind = u`, `uinp = v`, `udir = alpha E_d`, `udirp = alpha E_p`, `K = electric/dielec`, and `A = alpha_inverse − T_M`. The two permanent driving fields use the same selected direct kernel but retain their distinct native `dscale` and `pscale` exclusions. The mutual tensor uses `uscale`.

```
L(u,v;R) = K/2 [u^T A v − u^T E_p − v^T E_d]
A u = E_d;  A v = E_p
U* = −K/2 u^T E_p
dU*/dR = K/2 [u^T A' v − u^T E_p' − v^T E_d']
```

This dual functional is stationary and is a saddle in independent u and v. Positive definiteness is required of the symmetric physical response matrix A; it is not a joint-minimum claim for L. With fixed intrinsic polarizabilities, `A' = −T_M'`. Forces are minus this coordinate derivative. Both the permanent-field derivative and the mutual-tensor derivative must change consistently. No derivative of the converged u/v solutions is needed at stationarity; incomplete solves remain a numerical error.

## Exact existing routes

Line numbers refer to the unchanged source under `workspaces/.../explicit_model_v1/tinker_source/source/`; baseline and modified hashes are in the prototype source manifest.

| Native routine / lines | Actual role | Prototype requirement |
| --- | --- | --- |
| `polar.f` | `udir/udirp`, `uind/uinp`, `pdamp`, compact `jpolar`, pair `thlval/thdval` | Preserve all original arrays; add separate symmetric opt-in pair metadata. |
| `kpolar.f:381–519` | Reads original `POLPAIR`, resolves type-pair `thlval`; `pdamp=alpha^(1/6)` | Resolve gamma0 first; do not replace it globally or reinterpret `thdval`. |
| `induce.f::induce0a:238–248,402–405` | Builds the two direct dipoles, solves the same A against two fields | Preserve both PCG channels, guesses, stopping rule and preconditioner. |
| `induce.f::dfield0a:866,922–925` | `damptholed` orders 3/5/7; field versus fieldp apply dscale versus pscale | Selected D kernel applies to both fields. |
| `induce.f::ufield0a:1394` | `dampthole` orders 3/5; identical tensor acts on uind and uinp | Selected M kernel; preserve symmetric pair interaction and uscale. |
| `induce.f::uscale0a:7269` | PCG preconditioner calls `dampthole` | The same M lookup automatically reaches this path. |
| `epolar.f::epolar0e:1837–1852` | Energy `−K/2 sum(uind*udirp/alpha)` | Algebra unchanged; receives the modified converged u and E_p. |
| `epolar3.f::epolar3a:267` | Analyze uses pairwise energy and `damptholed` | Same D lookup as the driving field; compare with the single-loop energy. |
| `epolar1.f::epolar1a:438–527` | Direct factors and `rc3/5/7` derivative arrays | Use D factors and derivatives together. |
| `epolar1.f:553–620,849–869` | Permanent-frame torque accumulators and direct gradients | Keep D: dscale contracts with uinp; pscale contracts with uind. |
| `epolar1.f:887–951` | Mutual reset and mixed uind/uinp derivative contractions | Reset to M for **either** nonzero increment, including direct-only. |
| `epolar1.f:1212–1228,2417` | Pair-force virial and permanent multipole torque redistribution | Retain original accumulation; the changed D tensors reach these terms. |
| `epolar1.f::epolar1e:8445` | Gradient-side single-loop polarization energy | Same stationary energy as epolar0e. |

The stock alternate-direct reset uses `min(pti,thole(k))`, not the resolved pair `thlval`. The prototype never enables that global alternate mode; selected pairs instead use the same explicit M lookup in solve, preconditioner and gradient. Global alternate Thole also changes the exponent to 1.5 and has unsafe unmatched-pair semantics for this purpose; it is not the proposed repair.

## Kernel and shared pair selector

The keyword is `POLAR-PAIR-RESPONSE typeA typeB delta_direct delta_mutual`. Only unordered type pairs `{500,36}` (La/water O) and `{501,36}` (Dy/water O) are admitted in this prototype. This is an atom-type selector, not a class selector. It changes neither metal–H nor any water–water pair, including those inside a metal-containing system.

For the original resolved cubic pair coefficient gamma0:

```
gamma_D = gamma0 * (1 + delta_direct)
gamma_M = gamma0 * (1 + delta_mutual)
t = gamma_role * [r / (pdamp_i pdamp_k)]^3
s3 = 1 − exp(−t)
s5 = 1 − (1+t) exp(−t)
s7 = 1 − (1+t+0.6t²) exp(−t)
```

All deltas must be finite and greater than −1; resulting coefficients must be positive and finite. Duplicate/reversed-duplicate records, other type pairs and malformed records fail. Missing types are allowed for literal fragments. No numerical candidate increment is selected here; Planck owns the prospective finite qualification profiles.

Absent and all-zero records bypass new mode guards and retain literal original kernel arithmetic. The symmetric `pairresponse(i,k,role,gamma,selected)` is public for a small linked selector audit. `role=1` is D, `role=0` is M. On unmatched pairs it returns `selected=false`; original `dampthole`, `damptholed` and `epolar1a` branches then execute. `dampresponse` and `derivresponse` are the selected-pair shared factor and derivative helpers.

The derivative helper preserves native `rcN = grad(sN) * rrN`. In particular, the cubic recurrences are `rc3 = rvec * t exp(−t) rr5`, `rc5 = rc3 * 3t/r²`, `rc7 = rc5 * (−1+3t)/r²`. D and M factors are recomputed separately before their contractions. No independently chosen order-3/5/7 factors are introduced.

## Why a metal–O-only mutual change can alter collective response

Partitioning the response matrix into the metal and water sites gives the water Schur complement

```
(A_W − T_WM alpha_M T_MW) u_W = E_W + T_WM alpha_M E_M.
```

The metal-mediated cross-water coupling has rank at most three. Both frozen metal polarizabilities are nonzero (La 1.6597 and Dy 0.7280 Å³), so this is a real collective route even while A_W's original water–water kernel stays fixed. It does not guarantee the required sign or size, and the alpha ratio alone predicts neither.

The admitted W8 metal induced vectors cancel to printed zero by symmetry. W9 metal induced z components are La +0.374266 and Dy +0.198545 Debye. A symmetry-preserving mutual-only change cannot create a W8 central-metal vector in the exact symmetric limit; W8 B can nevertheless change through the subtracted isolated pairs. W9 has a nonzero full-shell route. These limits prevent claiming that this one extension must repair every full-shell discrepancy.

## Executable scope and qualification boundary

The isolated prototype supports vacuum, nonperiodic, no-list, original-Thole MUTUAL PCG only. For a nonzero declaration, its induction-entry guard rejects periodic/replica/polymer, Ewald, multipole/preconditioner lists, solvent, global alternate Thole, CHGPEN, EXCHPOL, charge flux, external fields, prediction and finite multipole cutoffs. Current positive finite site alphas are checked. Unsupported paths are not silently approximated. The guard runs after mechanic initialization.

Qualification must compare original/absent/zero behavior, audit unchanged pair factors inside metal-containing fixtures, compare analysis and both energy/gradient routes, verify analytic Cartesian and rigid-water torque derivatives, inspect both response channels and A stability, and retain exact native input/site/induction provenance. The finite call inventory and test increments live in `../metal_pair_response_v1/QUALIFICATION_PLAN.md`, owned independently. No build, native call, coefficient fit, GPU/periodic qualification or production change is claimed by this source preparation.

Native moment reporting uses uind, not the dual average. For general exclusion patterns with E_d != E_p, an electric-field derivative of this energy involves both channels. The current fixed-water tests must explicitly retain/check this distinction; protein electrical-response qualification cannot be inferred from the water fixture result.
