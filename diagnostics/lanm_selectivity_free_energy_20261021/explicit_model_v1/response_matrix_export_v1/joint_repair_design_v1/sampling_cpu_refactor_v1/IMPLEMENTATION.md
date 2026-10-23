# Vacuum CPU refactor implementation

The isolated source tree changes only `induce.f` and `polar.f`. The current
fit/native binaries, original source and separate donor extension are untouched.

## Common fields

`dfield0` and `ufield0` preserve the original Ewald/list/vacuum dispatch order,
then add h to both channels or subtract H times the corresponding u/v exactly
once. The original static external field remains after h. All four solver
field-dispatch sites and the exporter/dense-SPD probes use the common routines.
The old additions inside dfield0a/ufield0a are removed. The new enabled guards
still reject Ewald/list/periodic use before these routes execute. This is no
periodic permission or qualification. Other utilities and alternative solvers
remain outside this analyze/testgrad qualification.

## O(Q) geometry and force contraction

Selected-pair storage replaces `lrdh(3,3,n,n)`. Original radial/contact formulas
and beta derivatives remain; prefix/suffix sums provide each other-water sum
without subtracting its own potentially dominant contribution. Pair products
for q are accumulated in one pass. This preserves the exact algebra while
allowing normal floating-point summation changes for nonzero coefficients.

With w_i=u_i+v_i, z_i=3 phi_i r_i·w_i and pref_i=a+b eta_i, each pair stores
T_i=b f_i z_i/(1+other_i)^2. A second pass uses
C_j=pref_j z_j+sum_(i!=j) T_i and
G_j=K/2[(u_M·v_M) grad_j beta−C_j grad f_j
−3 f_j pref_j(phi_j w_j+(phi'_j/r_j)r_j(r_j·w_j))].
The contact limit of phi'/r is zero. Add G_j to water O, −G_j to the metal,
and sym(r_j outer G_j) to the original virial. Native local-frame torque paths
and reduced/analysis energies remain byte-identical and use the new dual solve.

All original parser/domain/topology/mode guards and 28-site dense normalized
SPD check remain. Dense-SPD work is still deliberately nonscalable; this change
removes the coordinate-Jacobian bottleneck only. It does not enable bulk MD.

## Evidence and next executable stages

77 bounded synthetic/source checks passed, including unequal dual dipoles,
contact/far/pair/collective cases, dense contraction equality, fixed dual-energy
finite differences, translation/rotation/permutation and symmetric strain.
They make zero native calls and are not actual native qualification.

The fresh57 runner is the original full-panel runner with its single startup
shared-dictionary publication changed to deepcopy, matching the admitted v2
persistence repair. Its final publication already occurs under lock. There is
no old27-output import or recovery gate: all57 programs/636 internal finite-
difference energies must execute against the newly built binary. All original
nine rows,57 literal input triples, numerical gates and helper code are retained.
The prepared manifest references original immutable input files; runtime copies
those inputs into fresh output directories. New executable/build/source-review
fields remain unbound.

BUILD_LOGIN.py preserves the original build command AST, compiler, flags and
targets, with root's separate BUILD_CONTEXT.json changing only the execution
resource envelope to truthful login/2 workers/null job ID. It records actual
controller/child CPU time and no scheduler allocation. No build or numerical
qualification has run; root owns both after independent review.
