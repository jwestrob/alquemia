# Four common donor coordinates: preparation and gated execution

Root authorized this bounded implementation on 22 September after the archived
force diagnostic. This task prepares/tests the adapter; it does **not** launch
optimization. Execution follows parent review and completion of the separate
four-context shared-pool pilot. No further user permission loop is implied.

## Scope and scientific rule

Use exactly 1H4I, 4MAE, PQQSEQ_83440678cbbd658047c9 and
PQQSEQ_07ab500e3df76b30d71c from the original30 sample0 proposal archive. Both PLM
labels remain unknown. Eight endpoint optimizations, one q0 start each. Keep
the original atoms, charges, protonation, source preparation, native OMOL model
and water/cofactor states. No new metal, PQQ, water or scaffold modes; no new
DFT, SCF rescue, entropy, calibration or automatic expansion.

Select four common angular primitives at the matched original q0 with the
already recorded alternating differential/individual, physical-heavy-normalized
Gram-Schmidt rule. Both metals use the exact same selected IDs. Optimize their
primitive q values through the existing nonlinear Kinematics mapping and actual
cap-aware Cartesian-force Jacobian, not the linearized orthogonal ranking basis.
Existing native origin replay tolerance remains 0.01 kcal/mol.

## Physical bounds and solver

The admissible **final candidate** displaces no physical source heavy atom by
more than 0.8 Angstrom from that endpoint's unchanged q0. Each selected angle
also lies in [-0.8,+0.8] radians; the angular box alone is insufficient for
proximal chi1/chi2. With delta_i(q)=x_i(q)-x_i(0), enforce the coupled inequalities

`c_i(q) = 0.8^2 - delta_i(q)·delta_i(q) >= 0`

and their analytic Jacobian `dc_i/dq_k = -2 delta_i·J_ki`.
This is the same explicit physical constraint used by
`coordination_preparation_context.py`, with the declared new extent and selected
columns. Use its SLSQP settings: maxiter200 and ftol1e-9 in eV, objective shifted
by the fixed original energy. Convert inherited kcal/mol mode gradients back to
eV for this optimizer, exactly once. There is no spring, penalty energy,
coordinate clipping, projected trial objective or alternate start.

Root explicitly resolved trial semantics before preparation: **SLSQP may evaluate
infeasible intermediate geometries; 0.8 Angstrom is not an all-evaluation domain.**
Angles remain bounded throughout. Retain each trial's actual extent, completed
model-call status and the count/maximum of infeasible requests. Existing finite
energy, source/cap/bond and severe-overlap safeguards remain. Final infeasibility
fails the endpoint without a rescue. Final tolerance is the existing 1e-7-Angstrom
numerical bound tolerance; fixed-source roundoff policy remains1e-12 Angstrom.

Require optimizer success, valid final geometry and no native energy increase
beyond the old1e-7-eV numerical allowance. Record final normalized selected/omitted
loads, physical-boundary atom count and the existing0.02-radian angular boundary
flag. These are diagnostics, not a claim of an unconstrained physical minimum.
Large residuals or boundary contact do not automatically authorize more modes.

SLSQP replaces the previous L-BFGS-B box optimizer because that implementation
cannot enforce coupled physical displacements. The older broad SLSQP pilot lowered
energies but worsened the tested separation and hit every0.20-Angstrom boundary.
This four-mode experiment is a new contained candidate-generation test, not a
claim that the previous negative result was fixed.

## Execution, output and comparison

Use the existing warm native GPU worker, oneH200/32CPU/200000MiB allocation,
immutable task manifest/cache keys and separate own lock/receipts. Preparation,
dry-run and collection make no model calls. Runtime requires the actual compatible
four-case shared-pool collection with every required pool cell available.

Export all eight candidate statuses, source geometry/map/state pins, actual MACE
receipts and residual/boundary traces. Missing candidates remain unavailable;
uncomputed GFN2/composite fields remain null. Parent's pool adapter evaluates both
metals on the same retained physical geometry set and retains the origin/prior
candidates. This proposal adapter does not choose a discriminatory state or launch
solvent scoring. Old/new descriptor/reference records remain distinct.
