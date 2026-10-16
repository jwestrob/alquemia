# Bounded intact-peptide accommodation pilot

Declared 16 October 2026 after force profile12512 and before any optimization.
Execution requires the complete two-step derivative check in
LBT3_RESPONSE_DESIGN.md to pass its unchanged numerical criterion.

## Scientific question

Can physically admissible response under the same solvent-containing energy
reduce the large starting-structure dependence, while permitting different
metal-specific conformations? This tests local accommodation in an intact
scaffold, not binding free energies or equilibrium conformational populations.

Both actual sources7CCO/7CCN and both metals remain in the denominator. Preserve
all256 atoms, charge−2, proton/explicit-water inventory, conditional termini,
MACE physical La1/Dy6 and GFN2 original f-in-core effective singlet. Use exactly
the qualified nativeMACE plus standaloneGFN2(ALPBwater−vacuum) energy and its
analytic derivative from12500/12510; no additional force-field energy or spring.

## Common coordinates and fixed search rules

Use the seven names selected by12512 for every source and metal: target x/y/z,
Asp7 CA–CB and CB–CG, Asp3 CA–CB, Glu11 CA–CB. Retain whole donor groups and
attached hydrogens. The whole peptide stays in the energy even where it is fixed.
Unused backbone motions remain available for a separately diagnosed extension.

At each source, form the seven-column physical Jacobian J0 over all256 atoms.
Let D contain its column norms and G=(J0 D^-1)^T(J0 D^-1). Require rank7 using
the existing relative1e−10 cutoff. Define T=D^-1 G^-1/2 and q=T u, making J0 T
orthonormal in the physical-displacement metric. This is numerical coordinate
scaling, not a stiffness or energy correction. At each evaluated geometry the
gradient is exactly T^T J(q)^T g_cartesian.

- Each torsion stays within ±0.20radian of its actual source.
- Target displacement norm and maximum displacement of any of the256 atoms
  stay within0.20angstrom of that source.
- Each trial moves any atom by at most0.05angstrom from the accepted point.
- Reuse the existing L-BFGS two-loop direction, history8; reset to negative
  gradient if it is not a descent direction. Scale by the current Cartesian
  tangent and verify the actual finite displacement before energy evaluation.
- Armijo acceptance uses the complete own-metal energy, coefficient1e−4, with
  at most four trials at step multipliers1,1/2,1/4,1/8. No difference objective.
- Each search has at most12 accepted steps and24 new distinct composite
  evaluations. These bounds define one finite experiment, not a project budget.
- Stop earlier for no admissible decreasing trial, an active-span projected
  maximum atom force≤0.1kcal/mol/angstrom and Cartesian RMS≤0.03, or a backend
  failure. No SCF restart, enlarged numerical tolerance or automatic rescue.

Re-use the source forces/energies exactly. Geometry-only rejected trials make
zero molecular calls. Cache identities include coordinates, state, both model
recipes, preparation and protocol. Count all attempts and failures.

## Physical admission and response diagnostics

Preserve existing source-relative bonds, stereochemistry, peptide basin,
protected-group geometry and stationary atoms. Check relative severe contacts
using current coordinates of ALL256 atoms, including the mobile metal. Keep
the original source as the comparison baseline; rebuilding it at each trial
would conceal accumulated deterioration. The existing peptide-only Geometry
helper alone is insufficient for a metal moving against fixed peptide atoms.

Physical rejection triggers backtracking, not relaxed admission. Backend,
electronic-state or gradient failure makes that search unavailable; any earlier
valid points remain archived separately. Reaching a search bound or finite limit
may yield a valid last accepted point, including the origin, but is never
advertised as a stationary or equilibrium minimum.

Report active-span projected forces, full forces, residual outside the active
span, and all64 current mode loads including omitted backbone modes. A small
active-span residual does not establish full-system stationarity. Record which
bounds are active and whether strain points into omitted motions.

## Staged execution and fair comparison

Run7CCO/La first on one warm allocated GPU with CPU solvent workers inside the
grant. If execution and physical admission succeed, run the other three with
identical settings, regardless of improvement magnitude or preference sign.
No scientific criterion selects a favorable source or metal.

Pool both original source geometries and all four admitted terminal search
geometries, deduplicating exact physical coordinates. This across-source pool is
permitted here only because the complete peptide composition and atom mapping
are identical. Evaluate both metals on every admitted candidate. Reuse exact
origins and own-metal endpoint evaluations; at most four opposite-metal endpoint
evaluations are newly required. Missing required cells leave the new pool result
unavailable while retaining the separately named original result.

The four searches plus cross-evaluation admit at most100 new MACE and200 new
GFN2 evaluations, excluding the previously declared derivative check. They are
planned maxima, not implied execution. Record actual latency, calls, allocation
and observed memory. Use current shared Slurm limits and completion wakes.

Report source-specific own-metal works, complete paired matrices, original and
new common-pool choices, change in paired source spread, physical checks and
cost. No probability, Kd, label threshold, individual Mex-site affinity or new
experimental validation follows. PQQ, original failures and all old LanM records
remain unchanged. The larger LanM profiling goal remains active after this pilot.
