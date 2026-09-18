# Frozen-density / AMOEBA-GK / MACE short-context pilot

Declared after all eight density-multipole numerical checks and all sixteen
native boundary preparations completed successfully. This is the first complete
score from those components. It is a new opt-in model, with no baseline/default
change, inherited zero or calibrated bands. All four representations and both
biological groups are already consumed development evidence.

## Fixed energy and scientific state

Protocol `vacuum_density_direct_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1`.
Use exactly the eight normalized vacuum native r2SCAN-3c endpoints, the paired
charge-closed AMOEBA boundary in density_gk_boundary_v1, the actual fine quantum
fields and fine spatial potential Hessians, and the qualified existing MACE-POLAR
medium short full/core readouts. No new DFT, charge fit, MACE, force, geometry,
protonation, water, assembly or cofactor change. Preserve all electronic/core
source hashes and raw energies; a prior *total* hybrid is never reused as a term.

For a physical cavity B and a frozen Q-proxy charge distribution Pq, define:

```
V_QE = sum_E [q_E phi_Q − mu_E·E_Q + Q_E:Hessian(phi_Q)]
G_QE = es_static(Pq,E;B) − es_static(0,E;B)
I(F) = −0.5*C*sum_E(mu_solv,d(F) · F_solv,p), C=electric/dielec

DeltaU = V_QE + G_QE + I(F_density+proxy_RF) − I(F_environment)
A_M = E_DFT,vac(core,M) + DeltaU_M + T_short(full,M) − T_short(core,M)
R = A_Ca − A_La
```

`es_static` is the actual native GK-plus-nonpolar component with all induced
dipoles zero. Same physical cavity/radii/nonpolar parameters make nonpolar and
environment-only terms cancel. The difference retains Q reaction self energy
and Q–environment reaction-field terms. Do **not** include native permanent
vacuum `em`: that contains proxy Q–Q and Q–E electrostatics, already replaced
by DFT and the actual density coupling. Report it only as an audit component.

The actual native four-field arrays give

```
F_Q_RF,d/p = (F_static(Q,E),solv,d/p − F_static(Q,E),vac,d/p)
            − (F_static(0,E),solv,d/p − F_static(0,E),vac,d/p)
F_total,vac,d/p  = F_environment,vac,d/p  + E_Q,density
F_total,solv,d/p = F_environment,solv,d/p + E_Q,density + F_Q_RF,d/p
```

Keep d/p channels distinct. Actual density fields in atomic units become native
e/Å² by division by bohr_to_A². The supplied-field native solver uses the
unchanged MM mutual/GK response operator and source-frozen mask. Quantum fields
are required on every permitted environment site; Q-site field entries are
inactive coordinates, not invented missing corrections. Permanent environment
charges and moments are identical for paired endpoints. The zero-source Ca/La
preparations already prove one environment reference can be reused per cavity.

No CPCM or second full solvation energy is added. The GK source remains an
approximate projected-monopole representation; exact direct coupling does not
make GK an exact quantum-density functional. Native covalent/damping rules
remain for MM–MM terms. QM–MM direct coupling is explicitly unscaled after the
declared boundary-charge preparation. Fixed cap density and its physical source
support remain an approximation to test, not silently modified geometry.
MACE's jointly trained local readout is not uniquely non-electrostatic; overlap
is a known hybrid approximation. No combined gradient, mechanics or entropy.

## Finite actual calculation inventory

Prepare native tasks for the four representations, with source states
environment-only, Ca and La. Use the common Ca2018-derived metal cavity policy
already qualified. Execute only these named variants:

1. Primary physical geometry and radius:12static energies+12direct-field queries,
   followed by12supplied-field responses at1e−7Debye and12at1e−9Debye.
2. The previously used proper rigid rotation and translation of every atom,
   primary radius:12static energies+12field queries+12responses at1e−9Debye.
   Rotate the stored actual density vector/tensor observations consistently;
   scalar potentials and physical IDs are unchanged. This is a full native
   geometry replay with transformed fixed-density data, not a new quantum
   evaluation or independent validation structure.
3. Common metal GK-radius factors0.95and1.05, original geometry:24static energies,
   24field queries and24responses at1e−9Debye. Change only the source-metal GK
   SOLUTE diameter identically for Ca/La; hold its descreen/other parameters
   fixed. These are prescribed boundary-sensitivity variants, not alternative
   parameters to select based on a sign.
4. Three actual vacuum identity energies on the real GGR-extended geometry:
   environment-only, Ca and La. Set all environmental permanent moments and
   all induced variables to zero and disable GK/nonpolar solvent. Retain real
   physical positions/source moments. The target environment now equals the
   isolated vacuum reference. Verify DeltaU=0; native source internal Coulomb
   may remain nonzero but is excluded from this correction. No response solve
   is needed when the entire response subspace is explicitly disabled.

Total:51native static `energy()` calls,48native direct-field queries and60native
coupled response solves, with no new high-level/model evaluations. Standard
and tight responses reuse the corresponding static fields; no duplicated
static query merely to change a solver tolerance. Partial failures remain
visible and successful outputs are reused only with matching settings/receipts.

Use the existing64CPU standard/memory allocation,64GB and node exclusions.
Previous controls suggest minutes of total kernel work; measure actual build,
preparation, execution, failures and allocation cost. No project CPU/time budget.
Ordinary scoring would need one environment reference plus two endpoints;
these numerical/sensitivity controls are one-time development work. Keep full
upstream DFT/charge/density/MACE costs separate and do not claim matched
production affordability from this incremental solver job alone.

## Frozen numerical, physical and predictive checks

- All native states/cavities/charges/masks match their declared preparation;
  active-site density inputs cover every environment site. Frozen induced
  components are exactly zero. No guessed missing component or baseline fallback.
- Primary1e−7versus1e−9 response refinement changes each endpoint environmental
  correction and every Ca−La contrast by≤0.01kcal/mol. Retain the native unrounded
  final SCF residual in addition to its rounded display; adding this diagnostic
  must not change solver equations or convergence rules. Maximum100iterations
  remains; nonconvergence is a failure, not an increased limit after inspection.
- Rigid full-model endpoint/component and Ca−La changes≤0.01kcal/mol. Identity
  corrections≤1e−8kcal/mol. Direct algebra closure≤1e−7kcal/mol. Preserve actual
  native permanent, GK/nonpolar and induced terms separately.
- Report maximum environment induced-dipole magnitude. Values>1eÅ are flagged
  as very large atomic response requiring separate justification, not clamped
  or silently removed. This diagnostic is not a fitted universal physical limit.
- Radius variants must change each alpha-minus-GGR contrast and the GGR
  partition shift by≤1kcal/mol relative to primary. Raw per-site R may also
  shift by a common reference offset; report it, without mistaking an offset
  for a change in relative predictions. The primary radius stays fixed even
  if a variant appears favorable.
- Primary GGR connected-minus-extended R must satisfy the existing2kcal/mol
  partition gate. The two alpha structures must each exceed both GGR R values
  under the existing direction rule (four differences>0.02kcal/mol, the unchanged
  numerical margin in `mace_explicit_field_short.TOL`, not a newly fitted margin).
  Report all failures and all denominators. GGR is direct same-assay evidence;
  alpha is qualified cross-study evidence. These are two consumed biological
  groups, not four independent or blind affinity observations.

Use the1e−9 primary result for the reported candidate, with1e−7only a numerical
comparison. A pass on these development cases warrants wider testing, not broad
validation, an aquo reference, inherited PQQ bands or automatic promotion.
A failed candidate remains a result; do not retune its radius, boundary, density,
QM region, water, label or threshold to make it win.
