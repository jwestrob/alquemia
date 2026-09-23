# Six common angles from existing four-angle proposals

**Status: feasibility verified on real archives; optimization and scoring not
executed.** Nine consumed development sources, one deterministic proposed rule.
This is a source-robustness/mechanism test with old-band transfer, not a new
calibration or independent biological validation.

## Fixed input population and selection

All five C5AXV8 La-conditioned samples, including canonical sample 1, plus
A0A3F2YLY8 Ca sample 1, A0ACD6B9F2 Ca sample 4, and A8R3S4 La samples 1/3.
No unsupported Ca-conditioned C5AX structure is relocated or substituted.
All nine have exact prior successful paired proposals and complete strict-native
three-geometry pools. Source, protonation, water inventory, charges, metal/PQQ
positions, context membership and caps remain their own original values.

Extend `adaptive_force_diagnostic.preview(..., maximum=6)` at the **common q0**,
using the actual archived Ca/La Cartesian MACE forces. This preserves its
alternating differential/individual normalized-load ranking, two-pass
Gram–Schmidt independence rule (`rank_tolerance=1e-5`) and deterministic ID ties.
The first four IDs reproduce the old selection in all nine cases. The new
primitive set is common to both metals; it does not depend on a desired class
or require any residue by name. Each source has seven or nine eligible angles.

| Source | Added fifth mode | Added sixth mode |
|---|---|---|
| C5AX La0 | Glu198 χ1 | Glu198 χ2 |
| C5AX La1, canonical | Glu198 χ3 | Glu198 χ1 |
| C5AX La2 | Glu198 χ1 | Glu198 χ2 |
| C5AX La3 | Asp319 χ2 | Asn275 χ1 |
| C5AX La4 | Asp319 χ2 | Glu198 χ2 |
| A0A3 Ca1 | A/197 χ3 | A/285 χ2 |
| A0AC Ca4 | A/204 χ1 | A/204 χ2 |
| A8 La1 | A/300 χ2 | A/300 χ1 |
| A8 La3 | A/300 χ1 | A/350 χ1 |

Minimum new-mode independent geometric fraction is 0.76996, far above the
existing rank tolerance. Selection is viable without near-degenerate additions.
For C5AX La3, added-direction normalized gradients at the actual old proposals
are Ca (−2.047, +6.591) and La (−4.849, +5.896) kcal/mol/Å. Those residuals
justify asking about missing accommodation; they are vacuum-MACE gradients,
not solvent-aware selectivity gradients or predicted score improvements.

## Warm start and admissible region

Each metal starts once from its **own** exact archived four-angle proposal,
with the added two primitive coordinates initially zero. All 18 mapped start
coordinates replay the archive within 1e-12 Å. Reproject their saved full
Cartesian forces into all six directions; padding the old four-gradient vector
with zeros would be incorrect. The original q0 energy/forces are also cached.

Keep the source-q0 reference for every bound: each angle in [−0.8,+0.8] rad,
maximum physical source-heavy displacement 0.8 Å, existing final numerical
tolerance 1e-7 Å. Do not recenter the trust region around the warm start.
Four warm starts already touch this boundary (C5AXLa3 Ca, A0A3Ca1 Ca, and
both A8 La endpoints); the others have extents 0.235–0.745 Å. SLSQP intermediate
points may violate the final displacement constraint, as in the existing
protocol, but actual infeasible trials and geometry failures remain recorded.
No clipped point, fake energy/gradient, relaxed geometry guard or automatic
retry can become a successful candidate.

Native-MACE objective/analytic gradient retain Hartree-equivalent scaling,
SLSQP `ftol=1e-8`, `maxiter=200`, source/cap/bond/overlap checks, and the same
maximum final energy increase rule. The warm-start reference additionally
allows direct reporting of whether each successful search actually improves
its prior native energy. All endpoint statuses and residuals remain explicit.

## Shared pool, exact reuse and finite scope

Score the same five geometries under both metals: q0, prior Ca proposal,
prior La proposal, new six-angle Ca proposal, new six-angle La proposal.
Use native MACE plus the exact qualified fresh native GFN2
`ALPB(water) − vacuum` recipe: TolE1e-10, NoAutostart/no seeds, 300 K electronic
smearing, MaxIter500, rank1. Both media use the same physical state.

- **108 existing strict scalar cells** and all existing MACE cells are reusable.
  Actual receipt/output/manifest pins and state/coordinate/MACE consistency were
  checked for all nine prior pools.
- **18 single-start searches**, not 18 total MACE evaluations. Search evaluations
  are data dependent, with the existing 200-iteration endpoint stop.
- At most **18 cross-metal MACE evaluations** and **72 fresh strict GFN2 cells**
  (two new geometries × two metals × two media × nine sources).
- No q0 recomputation, new DFT, new solvent forces, new starting structures,
  entropy terms, chemistry states, reference fit or full-cohort rescore.

Retaining old candidates makes the mathematical endpoint row minima no worse;
it does **not** guarantee an improved R or class. Operational candidate choice
must retain the existing tolerance/tie rule separately from the mathematical
minimum. Missing new candidates remain failed/missing extensions; an unchanged
old pool must not be reported as successful six-angle optimization. Report all
nine outcomes, old and expanded energies, component work, choices, boundaries,
displacements, residuals and execution cost. Transfer the existing strict32
bands only as a labeled development comparison. Do not tune these bands.

## Implementation and cost feasibility

`Kinematics` and `accommodation_proposals.MACEProposal` already support arbitrary
active-index lists and correct cap Jacobians. However,
`adaptive_angular_proposals.full_q` and the completion optimizer hardcode four
coordinates. A separate private six-angle adapter is therefore required for
shape checks, start vector, six-component cached gradients, constraints and
collector validation. This should reuse the same physical mapping and SLSQP
constraint formula, not mutate pinned production or current four-angle modules.
Verify the six-column constraint Jacobian on actual warm-start geometries before
any launch. The physical inequality remains
`0.8² − |r_i(q) − r_i(q0)|² >= 0`, with its exact mapping derivative.

The existing 18 four-angle searches took **65.55 summed worker seconds** and
260 objective evaluations; their maximum endpoint time was 6.07 s. The 108
reused rank1 strict scalar receipts sum to **5001.52 elapsed seconds**, median
23.08 s, maximum 141.82 s. Scaling that observed cell mean to 72 new cells gives
about 3334 summed rank1 seconds, before preparation/collection/allocation idle
time. Warm six-angle search may take more or less work; this is only a cost
class estimate, not a guarantee or scheduling estimate. Existing allocation
conventions would be one H200/32 CPU/200000 MiB for the warm proposal worker,
and CPU-only 32 one-rank workers/32 CPU/64 GiB for strict scalar scoring.

## Reproducible feasibility inventory

`workspaces/c5ax_response_20260923/SIX_MODE_FEASIBILITY_v1.json`, SHA256
`d580f4bf6ede0b5fd312602475832b4e29c6b481c7e5475dcbac26702e6b18af`, contains
the exact nine source IDs, selections/indices, source tasks, all 18 archived
warm-start receipts/full q/forces, strict 108-cell receipt pins and future count.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
scripts/audit_six_mode_feasibility.py \
  --c5ax workspaces/c5ax_response_20260923/AUDIT_v2.json \
  --strict225 workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json \
  --output workspaces/c5ax_response_20260923/SIX_MODE_FEASIBILITY_replay.json
```

This command projects archived forces and verifies mappings only. It cannot
submit a job, optimize a structure, call MACE/GFN2 or change a score. The next
decision is whether to implement this contained experiment, not whether to
promote a new scorer.
