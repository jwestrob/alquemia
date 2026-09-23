# Scaled four-angle completion: numerical result

The same four-angle physical model now produces **60/60 valid original30
candidates**, including the prior MMOL1770 La search failure. These30 include
25 canonical references,3 consumed crystal controls and2 unlabeled PLM cases.
This completes
the numerical prerequisite for a fair, separately calibrated adaptive-family
comparison. It does not yet establish improved discrimination. Root owns the
subsequent common-pool scoring and canonical-only calibration.

## What changed

Protocol `common_four_angular_native_OMOL_SLSQP_hartree_units_v2` expresses the
shifted native-MACE objective and its analytic gradient in fixed
Hartree-equivalent units instead of eV. The same positive conversion multiplies
the previous objective stopping tolerance. It changes SLSQP's numerical path,
not the physical energy, selected coordinates or feasible set. The single
SciPy tolerance also makes some internal non-energy checks stricter; this is
an explicit numerical version, not a claim of identical optimizer trajectories.

The factor is 0.03674932217934773; the scaled stopping tolerance is
3.6749322179347735e−11. All original states, paired four-mode selections,
analytic forces, source atoms, charges and frozen coordinates remain identical.
The unchanged constraints are ±0.8 rad per selected angle and at most 0.8 Å
final source-heavy-atom displacement. Unchanged chemical/overlap guards apply
to intermediate evaluations. There is no alternate scale, restart, extra
starting geometry, clipping, synthetic barrier energy or relaxed guard.

## Actual execution

| Job | New searches | Allocated wall seconds | CPUs | Allocated CPU-s | Allocated GPU-s |
|---|---:|---:|---:|---:|---:|
| 1209963, pilot | 8 | 69 | 32 | 2,208 | 69 |
| 1209968, completion | 52; 8 pilot results reused | 444 | 32 | 14,208 | 444 |
| Total | 60 | 513 | — | **16,416** | **513** |

Both allocations completed. Native q0 energies and Cartesian forces were reused
for all 60 endpoints. No new GFN2 or DFT call belongs to this proposal adapter.
Worker intervals were 50.180121283978224 and 419.66119112074375 seconds; scheduler
totals above include launch/collection overhead. Per-step measured MaxRSS was
1,171,388/1,174,884 KiB for the GPU worker steps and 82,332/80,936 KiB for the
batch steps. These are scheduler process-memory observations, not GPU VRAM
measurements or whole-allocation peak-memory sums.

All 60 optimizer exits report success, and all native proposal energy changes
are nonpositive. Maximum iteration count is 60. Final heavy displacement reaches
0.800000000020692 Å, within the unchanged 1e−7 Å acceptance tolerance.
**27/60 candidates reach a displacement or angular boundary.** They are bounded
geometry proposals, not demonstrated unconstrained minima or equilibrium states.

The trace contains 205 completed model requests outside the final 0.8 Å domain;
some requests can reuse an already evaluated point. Maximum intermediate heavy
displacement is 1.6336182181699836 Å. The original plan explicitly constrains
final admissibility, not every intermediate SLSQP evaluation. All evaluated
points retained the independent chemical/overlap guards. No invalid final
candidate was projected into acceptance.

The MMOL1770 Ca/La native works are −15.82932377784/−21.96987040521 kcal/mol;
both end on the 0.8 Å boundary. The old unscaled La failure is preserved.
These works alone are not composite scores or evidence of classification gain.

## Traceable artifacts

Under `workspaces/adaptive_completion_20260922/`:

- `FIRST_STEPS_v1.json`: pre-execution geometry-only probe of all60 actual q0
  gradients; zero nonzero-point energies/gradients returned. All60 first trial
  geometries pass chemical guards;45 lie inside the final displacement domain.
- `original30_v1/manifest.json`: SHA256
  `51e8d8f458a9cc06b7b3bfb309ca19193d9c3c7d213c55de4629e14a57ca80f7`.
- `original30_v1/after_proposals_1209963.json`: immutable partial8/60 collection.
- `original30_v1/after_proposals_1209968.json`: actual60/60 final collection.
- `original30_v1/TERMINAL_SUMMARY.json`: all60 terminal rows, residuals,
  boundaries and native work; `ACCOUNTING.txt`: actual scheduler receipts.
- Each endpoint retains its immutable origin, proposal, evaluations, constraints,
  analytic force arrays, optimizer exit and intermediate request trace.

The older terminal-only proposals, unscaled angular results and separate
joint-metal experiment remain distinct. No production default changed.

## Primary225 transfer status

Preparation only is complete. The same optimizer and physical selector have
410 executable tasks for205 base-available cases. All225 declared cases and450
endpoint statuses remain in the manifest/collection. The20 inherited unavailable
cases (17 source preparation,3 numerical base failures) get no new search and
remain in the denominator. No fold/source was selected by score or class.

`primary225_v1/manifest.json` SHA256:
`cbe7cc25c62bd4f2012488f1db31ae0263aa5f919418545d563ff20b001e8751`.

`primary225_v1/collection_unrun_v1.json` explicitly records zero candidates.
Execution requires the parent's separately frozen scaled-family reference.
The old terminal-only q0 forces are reprojected through the same physical maps
onto the selected four angular primitives; their old1/2-component gradient
vectors are not misused as four-component gradients. Both original force
records and the exact new projection are retained. These reference folds are
already consumed development evidence, not prospectively blind tests.

## Tests

Seven real-artifact tests pass in62.448 seconds, zero skips
(`TESTS_v2.txt`). They check the actual30/60 source population and fixed pilot,
positive objective/gradient/tolerance unit algebra, all60 geometry-only initial
steps, actual MMOL1770 q0 reuse without a model call, full primary225 preflight,
all410 unchanged-force reprojections, and explicit225/450 unavailable/unrun
statuses. These replay tests make no new molecular calls. The independently
executed60 native-MACE searches are scientific results, reported above.

**Recommendation:** finish the separately calibrated common-pool comparison
before interpreting or expanding this family. Numerical completion is now
successful; accuracy, transfer value and composite scoring cost remain separate
questions. The baseline is unchanged.
