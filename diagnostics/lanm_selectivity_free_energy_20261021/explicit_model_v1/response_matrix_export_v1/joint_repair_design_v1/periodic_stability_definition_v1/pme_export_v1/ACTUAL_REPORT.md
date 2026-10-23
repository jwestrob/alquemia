# Initialized-PME export: build and finite native qualification

The isolated import-only recovery builds successfully, and the exact two-program
qualification passes. The original failed build and source remain preserved.

## What passed

The unchanged consumed neutral six-atom water dimer was evaluated with the
predeclared 44 Å orthorhombic conducting cell, 32³ FFTW grid, order-5 splines and
Ewald alpha 0.3 Å⁻¹. Exporter OFF and ON give byte-identical full-precision energy
and all 18 Cartesian gradient components (energy 5.27155497953516505 kcal/mol;
gradients in kcal/mol/Å). The initialized transform, positions, identities,
polarizabilities, grid indices, spline values and every qfac value used by the
native reciprocal applications exactly match the snapshot.

One standalone harness link and two native programs completed. Each program made
one analytic gradient call. ON then made exactly 54 internal reciprocal
applications: 18 first-channel basis columns, 18 second-channel basis columns,
and 18 unequal packed-channel columns. These applications are separate from the
two energy/gradient evaluations. Model induced dipoles, qfac and splines remained
unchanged. No displaced-energy calculation, optimization, full-box call or
periodic-response execution occurred.

All 54 applications agree with the complex FFT replay of the exported actual
weights. Maximum component difference is 1.0408340855860843e-17 against the
predeclared 1e-11 scale-one tolerance. Channel leakage is 5.7959690062212124e-18.
The replay's 66,417 checks are literal scalar/input joins, not independent
physical tests. The harness tests positive normalized reciprocal R; native
umutual1 subtracts its dipole field, giving the −R field contribution.

## Initialized-state norm evidence

| Quantity | Native/replayed value |
|---|---:|
| Maximum absolute scatter column sum | 0.8871534942241106 |
| Maximum absolute scatter row sum | 0.5726795537472118 |
| Squared Frobenius scatter norm | 0.5314889363107518 |
| qfac minimum / maximum | 0 / 0.0069466602009871535 |
| qfac zero mode / negative entries | 0 / 0 |
| Maximum q(k)−q(−k) magnitude | 3.0357660829594124e-18 |
| Nominal U_R = M max|q| min(c1*c∞,cF) | 115.64754981381782 |
| Nominal evenness perturbation bound | 2.5269618719241316e-14 |

The native five norm quantities and five qfac summaries pass the declared 1e-12
scale-one replay tolerance. No multiplier was clamped or symmetrized. These are
ordinary-double numerical estimates for this exact initialized state. Parity is
not a PSD certificate or a full D stability test, and the tiny nonzero stored
asymmetry remains explicit. This bound is not a value for the large BIO/salt
system or a future geometry. Both periodic-response fatal guards remain intact.

## Preserved failure and costs

The first compile failed because the exporter imported `type` from nonexistent
`atmtyp`; the native authority defines it in `atoms`. Root separately authorized
one recovery changing only that import. source_v2/build_v1 and their raw failure
are unchanged; source_v3/build_v2 hold the corrected isolated copy. The earlier
pre-copy preparation KeyError is also preserved; it made no compiler/native call.

| Phase | Wall s | Child CPU s | Controller CPU s |
|---|---:|---:|---:|
| Failed first build | 181.919246 | 169.491973 | 0.428984 |
| Successful recovery build | 281.334379 | 265.240784 | 0.748780 |
| One link + OFF/ON + saved replay | 4.945719 | 4.099371 | 0.481389 |

All phases used one worker on biotite with truthful job_id:null and no scheduler
allocation. The table totals 440.491281 CPU seconds; small source preparation and
source-control costs are separately recorded in their immutable receipts. Peak
child RSS was 306,952 KiB in the successful build and 351,232 KiB in the finite
qualification (including compiler and saved replay). Every native/link process
completed under its 120-second wall limit. No retry, second candidate, or
production replacement followed.

## Evidence

Paths below are relative to the project root; ACTUAL_RECEIPT.json pins them all.

- W/J/periodic_stability_definition_v1/pme_export_v1/build_v1/execution/RESULT.json
  (`d4d15cda…`), preserved failed build.
- W/J/periodic_stability_definition_v1/pme_export_v1/source_v3/MANIFEST.json
  (`d54c0964…`) and build_v2/execution/RESULT.json (`4cf176c0…`).
- W/J/periodic_stability_definition_v1/pme_export_v1/qualification_v1/run_v1/
  ROOT_RELEASE.json (`87716e11…`), CONFIG_RUN.json (`999fcac4…`),
  execution/RESULT.json (`fc5e5cf4…`) and execution/replay/RESULT.json (`af200ca3…`).
- D/J/periodic_stability_definition_v1/pme_export_v1/qualification_v1/
  ROOT_SOURCE_REVIEW.json (`3a7caa84…`) and READY.json (`eaf2df29…`).

W = workspaces/lanm_selectivity_free_energy_20261021;
D = diagnostics/lanm_selectivity_free_energy_20261021;
J = explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1.
