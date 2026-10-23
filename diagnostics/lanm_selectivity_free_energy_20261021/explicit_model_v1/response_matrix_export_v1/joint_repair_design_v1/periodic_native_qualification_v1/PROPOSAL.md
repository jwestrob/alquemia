# Periodic native qualification: staged implementation proposal

Status: prospective source preparation only. Root intent `e4dc6021` authorizes the isolated operator harness; no build, native application, factorization, energy/force evaluation or guard removal is released. Exact authorities are in `SOURCE_AUTHORITIES.json`; Arendt owns `NUMERICAL_CONTRACT.md`.

## Fixed identity and question

Use only the old terminal La candidate `e41fcd3d…` (artifact `fdf42081…`), composed BIO manifest `d8622d40…`, and the literal 11958-site La geometry/cell from the admitted static input `4b166082…`. The current calibration is a different candidate. This is numerical qualification of a fixed Hamiltonian, with no claim of physical accuracy or molecular-dynamics readiness.

## Stage 1: complete operator parity, guards intact

Add a standalone program and one derived setup entry outside the immutable library. The setup entry preserves the checked `periodicresponsesetup` body through `periodicresponsegeometry`; its only semantic delta is returning to the explicit-operator harness instead of calling the backend guard. Preserve the original setup and both fatal guards byte-for-byte. This new entry must never call `induce`, `energy`, `gradient`, `analysis` or PCG. The normal application remains disabled.

Reuse the admitted PME-export library (`61e2fab0…`, source manifest `d54c0964…`), native `dfield0c`, `ufield0`, `grid_uind`, FFTW, qfac and `fphi_uind`. Export native h, beta, alpha, atom/type/xyz, active map, coefficients, direct RHS and initialized PME state before actions. Require correspondence with the exact static manifest/H and raw input lineage. No H imported into the native operator as a substitute for testing its calculation.

For C=diag(sqrt(alpha)), native `ufield0` returns (T−H)Cz. Export the normalized complete action z−C field and compare against S z + R z, where saved lower CSC S is real+self+H and R is the independent full-complex packed reciprocal replay. Do not symmetrize qfac or assume equal dual channels. The separate reciprocal reference uses the admitted primitive chain and its sign is positive in the normalized matrix.

Use the ten fixed probes from `BACKEND.ten_probes`: three uniform Cartesian, three metal coordinate basis, four fixed-seed Rademacher. For each probe j evaluate channel pairs (z_j,0), (0,z_j), (z_j,z_next), with cyclic next. Thus one program has **30 complete actions plus 30 separate reciprocal-reference actions**; the complete calls already contain 30 reciprocal applications, so the total is **60 reciprocal primitive applications**. No induced solves or energy/force calls. Probe arrays are serialized from the frozen helper, not regenerated with a Fortran RNG. Retain raw channel labels, normalization/signs and unit metadata.

Use the inherited sampled action gate 1e−11*max(1,max|predicted|), separately for each channel and comparison. Compare full complex reciprocal replay and sparse real action independently. Check persistent qfac/splines/alpha/h/H/maps, restore u/v and scratch between calls, and preserve failed output. These are finite probes, not a uniform application-error bound: e_app stays unknown.

## Stage 2: separate implementation/release after actual parity

A subsequent isolated finite-evaluation wrapper must bind a per-geometry certificate to literal inputs, native setup state (coordinates/cell, alpha, damping/scales, h/H, coefficients), PME split/grid/qfac/splines, source and runtime. Any changed state invalidates the certificate. The normal guard must reject missing/stale/mismatched bindings; PCG convergence is not a replacement certificate. Root must review this new guard implementation separately.

Reuse native `energy`/`epolar0e`, `epolar1` gradient and frame/PME derivatives, and the common local-response correction. The `epolar3` analysis route has a separate correction path and needs a declared same-state comparison if added; it is not silently included in the current 36-program count. A finite harness calls one gradient per program and saves its energy, all Cartesian gradients, virial, both unequal induced channels and direct RHS. Central readback verifies shared h on both RHS, complete-operator residuals, and the stationary dual energy contraction. These add reporting/checks, not programs. Do not use all-Cartesian `testgrad` on this full box.

The prospective La-only panel is 1 active center + 32 endpoints (8 directions × 2 steps × 2 signs) + 3 center baselines (admitted-parent absent, new absent, new explicit zero): **36 E/F programs**. Directions are three metal translations, one deterministic nearest-water rigid rotation, one fixed whole-box tangent, and three diagonal strains. Arendt freezes their literal geometry/tangent rules. Steps are ±1e−4 and ±5e−5 Å; coordinate and strain slopes use the inherited .001 kcal/mol/Å gate plus literal-token uncertainty. Only these projections are derivative-qualified despite full gradient output. Off-diagonal strain, Dy, trajectories and transfer are outside this panel.

Reuse the admitted center certificate; each of the 32 displaced states needs a new PME initialization/export, real+self+H assembly and one fixed nominal shifted factor with the existing ten-probe audit. With separate initialization processes this adds **32 operator-only native programs**, **32 assemblies**, **32 factor attempts and 320 audit solves**. Combined with Stage1 the two stages total **69 native programs** (33 operator-only +36 E/F), without hidden internal Cartesian FD. Failed certification stops that profile and does not establish instability. No automatic retries or alternate factorization settings.

## Cost and remaining implementation

The admitted center factor took 75.577 s wall and 486.180 CPU-s; its assembly took 13.464 s wall. A serial projection of 32 factors alone is about 40.3 minutes, excluding assembly, initialization/replay, I/O and E/F work. This is a finite qualification budget, not an MD rate. Root may schedule independent profile certificates in parallel under explicit resource limits after parity; do not make one thread wait serially by design.

Minimal source work now: isolated checked-setup derivative, standalone parity harness, immutable input/probe preparation, thin raw-output comparison reusing sparse and reciprocal readers, source controls and a compile/link recipe. No broad electrostatics reimplementation, dense matrix, new Hamiltonian or full-factor rerun is needed for Stage1. Stage2 requires a separately reviewed exact-state certificate gate and finite-gradient reporter before any E/F execution. Rigorous e_app remains optional under `NOMINAL_CONTRACT` and is not a new blocker to ordinary-double numerical testing.
