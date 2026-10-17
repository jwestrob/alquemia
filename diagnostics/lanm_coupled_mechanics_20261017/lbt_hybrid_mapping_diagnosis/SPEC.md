# Exact origin-map failure diagnosis

## Question and immutable evidence

Why does frozen hybrid evaluator `a482559c` pass zero-call preflight 12972 on
`node-64-768g-1`, then fail its exact mapping comparison before any molecular
call in scout 12975 on `node-224-2t-8gpu-1`? This diagnostic tests the actual
saved four origin partitions; it does not presume floating-point noise.

Bind the byte-identical hybrid config, both actual application results,
snapshots of terminal receipts, both launchers, the original partition result
12953, preparation implementation/config, exact cap-gradient/common helper,
and the definition-only GFN backend with its import-time thread environment.
Preserve the original evaluator, preparation maps and failure without edits.

## Finite computation and environment

Run only in Slurm on the failed host, using the failed launch's literal MACE
venv Python. Use OMP_NUM_THREADS=12 and OPENBLAS/MKL/BLIS/NUMEXPR_NUM_THREADS=1,
as in that launch. No model, molecular engine,
energy, force, optimization or GPU operation is requested. NumPy is imported
only during the Slurm run. Do not call preparation.initialize or instantiate
HybridEvaluator. Import the pinned preparation and cap-gradient definitions
(plus the latter's common dependency), then the pinned GFN backend after
NumPy, matching the evaluator's import order. The backend's module-level
THREAD_KEYS loop changes OMP and BLAS variables to 1; preserve this actual
side effect without invoking any backend function. Record the initial and
post-import environment separately. Twelve allocated CPUs match the failed
launcher; the actual mapping work requires no parallel molecular workers.

For each source/region, call the original `update_geometry` twice, injecting
the same NumPy and original `cap_jacobians` dependencies used by the evaluator.
Use all 256 exact archived physical coordinates. Save both complete generated
maps and enumerate every unequal or type-different leaf versus the complete
archived map. Report old/new values, exact float hex, absolute/signed difference,
and float64 ULP distance, with separate coordinate and Jacobian maxima.
Preserve both the evaluator's exact Python dictionary equality and stricter
typed leaf diagnostics. Repeat differences test process-local determinism.

For each of eight cap instances, preserve definitions, physical endpoint
coordinates, derived vector, archived/regenerated cap rows and norms computed
with the actual NumPy norm, dot, sum and Python fsum/hypot primitives. These
are explanatory measurements, not replacement maps or new admission rules.

Record host, Python binary/environment, NumPy module/version/build/runtime
dispatch information, CPU information/affinity and scheduler/thread variables.
All four rows remain in the denominator if any one diagnostic row fails.

## Interpretation and limits

This diagnostic makes zero molecular calls and changes no tolerances. A
completed diagnostic means all four comparisons were obtained, not that any
mapping or method has been admitted. Type, identity, ordering, coordinate,
charge or Jacobian changes remain distinct. ULP-sized differences alone do
not authorize a fix. Root reviews actual differences before a separate,
versioned recovery or additional cross-host test. Passing preflight on the
other host is evidence of the historical exact gate; this run does not
reproduce that other process environment.

## Commands (root owns execution)

Metadata only: `python DIAGNOSE.py --stage bind --output FRESH_STAGING`.

In the specified Slurm allocation and environment:
`MACE_VENV_PYTHON DIAGNOSE.py --stage run --config FRESH_STAGING/CONFIG.json --output FRESH_EXECUTION`.

No output path may already exist. Config/source/spec and all declared input
pins are checked before mapping arithmetic. Outputs are immutable new files.
