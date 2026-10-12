# MACE cross-context singlepoints — complete

All40 logical endpoint energies are available:24 exact archive reuses and16 new
cross-context calls. Four additional CPU replay calls match archived GPU origins
exactly (0.0kcal/mol error; frozen tolerance0.01). This qualifies the execution
route for this finite comparison, not biological discrimination. No optimization,
new weights, neutral-charge adapter, force calculation or threshold was introduced.

The original native OMOL checkpoint/head, float64 dtype, physical total charge
and multiplicity convention are unchanged. Reuse verifies complete model/state,
pinned receipts and numeric coordinates (maximum serialization difference
5.56e-17angstrom). Both metals use the same five physical geometries per context.

Actual20 CPU calls completed in113.884s including24.096s model load, with
590.696 measured process CPU-s and4,451,228KiB peak RSS; zeroGPU. A warm model
ran on an eight-CPU affinity with torch intraop8/inter-op1. The process created
15 OS threads (two pools), constrained to eight CPUs; numeric concurrency never
exceeded root's eight-thread allocation. Other agents used at most4+2+1 threads.

`RESULT.json` contains all40 unrounded energies and provenance. The final
context/geometry interpretation requires the independently checked GFN2 solvent
component and is reported by the integrator. MACE-only data do not repair the
failed corrected accommodated calibration.

Reproduce into a new output directory:

```bash
python scripts/context_factorial_mace.py prepare \
  --preparation workspaces/orca_parameter_audit_20261012/context_factorial_preparation_v1/PREPARATION.json \
  --output workspaces/orca_parameter_audit_20261012/context_mace_replay
OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/context_factorial_mace.py execute \
  --manifest workspaces/orca_parameter_audit_20261012/context_mace_replay/manifest.json
```

Only run under a current allocation or explicit login allowance, keeping aggregate
threads within that allowance. Existing output directories are never overwritten.
