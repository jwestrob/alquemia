# Eight-cell La/Dy runtime qualification: static review

16 October 2026. Independent source review only; no adapter execution, engine
calls, numerical analysis or submissions by this reviewer.

## Reviewed source

- `DY_QUALIFICATION_PLAN.md`: exact real 48-atom 7CCO La/Dy pair; eight
  native/standalone × vacuum/water endpoints; absolute energy and full-population
  gates retained separately from solvent-component agreement.
- `DY_QUALIFY.py`, SHA256
  `71c16ea9d4f7978c2c4e3f6cc855d8522b7ebc0daecbf2a5fb98eb2d9fa5c8a8`.
- `DY_RUN.py`, SHA256
  `98a81051081333512e99c0e07e2adebb8930f809ce9dab07651360da12cc1438`.
- `DY_RUN.sbatch`, SHA256
  `9b9d59cb4b5f0d1a89b8bf166c81e9c07ba22f4f024f8f56595482be0647cc67`.

No remaining source blocker found for the declared matrix. Configuration and
actual execution remain root-owned and are not accepted by this static review.

## Findings resolved before execution

The initial runner made standalone electron/spin fields optional and checked
only generic convergence and populations. The reviewed version requires actual
GFN2 method/version, effective electron/spin state, printed net charge, smearing,
accuracy, iteration limit, no restart, no external point-charge potential and
the declared solvent/reference settings. It checks normal termination across
stdout and stderr; the pinned helper writes xTB termination to stderr.
It retains high-precision stdout energy and applies the existing 1e−7-Eh
JSON/stdout consistency gate because installed xTB JSON energy is printed at
lower precision. The 1e−6-Eh native/standalone absolute gate is unchanged.
The misleading zero `model_calls` receipt field was renamed `new_MACE_calls`;
actual endpoint attempts remain explicitly counted.

## Contracts checked

The native raw parser's atom/charge/multiplicity/electron/basis schema matches
the runner. Full-precision property populations, all three printed SCF criteria,
normal termination, property/stdout agreement and charge closure are required.
Native diagnostic preflight/postflight verifies the full task-local La+Dy
parameter file, exact runtime import evidence and pinned inputs/dependencies.
The guard remains diagnostic and does not confer research-registry admission.

The standalone helper accepts the staged `atom_count`/charge task schema and
uses the actual `xtbout.json` filename. Both backends preserve the same atom
order and paired source coordinates. GFN2 effective multiplicity one is kept
distinct from physical Dy(III) multiplicity six.

Sequential staging precedes at most four one-thread workers. The native process
environment is set once before worker creation; per-task subprocess environments
do not mutate it. MaxCore is derived from the actual allocation with 25 percent
headroom. The wrapper requests four CPUs without memory or time overrides.

Failures remain in the eight-cell denominator; no retries, offsets, relaxed
parity gates, automatic registry updates or affinity claim are introduced.
Resource receipts preserve each process's wall/CPU time and peak RSS. Existing
helper/schema comparisons were by source inspection, not by running a mock or
new scientific calculation.
