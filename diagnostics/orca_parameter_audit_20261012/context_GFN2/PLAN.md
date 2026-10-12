# Linux GFN2 context/motion factorial — frozen before execution

12 October 2026. Question: does context expansion modify the solvent-transfer
response to the same archived donor motions, and does the current Linux runtime
reproduce the corrected boundary evidence? This is a two-source development
factorial and runtime transfer, not a new 25-reference calibration or repetition
of the completed 2,796-point original audit.

Inputs are the admitted 20 context/geometry rows in
`workspaces/orca_parameter_audit_20261012/context_factorial_preparation_v1/PREPARATION.json`.
Original and expanded contexts each receive the common origin and four archived
Ca/La donor proposals. Shared coordinates match exactly; both endpoint metals
use identical coordinates. No additional geometry, optimizer, water, proton,
PQQ, charge, spin or cavity-parameter change is permitted.

The selected synced archive lacks the corrected per-candidate raw receipts, and
its original `/Users/jacob/Documents/Codex/...` run paths are absent here.
Accordingly the new Linux matrix is explicitly **80 xTB points plus 80 tblite
checks maximum**, including and reusing the four scout points. This supplies
20 geometries × two metals × vacuum/water under a pinned Linux runtime.
Historical corrected scores remain preserved and distinct. Selected archived
energies serve as numerical comparisons only; they do not fill missing points.

## Scout and frozen numerical rules

First run q9z4j7 original context at its original adaptive_La geometry for both
metals and both media. These four exact scored XYZ hashes match the published
edge table. Proceed to the remaining 76 points only if every scout passes fresh
convergence, bookkeeping, xTB–tblite parity, and archived-energy comparison.

Primary: `/home/jwestrob/miniconda3/envs/lanm_bench/bin/xtb`, version6.7.1,
native GFN2 implementation (never `--tblite`). Independent: the isolated local
`tblite_0_7_0/venv/bin/python`, verified package version0.7.0. The compiled
libtblite0.4 linked by xTB is not substituted for this independent implementation.
`RUNTIME.json` pins executable, external parameter snapshot, Python distribution
files, extension modules and dynamic libraries before molecular execution.

Both backends: chargeCa−3/La−2, physical multiplicity1 / unpaired count0,
300K, accuracy0.0001, max500iterations, no restart, named-water ALPB with default
gsolv reference. tblite receives positions in Bohr and temperature
`300*3.166808578545117e-6` Hartree; local0.7 API and official Python documentation
were inspected. Original GFN2 remains spin independent. Electron checks use the
declared GFN2 valence space, not full atomic electron totals.

Acceptance requires process success and solver convergence, finite energies,
one finite charge per input atom, charge closure within1e−5e, and tblite
occupation sum matching the declared valence electron count within1e−5e.
Per-input backend energy difference must be ≤1e−6Eh and maximum all-atom charge
difference ≤1e−5e. Where the exact archived standalone selected energy exists,
the primary Linux energy must also agree within1e−6Eh. These tolerances precede
all new molecular outputs. Continuous residuals and raw output are retained.

Failed primary or parity checks remain unavailable; no retries, relaxed
tolerances, alternate occupation guesses or hidden successful fallback.
The scout and final collections keep every declared cell. The root agent owns
MACE evaluation and composite integration; this module supplies GFN2 components.

## Execution and outputs

Jacob explicitly authorized ≤16 aggregate login numerical threads; this agent
owns **two total**, implemented as at most two simultaneous workers each running
one one-thread xTB or one one-thread tblite child. The two backends run serially
within a worker. OMP/BLAS/MKL/NumExpr/VECLIB are set to1, nested threading disabled,
and each owned process group must disappear before a worker slot is reused.
This allocation is not a Slurm template. Current CLAUDE.md/cluster instructions
were read. No Slurm submission, MACE, ORCA, optimization or DFT is performed here.

Every point retains exact input, commands, full stdout/stderr, actual wall/CPU/RSS,
process-group completion, energies, charges and parity results. `manifest.json`
and `RUNTIME.json` are under
`workspaces/orca_parameter_audit_20261012/context_gfn2_v1/`.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/context_factorial_gfn2.py execute \
  --directory workspaces/orca_parameter_audit_20261012/context_gfn2_v1 --stage scout
```

The same command with `--stage remaining` refuses execution unless the saved
four-point scout passes. Existing point results are immutable. Final numerical
availability and scientific class validity are distinct; these two inspected
sources alone cannot define a calibrated discriminator.

Documentation: https://tblite.readthedocs.io/en/latest/api/python.html .
