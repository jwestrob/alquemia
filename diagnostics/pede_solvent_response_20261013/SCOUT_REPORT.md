# First DFT endpoint converged; collection requires a parser repair

The sample1 La/water molecular calculation succeeded. Job3354 failed afterward
because the collector looked for `Last Orbital Gradient` in a summary that this
DIIS→TRAH output does not print. The active TRAH header and `TolG` both specify
1e-5; the final Newton–Raphson macro residual is5.663e-7. The intended criterion
is satisfied. No tolerance, energy, coordinate or electronic state is changed
to recover collection, and no repeat molecular calculation is needed.

The original job, failed collection and native output remain immutable.
Collection-only recovery3370 passed15 checks of the actual output with and without the
optional summary, older real TRAH output, and explicitly corrupted parser
fixtures. Missing/conflicting tolerances must still fail. All other original
admission checks must remain identical. See the versioned recovery workspace
linked in [EXECUTION.json](EXECUTION.json).

## Actual molecular evidence

- Native ORCA6.1.1 returned0 and terminated normally after101 SCF cycles.
- The final SCF energy is−6072.995733199557Eh; D4−0.133389313 and gCP+0.165475762Eh
  are printed once before final singlepoint−6072.963646749778Eh.
- Charge−2, restricted singlet,192 atoms,828 explicit electrons, La ECP46 and
  native r2SCAN-3c basis agree with the frozen preparation. The neutral atomic
  La initial-guess subcalculation's multiplicity is not the molecular state.
- Both complete192-atom populations and the GBW/density files were retained.
  The printed CPCM dielectric energy−0.59967299838052Eh is already included;
  it is not `E_water−E_vacuum`. That solvent transfer is still unavailable.
- Some printed one-/two-electron decomposition totals do not reconstruct the
  final SCF total in this output. Preserve that reporting discrepancy. The final
  macro, SCF/property total and virial total agree; the planned comparison uses
  complete final singlepoint differences, not the inconsistent component block.

This is one admissible endpoint, not validation of charge response, orbital
stability, biological discrimination or a repaired solvent model. Population
schemes differ across methods; similar group charges would not prove agreement.

## Measured cost and continuation

Scheduler wall9580s on12 CPUs gives114960 allocated CPU-seconds. Slurm reports
TotalCPU1day07:42:29 and sampled peakRSS18172MiB. ORCA wall is2h39m32.091s.
Actual allocation48000MiB and runtime MaxCore3100decimalMB/rank were sufficient.
The computation was expensive but used its assigned CPU allocation effectively;
this is a development reference, not a proposed routine inference cost.

The collection regression passed. Root submitted the remaining seven declared
endpoints in four independent12-MPI lanes: original sample1 La vacuum, sample1
Ca water/vacuum, sample3 La water/vacuum, and sample3 Ca water/vacuum. Within
each pair, a failed first cell stops the lane before its partner. This is the
same frozen eight-cell experiment, not an expansion; the completed water scout
is reused in place. Exact source/input/protocol hashes match across the separate
workspace manifests. Jobs3371–3374 are running; `REFERENCE_DISPATCH_v1.json` maps each of the eight cells
to exactly one workspace and preserves missing cells.

Four lanes require48 CPUs, below the current60-user limit together with the
three preserved scheduler-listed jobs. Recheck actual allocations before
submission. No memory/time flags, priority changes, exclusive bypass or old
job cancellation. Root arms completion and health wakes for every owned lane.
Production, calibration and earlier failures remain unchanged.
