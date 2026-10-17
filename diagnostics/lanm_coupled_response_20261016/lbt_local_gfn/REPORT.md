# Local LBT source response: MACE and standalone GFN2

Job **12825 completed the four-cell comparison** using exactly two new 7CCN
GFN2 vacuum single points. Both 7CCO GFN cells from job 12488 and all four
matched-core MACE cells from job 12471 were reused unchanged.

All cells use the same prepared 48-atom cores, charge −1 and paired La/Dy
coordinates. S = E(7CCN) − E(7CCO), calculated within each metal before unit
conversion. Values below are copied from the saved result, in kcal/mol.

| Method | S(La) | S(Dy) | S(Dy) − S(La) |
|---|---:|---:|---:|
| MACE | +77.640037 | +3.950156 | −73.689881 |
| GFN2 | −14.563185 | −22.508206 | −7.945021 |

**The models disagree on the source ordering for each metal, especially La.**
Both methods give a negative metal-dependent strain difference, with very
different magnitudes. Neither the common sign nor the flatter GFN response
establishes accuracy. These are local vacuum source strains, not binding
affinities or the whole-peptide composite exchanges. The independent four-cell
PBE0 reference remains a separate running experiment; no PBE0 conclusion is
available here. Its interpretation is declared in
[LOCAL_REFERENCE_INTERPRETATION.md](../LOCAL_REFERENCE_INTERPRETATION.md).

## Actual admission and cost

Both new endpoints returned exit 0, normal xTB termination and empty owned
process groups. The native outputs identify GFN2-xTB 6.7.1, 144 valence electrons
and zero unpaired electrons. Saved charge-closure residuals are 0 for La and
approximately 2e-8 e for Dy, within the unchanged 1e-5 e gate. The source,
runtime, parameter and paired-state preflight passed. Both stderr files retain
IEEE underflow/denormal notices alongside normal termination; no retry or
setting change was used.

Physical metadata remains La multiplicity 1 / Dy multiplicity 6, separately
from GFN's effective singlet. No solvent, optimization, gradient request,
new MACE call or new DFT call was made by this job.

Slurm reports **COMPLETED, exit 0:0**, **6 seconds on 2 CPUs**, no GPU, and
5.812 seconds of measured CPU time. The two molecular wall times were
0.493 seconds for La and 0.512 seconds for Dy. Application time was 5.338
seconds. Exact per-process CPU/memory measurements and scheduler records are
in [ACCOUNTING.json](ACCOUNTING.json). [RESULT.json](RESULT.json) is a byte-for-byte
copy of the complete execution result, with all four cells and source pins.
