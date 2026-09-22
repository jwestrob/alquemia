# The unchanged native solver converges with four additional iterations

**Three real native GFN2 single points passed the declared diagnostic.** The
failed 4MAE La-at-adaptive-Ca vacuum cell converges after 129 cycles when only
`MaxIter 500` is added. The original 125-cycle failure remains intact. Both exact
origin controls reproduce the primary energies far inside the frozen 0.05
kcal/mol tolerance; the native parameter exports and electronic states agree.

| Vacuum endpoint | SCF cycles | Energy / Eh | Control difference / kcal mol−1 |
|---|---:|---:|---:|
| 4MAE La at adaptive-Ca candidate | 129 | −351.932114159458 | unavailable: primary failed |
| 1H4I La origin | 49 | −266.813855889963 | 0 |
| 4MAE La origin | 56 | −351.931156450368 | +0.0000000018905 |

Each calculation terminated normally with the native GFN2 LEANSCF driver, special
xTB mixer, fresh SAD guess, 300 K smearing, unchanged NoAutostart and TolE=1e−6 Eh.
All source XYZ files are byte-identical, and removing the single added MaxIter
line reproduces each original input exactly. State/valence accounting, atomic
charge closure and exported GFN2 parameter equality pass. Singlet charges remain
−1 for 1H4I and −2 for both 4MAE geometries. No alternative orbital guess,
Hamiltonian, geometry, tolerance or smearing was used.

## Meaning and limits

The original failure exhausted its iteration allowance just before convergence;
it was not an unsupported-element or input-state failure. Its last unconverged
iterate was never used as an energy. This result supports considering an explicit
numerical-recovery overlay for this one missing cell. **No overlay or primary
collection was changed by this diagnostic.** The parent owns that decision and
the resulting common-pool comparison. These three points do not establish a new
classifier's accuracy or generally guarantee SCF convergence.

The 0.05 kcal/mol control tolerance was frozen before execution from the existing
four-term 0.20 kcal/mol numerical score policy. Control reproduction does not
independently bound error at the formerly failed geometry. There was no further
iteration-limit escalation, new DFT call, geometry proposal or calibration.

## Actual execution and records

Job **1209876**: three concurrent eight-rank calculations on
`node-224-2t-8gpu-1`; 24 allocated CPUs, 64 GiB host request, zero GPUs.
Scheduler wall time **57 s**, **1,368 allocated CPU-s**, zero GPU-s; actual total
CPU time 679.146 s. Runner interval 54.675563 s (1,312.213520 allocated CPU-s),
with per-task execution intervals 54.572506/35.503558/41.442373 s in table order.
Reported batch MaxRSS is 2,181,704 KiB; this scheduler field is not asserted to
be a simultaneous sum of memory across all MPI ranks.

All three actual-source recipe/state/parser checks and the existing runner's
path/pin/MPI preflight passed. These are real executed scientific checks, with no
fabricated outputs. Full collection, receipts and original failure pins:

`workspaces/adaptive_maxiter_20260922/prepared_v1/collection_1209876.json`

Compact exact data: [RESULT.json](RESULT.json). Frozen method:
[PLAN.md](PLAN.md). Reproducible read-only validation/collection:
[COMMANDS.md](COMMANDS.md). Baseline, production default and original adaptive
pool remain unchanged. Do not rerun this completed three-point diagnostic.
