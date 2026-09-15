# Completed environmental pilot: current challenger fails physical checks

**All scheduled calculations have been attempted. The frozen APBS challenger
does not meet its physical acceptance criteria.** This is a measured result,
not the earlier pending or budget-stopped status. The baseline remains default.
This finding concerns this implementation, not every environmental model.

## What ran

- Six new native r2SCAN-3c/CPCM endpoints finished; six MBIS charge distributions
  passed charge closure and the prescribed actual-ESP quality check.
- The continuation reused these outputs and the one completed APBS state.
  It launched **102 independent charging solves with 102 workers**, with no
  compute budget or time limit. No additional high-level endpoints.
- Job **1198968** finished in **171 seconds**. There were 100 finite new
  charging outputs and two nonfinite outputs, both in the identity state's
  all-zero-charge environmental terms. Including the reused state, **17/18
  APBS states have complete numerical results; identity failed explicitly**.
- APBS printed `-NAN` for both zero-charge terms despite exit status zero.
  Their values remain unavailable. They were not replaced with invented zeros
  or reported as a passing identity test.

## Physical results

Tolerances were frozen before these outputs; none was adjusted afterward.

| Check | Observed movement, kcal/mol | Frozen tolerance | Result |
|---|---:|---:|---|
| Baseline contrast reproduction | ≤0.000029 | 0.5 | Pass |
| Common physical cavity across endpoints/partitions | Identical serialized cavity hashes | Identical | Pass |
| Repeat, including independent-process execution | 0.00000000000182 | 0.5 | Pass |
| Box extension, endpoint corrections | 0.0144–0.0171 | 0.5 | Pass |
| Rotation, qm33 endpoint corrections | 0.0609–0.2274 | 0.5 | Pass |
| Translation, qm33 endpoint corrections | 3.3653–4.2793 | 0.5 | **Fail** |
| Grid refinement, endpoint corrections | 2.8795–4.6474 | 0.5 | **Fail** |
| Analytic versus grid Coulomb cross term, primary | 1.2200–1.4126 | 0.5 | **Fail** |
| Analytic versus grid Coulomb cross term, refined | 0.7113–0.8318 | 0.5 | **Fail** |
| Corrected-score partition movement | **10.2341** | 2.0 | **Fail** |
| Identity | Two nonfinite APBS outputs | 0.01 | **Unavailable / solver failure** |

Some numerical error cancels between metals, but insufficiently: the La/Ca
correction changes by **0.6576** (qm33) and **0.8604** (qm36) under refinement,
and **0.9140** under the prescribed qm33 translation. The partition movement
is still **10.0313** on the refined grid. These are gauge-independent
differences; no new threshold or aquo reference was fit.

The energy components show substantial cancellation, not a unique diagnosed
cause. For the score contribution `DeltaU_Ca - DeltaU_La`:

| Partition | Direct core/environment term | Total reaction-field contrast | Subtracted core-reference reaction-field contrast | Net correction |
|---|---:|---:|---:|---:|
| qm33 | +58.6201 | −109.2558 | −67.9992 | +17.3635 |
| qm36 | −0.4452 | −111.2876 | −137.1264 | +25.3936 |

All entries are kcal/mol; net = direct + total reaction field − reference
reaction field. Environment-only reaction-field terms are identical between
metals within each partition and cancel. Raw APBS terms, mappings and receipts
are retained. The component changes do not establish which approximation is
responsible without a separate investigation.

## Cost and usefulness are separate judgments

**Numerical credibility:** fails the prescribed checks. **Predictive usefulness:**
not established; only consumed 1H4I development preparations were used, and
the failed checks prevent interpreting the correction as a validated score.
**Affordability:** actual parallel execution is short, but no matched production
baseline/challenger ratio has been established. There is no CPU/GPU spending
threshold stopping development.

| Job | Purpose | Wall seconds | Allocated core-seconds |
|---|---|---:|---:|
| 1198934 | Six quantum endpoints | 1491 | 95,424 |
| 1198939 | Parser-failed follow-on | 2 | 688 |
| 1198958 | ESP and serial APBS attempt | 166 | 57,104 |
| 1198968 | Parallel completion | **171** | **58,824** |
| Total | All attempts retained | — | **212,040** |

The parallel APBS processes used **6381.92 user + 1196.93 system CPU seconds**;
summed process wall time is 7612.57 seconds, not job wall time. Largest single
process RSS was 7,684,672 KiB; batch MaxRSS was 408,080,160 KiB. The earlier
one-CPU/full-node execution was wasteful and remains in the cost ledger.
Small preparation costs are recorded separately; prior uninstrumented work is
not silently assigned zero cost. No GPU allocation was used.

## Code, tests and remaining scope

Compute-budget stopping rules were removed from both development runners.
Original budget-bearing manifests remain historical records. The concurrent
runner checks frozen state hashes, solver/input/output provenance and cache
compatibility, records partial failures, and uses resource-aware concurrency.
Vectorized preparation reproduces original PQR bytes and Coulomb values to
better than 1e-10 kcal/mol on the real repeat case.

**23 software tests passed, none skipped**, including real MBIS integration,
split-input equivalence, unchanged state hashes, source-geometry repairs and
honest handling of the actual nonfinite identity outputs. This does **not**
mean 23 scientific acceptance checks passed: the physical failures above are
preserved and the reported environmental numerical score remains null.

Recommendation: **retain the baseline; do not promote this frozen challenger**.
A new numerical or boundary model would require a new version and its own
tests, without treating these consumed cases as blind. No further parameter
variants were launched to make this one pass.

The separate peptide-amide v3 repair remains prepared and verified on six
real sites; its energetic comparison has not run. Gradient mapping support
remains experimental; mechanical and entropy scores remain disabled as
`response_model_not_validated`. The fixed-core whole-chain environmental
preparation remains unsupported at terminal Lys. These are explicit limits,
not completed scientific validations.

Exact results: [COMPLETED_PILOT_RESULT.json](COMPLETED_PILOT_RESULT.json).
Paired ledger: [comparison_completed_checks/comparison.json](comparison_completed_checks/comparison.json).
Next read-only command: `cat diagnostics/affordable_challenger_20260915/COMPLETED_PILOT.md`.
