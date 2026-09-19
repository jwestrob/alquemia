# Completed PQQ utility comparison

**Masked MACE preserves the prepared reference task and is faster. An accuracy improvement has not been demonstrated.** Both fresh MACE runs return the correct 25/25 canonical calls and exactly reproduce the archived scores. The supported crystal transfers remain 2/2; 1KB0 remains unsupported, so coverage is 2/3 of the complete transfer set.

The cached workflow takes 129.00 seconds median versus 239.84 for DFT (1.859x); the original MACE workflow takes 154.72 seconds (1.550x). These are complete prepared-input-to-score timings on the recorded hardware, not inference-only timings.

**Recommendation:** retain the DFT production default. Keep frozen MACE available as an optional PQQ research scorer; do not claim a fully qualified replacement or broader accuracy gain. The strict combined gate fails solely because one DFT rerun differs from its archive by 0.12477 kcal/mol, exceeding the declared 0.01 reproduction tolerance while retaining its Ca classification. MACE's fidelity and the measured speed advantage are separate positive findings. See [the endpoint investigation](DFT_REPRODUCIBILITY.md). No scientific rerun or parameter change was used to remove this failure.

The generated report below preserves the original automated gate and recommendation verbatim. Its generic negative recommendation must be read with the specific DFT reproduction failure above. [Commands](COMMANDS.md) reproduce the complete comparison without new energy calculations.

---

# PQQ reference fidelity and measured MACE utility

Recommendation: **do_not_adopt_for_claimed_PQQ_speed_utility_on_this_evidence**. The production DFT default is unchanged.

## Prediction fidelity

| Reference class | Cases | DFT correct | Original MACE correct | Cached MACE correct |
|---|---:|---:|---:|---:|
| Ca | 14 | 14 | 14 | 14 |
| La | 11 | 10 | 11 | 11 |

Archived crystal transfers: DFT 3/3; MACE 2/2 scored, with 1KB0 unsupported (2/3 of the complete transfer set). Reference replay preserves performance; it does not demonstrate improved accuracy on unseen proteins.

## Full workflow timing

| Workflow | Median seconds / protein | Total case seconds | Fresh correct / 25 | Scores reproduced / 25 |
|---|---:|---:|---:|---:|
| DFT | 239.84 | 6806.86 | 24 | 24 |
| Original masked MACE | 154.72 | 3877.83 | 25 | 25 |
| Parser-cache masked MACE | 129.00 | 3232.25 | 25 | 25 |

Ratio of median times versus DFT: original **1.550x**, cached **1.859x**. The predeclared engineering target was 1.5x plus lower total time.

## What this establishes

- The frozen masked-MACE descriptor retains the labelled PQQ reference task; the cache changes only repeated coordinate parsing.
- Archived calibration: 25/25 for both methods. Consumed crystals: DFT 3/3; MACE 2/3, with 1KB0 unsupported.
- The benchmark uses each method's own frozen calibration. These are consumed references, not new independent biological validation.
- PQQ functional association is the target. Broad affinity discrimination, physiological occupancy and within-lanthanide preference are not established.

## Numerical boundary and coverage

Literal decisions and score-reproduction checks are reported separately. Small DFT rerun changes can cross a calibration extremum; the full table retains every inconclusive call. No threshold was moved. 1KB0 remains unscorable by whole-chain MACE, with no baseline substitution.

## Measured resources

| Job purpose | Allocation seconds | Allocated CPU-seconds | GPU allocation seconds |
|---|---:|---:|---:|
| DFT | 6808 | 217856 | 0 |
| MACE_original | 3879 | 62064 | 3879 |
| MACE_cached | 3233 | 51728 | 3233 |
| failed_DFT_startup | 3 | 96 | 0 |
| report_only_check | 107 | 107 | 0 |

DFT used 32 Slurm CPU units; each MACE run used one A5000 and 16 CPU units on the same host. Each main job requested 64474 MiB. CPU units are scheduler units, not an asserted count of physical cores. A GPU-second is not equated with a CPU-second or a monetary price.

Timing includes validation, loading, execution and reporting from existing prepared inputs. It excludes queue waiting and original folding/protonation. Each workflow has one fresh run per case; shared-host variability was not estimated through repetitions.

All 150 successful fresh endpoint receipts, four failed MPI startups, exact scheduler costs, local development receipts, and original/cached comparisons remain linked in result.json. No new scientific model, reference, threshold or production default was introduced.
