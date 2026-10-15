# Native protein-exclusion scout: mainly a common shift

All 12 declared native cavity updates completed. The two origin sources receive
almost the same metal contrast. This does not explain their different classifier
behavior and does not justify adding this component to routine scoring.

| Source | Ca minus La transfer (kcal/mol) | Shape part | Cap-volume control |
|---|---:|---:|---:|
| sample0 control | 5.031545227 | 0.544100961 | 0.003540124 |
| sample3 ambiguous | 4.936869274 | 0.542446286 | 0.003588317 |

Raw sample3-minus-sample0 contrast: -0.094675953 kcal/mol.

## Numerical qualification

Archived replay is exact; translation, rotation and permutation pass the frozen
1e-9 Eh electrostatic tolerance. All native block/repeat and endpoint precision
checks pass. Both paired contrasts pass their precision flag. The aggregated
source-contrast flag is **false**: its conservative summed charge bound exceeds
the frozen 1e-6 Eh ceiling. Keep the raw contrast with this qualification; no
tolerance change or rerun is justified to rescue that flag.

## Scope and next decision

Close this origin-only exclusion hypothesis as a remedy for the source gap.
The test holds archived vacuum populations fixed and excludes permanent protein
electrostatics. Synthetic caps remain native-radius centers. The cap-volume
control is small but does not show that charged caps are negligible; their
real/cap/cross blocks remain in RESULT.json. Native hydrogen-bond burial changes
are reported separately and excluded from T. This is not total solvation, an
affinity, a calibrated correction, or a test of accommodation-dependent burial.

## Actual execution

Worker4562: 136 scheduler seconds, four allocated CPUs (544 allocated CPU-s),
133.610 measured TotalCPU seconds; scheduler peak RSS 7640020 KiB. The serial
engine used four shared slots for default memory. Runner/child timing boundaries
are separately retained in RESULT.json. Preparation4539 used7s/1CPU; build4543
used4s/1CPU. No new SCC, MACE, DFT, optimization or folding. Production unchanged.
