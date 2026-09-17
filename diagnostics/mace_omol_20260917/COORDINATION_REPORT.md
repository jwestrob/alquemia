# Matched OMOL coordination descriptor

Numerical gate: True. Canonical gate: False. Non-PQQ primary/robustness: False/False.

Canonical calibration gap: 85.73797606923159 kcal/mol.

| Case | Coordination contrast, kcal/mol | Own calibration decision |
|---|---:|---|
| 1H4I | 39.74891207321713 | Ca-supported |
| 4MAE | 147.18627726485983 | La-supported |
| 1KB0 | 55.31731415524639 | inconclusive |
| GGR_extended | 67.69455916106968 | None |
| GGR_connected | 54.129509756126254 | None |
| ALPHA_1F6S | 54.43882936400771 | None |
| ALPHA_6IP9 | 51.51222152483938 | None |

| Non-PQQ comparison | Difference, kcal/mol | Pass |
|---|---:|---|
| ALPHA_1F6S minus GGR_extended | -13.255729797061967 | False |
| ALPHA_6IP9 minus GGR_extended | -16.182337636230294 | False |
| ALPHA_1F6S minus GGR_connected | 0.30931960788145574 | True |
| ALPHA_6IP9 minus GGR_connected | -2.617288231286871 | False |

These are consumed development comparisons. Both crystals overlap calibration accessions, two alpha structures are one biological group, and both GGR representations are another. PQQ composition confounding remains.
The finite-range coordination descriptor cancels geometry-independent terms but does not certify fragment ionic states or physical dissociation energies. No solvent, aquo offset, entropy, training or fitted correction is added. Production remains unchanged.

## Decision and measured cost

Reject this version as a replacement for the original OMOL descriptor: one PQQ
transfer loses support, both primary affinity comparisons remain wrong, and
representation robustness still fails. The larger calibration gap is not an
accuracy improvement. Preserve the original OMOL candidate and the production
baseline; continue the broader research goal.

All32 scores are available. Qualification1200803 ran10calls and passed19checks;
benchmark1200804 ran60newcalls and passed60isolated-metal checks, reusing four
qualified references and64original bound endpoints. No scientific failures,
newDFT, solvent calls or training. Total695GPU allocation-seconds,
11,120allocatedcore-seconds and769.676actualCPU-seconds. Peak allocatedGPU
memory3,099,778,560bytes. The median sum of archived four-endpoint evaluation
times is2.880806s; this excludes startup and is not end-to-end throughput.
The complete receipts include startup/collection; separate local reporting took
95.96wall seconds and84.70CPU seconds. Preparation/tests are not fully profiled.

Five real coordination software/integration tests pass(23.960s), including
actual four-endpoint algebra, completed numerical qualification, and rejection
of corrupted real geometry or incompatible bound caches. Earlier OMOL six tests,
readout four tests and legacy runner ten tests pass; two legacy kernel tests
were explicitly skipped because they require the isolated MACE environment.

An interpretation check against primary documentation is recorded separately in
[COORDINATION_NOTES.md](COORDINATION_NOTES.md). No test rules were changed after
seeing predictions. Raw workspace reports/results and the complete failed
candidate remain immutable and reproducible through COMMANDS.md.
