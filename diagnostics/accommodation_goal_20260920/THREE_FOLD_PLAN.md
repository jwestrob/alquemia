# Developmental three-fold robustness addendum

## Authorization and motivation

After the original frozen fold comparison completed, root explicitly authorized
this separate zero-calculation analysis under Jacob's overnight goal: actual PLM
inputs have three folds, so test whether the four-fold reference result plausibly
supports an efficient three-fold descriptor. The original La4/Ca5/balanced results
remain unchanged. This is already-consumed development data, not blind evaluation
or validation of the three actual PLM samples.

## Rule recorded before computing triples

Use the completed `folds_v1/comparison_v1.json` only. For each of the25 reference
proteins, enumerate **every one of the four possible three-member subsets** of its
four noncanonical La-conditioned source samples. Do not select a favorable triple
or remove a failed sample according to its score. Analyze both original core and
complete context, native MACE and the frozen composite separately. This produces
100 triples per representation, each with two method decisions (200 triple rows,
400 decisions). No new source, geometry, molecular calculation or coefficient.

For each method/triple require all three original contrast values; otherwise the
triple is unavailable. Use their ordinary median and the same exact existing
representation/method bands. Do not fit, widen or recalibrate a threshold.

Report all four triples per protein, missing members, each median/decision,
correct/wrong/inconclusive/unavailable counts, and explicit per-protein flags.
`all_four_correct` requires all four triples to be available and correct.
For a compact worst-outcome label use wrong first, then unavailable, then
inconclusive, then correct. Keep all individual flags/counts alongside that label
so the ranking cannot hide missing coverage. Also report the worst scored outcome
separately from missing coverage. Denominators remain100triples and25proteins per
representation/method; structural subsets are correlated, not independent biology.

No promotion or PLM rescore follows from this addendum. It tests sensitivity to
which existing La-conditioned folds enter an inexpensive descriptor, not thermal
populations, quantitative affinity or new protein-class generalization.
