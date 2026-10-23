# Independent aqueous namespace/numerical review

**The exact prepared mapping and two-cell numerical test pass source/metadata review.** No native calls or numerical control reruns were made. Root still binds the release.

The mapping preserves BIO parameter bytes, frozen17605 water parameters and exact17612 aquo8 geometry/topology/native lineage. Only metal types/classes900/901 and their selected water-O class73/type402 pairs are added. Water H and salt pairs retain native mixing; no donor parameter is bound. All four prepared keys have the same mapping. The full aqueous boxes remain unrun. Common metal mass138.91 is explicit and does not enter this static potential comparison.

Native initprm.f sets Urey unit1, cubic0 and quartic0; no key or parameter overrides them. eurey1.f implements k(r−r0)^2 without a half factor. Therefore the declared BIO-minus-SI correction and its gradient sign are correct for k=−7.60 and the1.5337→1.5537Å equilibrium shift. Oxygen/metal correction is zero; each H-pair gradient is equal and opposite. The actual water-model difference remains visible.

The runner requests exactly2×(analyze+analytic testgrad), using saved SI outputs. The existing strict parser requires all25 indexed finite gradients, finite consistent energies, native version and no failure marker. Final admission requires energy remainder≤1e−4, maximum Cartesian gradient remainder≤1e−3, matching component coverage/counts, Urey remainder≤1e−4 and every other component difference≤1e−4. Missing cells cannot pass. Root's12 arithmetic/failure controls and83 author metadata controls were reviewed without repetition; actual zero-call preflight is bound to the reviewed configuration.

One nonblocking runtime observation: a native text failure with exit0 is detected by the joint parser after both programs, so testgrad may already run. Such output still cannot pass. No criteria change is requested.

This validates readiness for a finite mapping test, not physical water populations, salt pairing, periodic electrostatics or selectivity. Exact pins and saved metadata checks are in the accompanying workspace REVIEW.json.
