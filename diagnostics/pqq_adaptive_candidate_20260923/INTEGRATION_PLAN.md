# Fresh two-source integration, separately authorized

Parent authorization 2026-09-23: "let's actually validate fresh source-to-score
operation before advertising readiness ... exactly two existing pinned released-
source requests, one Ca reference and one La reference (prefer 1H4I/4MAE) ...
No search variants/retries ... This finite 2-source validation is now authorized."

Use 1H4I and 4MAE from the released source inventory. The first is the exact
existing source example; the second retains the inventory's raw source, selected
assembly/model/chain/metal/PQQ and every role. Omit only archived comparison pins
in the fresh request, matching the released standalone 1H4I example. This permits
an independent preparation; it does not change any preparation parameter.
Compare resulting states/coordinates and original scores to archives explicitly.

Per source: one fresh preparation, paired native energy/analytic force calls,
one bounded SLSQP200 proposal for each metal, two cross-native evaluations,
four origin plus eight proposal GFN2 calls (MaxIter500, native mixer,300K).
Total: four origin MACE, four bounded searches, four cross MACE, at most24 GFN2.
No DFT, folds, method variations, extra starts or automatic failed-cell retries.

Frozen minimal canonical25 reference and separate original bands; no refit.
Original release and all historical results remain unchanged. Both sources are
consumed controls, not new independent predictive validation.

Resource allocation: one H200,32CPUs,200000MiB; four 8-rank ORCA endpoints in
parallel, GPU requests via the existing warm native worker. The earlier two-row
recovery (without preparation/origin calls) cost49 GPU-allocation seconds plus
56CPU-allocation seconds. This integration adds measured preparation/origins;
report its actual allocation and all model calls independently.

Tests and exact request/plan hashes are recorded before submission. Any changed
source, unsupported geometry, missing force, failed optimizer or quantum cell
stays explicit. No molecular retry or class-directed repair.
