# Source preparation checks

Before any fit release, use at most four local worker threads, each numerical thread control1. No native executable, model engine, reference calculation or reserved validation output may be invoked/read.

1. Parse all new source and verify final direct config pins.
2. Exercise the pure search callback on an artificial four-dimensional quadratic, including out-of-order completion. Require complete polls, stable selection, four sensitivity columns, and no candidate reselection from sensitivity.
3. In a separate test module instance, lower only the in-memory synthetic search cap to9 and total cap17. Require a complete first poll and refusal of the next incomplete poll; the on-disk scientific limits remain233/241.
4. Exercise a callback that fails on a noninitial point. Require global stop, preserved failed slots and no subsequent poll; no molecular evaluator is involved.
5. Run actual saved-only preparation using the frozen config: admit the four baseline inputs/raw results, resolve the family, compute the symmetric baseline objective, and replay Water48 directional summaries. Require zero native invocations. This preparation records actual local wall/CPU costs and supplies the concrete run-release input.

These checks test software behavior and saved-data admission, not new parameter outcomes or reference accuracy.
