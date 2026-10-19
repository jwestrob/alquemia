# Independent static review: development CPU continuation

19 October 2026. **PASS; no static blocker for the declared four-call CPU stage.** Review was source, saved metadata and file-hash inspection only, without scientific imports, numerical execution or submissions.

Frozen implementation `46c8d49c8a0fd25adc52853ce8a4bd218c42a331070a05cf7a72a2c864df7bdf`, SPEC_CPU `9bfe218684924c92ef7590f1b09521c419230a16f61ec952770eaf4f8f0b5fa2`, CONFIG_CPU `92fb3225b943e66b81b5e85f8ac7735cee95533d1a0f9999f165ffa547fe0504`, and RUN_CPU `44f7f28996f1c84fcaee9db308965d7197d2f964c1965b7d702f073eae37b9ce` were checked. All exact path/SHA config records verify. BASE_DELTA remains unchanged at `f760cd79…`.

## Actual parent binding

The config binds passed zero-call preflight 14799 (`40fa2210…`) and admitted one-call scout 14800 (`71cd1a1f…`). The stored scout declares native and cost gates passed, wall time 14.724745035171509 seconds, one whole-GFN attempt and no model load. The continuation checks these actual gates, replays the original saved preparation and all eight signatures, then re-admits the scout's raw native receipt `f3712a77…`. It matches the scout's physical source, IDs, XYZ hash and recorded native time. No scout evaluation is repeated.

## Calls, state and failure behavior

The fixed sequence is three whole-GFN calls (7CCO after-La, before-Dy, after-Dy) and one isolated regional-GFN call (7CCO after-Dy). No other cell/scope is allowed, request IDs cannot duplicate within the invocation, and the maximum new count is four. Successful completion records cumulative counts of five evaluations including the scout: four whole, one regional and zero MACE. The candidate-wide ten-call limit remains inherited; later GPU and challenge work require separate source and actual gates.

Each request stages the exact declared XYZ and complete parameter text in a fresh directory. The full system has 256 atoms/target 255; the capped region has 172/target 0. Original charges, physical multiplicity metadata, effective GFN UHF 0, 750/538 valence counts, default solver settings and no-field/no-solvent Hamiltonian remain unchanged. The original native admission and finite gradient checks are used.

All eight slots now exist before preparation. Native attempts are counted before calling the backend. A native or validation exception explicitly marks its result unavailable and exits the sequential loop; later calls remain not requested. The two review findings concerning early denominator retention and post-admission exception status were fixed before freezing. No fallback, retry, damping change or automatic next stage exists.

The launcher requests one CPU, no GPU and no requeue with default scheduler time/memory and one native thread. The Slurm guard precedes scientific helper use. Saved preparation/certificate replay remains free of molecular calls.

## Limits of this pass

This stage does not assemble candidate/base work or establish the both-metal development gate. The single missing regional MACE call, final assembly contract, cumulative six-call development total, saved comparison and any conditional four-call challenge still require their own reviewed source and actual admission. No derivative, physical-state, force-accuracy or optimization qualification follows from this review.
