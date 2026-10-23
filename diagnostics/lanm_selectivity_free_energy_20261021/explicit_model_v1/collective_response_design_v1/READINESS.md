# Prepared exchange-polarization discriminator

**Source and saved preparation are complete; no new molecular calculation has run.** The proposed native term changes the anisotropic induced response of water O near a metal, with all original parameters retained. It can change both isolated-pair and collective energies, which the fixed matrix measures separately.

- Final source: `RUN.py`, SHA `9eb552e6ec65da01b3e62ac803740405ee7cc05ab45d1371294fdea5b3ca422d`.
- Scientific declaration: `SPEC.md`, SHA `7e5f947165ef94693538db7a912249ab303b2c926c85361e8290df73c74d71fc`.
- Actual saved preparation: matching workspace `preparation_v3/RESULT.json`, SHA `880cfab900e427f2473fb527a17c407356ad4ecd7d4c57ece2498bae5db74912`; 19.935 s wall, 18.442 CPU-s, one thread, zero native calls.
- Focused controls: `source_controls_v1/RESULT.json`, SHA `3d7aa596b9386d44252cf46f5ac6cb238dd57e1b8a40b9bc5c820e9c5ef8ea14`; 229 checks pass, 1.759 s wall, 1.343 CPU-s, zero native calls. Includes exact inventory, unchanged input prefix, positive off-diagonal tensors, malformed induction summaries, hard limits, missing/nonfinite controls and actual worker ordering under a synthetic analyze failure.
- `READY.json` binds final config, manifest, four-CPU launcher, release template and exact command arguments. Earlier preparation/source revisions remain preserved.

The root release must bind `CONFIG_RUN_TEMPLATE.json` to a concrete release with status `approved_collective_response_exchpol_diagnostic`; the template is unbound. Launch order is numerical invariance controls, four native Cartesian derivative checks, then the remaining fixed matrix. Maximum cost is **116 single-thread native programs plus 636 internal finite-difference energies**, with four workers. No fit, candidate selection or automatic extension is implemented.

`admit_result(config_dict, manifest_dict, result_dict)` replays saved raw outputs, receipts, original scores and C/P/B without native calls. Output `summaries.q01` and `summaries.q05` contain the two complete ten-cell scores, component-resolved cycles and baseline changes. `summaries.scale_dependence` compares their normalized C/P/B responses; it does not select an amplitude.

The native CPU binary contains the EXCHPOL setup, tensor, damping and derivative symbols (read-only `BINARY_SYMBOL_EVIDENCE.json`). Actual force consistency and induction convergence still require the declared run. Site-tensor positivity is not a proof of full-system induction stability. Whole-protein zero-alpha handling, BIO mapping, periodic/GPU parity and physical accuracy remain unqualified. The common-basis quantum decomposition can change the interpretation of the remainder; it cannot erase the original full coordination error or uniquely identify it as missing polarization.
