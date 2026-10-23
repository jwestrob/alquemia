# Four fixed coupled water-parameter probes

## Scientific question and fixed matrix

Do the two saved bounded coupled directions actually improve the original four-parameter water objective? The local Jacobian predicts improvement, but its sampled nonlinearity precludes treating that as an observed result. This is one finite scout, not a fit continuation, parameter-family expansion or transfer test.

Start from the unchanged pilot17567 selected vector, bound by RESULT SHA256 `524d093fb920304e254369097a96843408ab21597aa7898f5d1a3947c64fe6b6`. Bind saved local diagnosis `7943ce93d34ce7434d84e836e8e13639435564ab4c86e18a472eb8a25f505647`. Freeze exactly four vectors in this order: half and full of `uniformly_scaled_Gauss_Newton.normalized_delta`, then half and full of `bounded_linear_least_squares.normalized_delta`. Convert each displacement using the original four interval widths and add to the original selected physical vector. Every normalized component is at most0.025 and all original bounds hold. The unrestricted Gauss–Newton point is excluded.

Each vector uses the same four water cells, fixed epsilons, original objective/weights and separate La/Dy/paired/energy controls. No coordinates, charges, polarizabilities, water/protein terms, reference values or thresholds change. Reuse the original EVALUATE.py and CONFIG_EVALUATOR.json byte-for-byte. The original pilot/source/result/ledger and held seven-parameter family remain immutable.

## Admission, cache and finite execution

Saved-only preparation pins all source/config/reference inputs and these four numerical vectors before native execution. Re-admit the original selected model from its exact completed native cells and replay its score. Prepare a **fresh** cache and zero-count ledger. Seed it with symlinks used read-only to the142 completed original native-cell directories, checking source authority, input/output hashes, return codes and effective cache identity; never copy or continue the old ledger. The evaluator's completed-cache path performs no writes to those directories. Verify the original ledger/result hashes again after the run.

Run the four vectors once with a global limit of **32 new native program invocations** (four vectors ×four cells ×analyze/testgrad), in addition to the original evaluator's241-vector/1920-program ceilings. The thin runner checks the scout cap before every native program through a scoped guard on the unchanged evaluator instance. Native physical-input cache reuse remains mandatory. The root release binds both the original evaluator contract and this narrower four-vector/32-program scope. No generic optimizer or fallback is invoked.

Use a new four-CPU root launcher, one thread per native process and at most four simultaneous vector workers. The original eight-CPU pilot launcher remains unchanged. First native/admission failure stops new scheduling globally; preserve any already-running siblings, all four vector slots and their four-cell availability. A failure yields no new selection or retry.

## Selection and interpretation

Only after all four vectors complete, compare their unchanged objective values in the frozen order, with the original selected candidate first. Keep the lowest value and retain the earlier point on exact ties; if none improves, retain the original. Report all separate control outcomes for every point, not just scalar loss or the selected point. Choosing a smaller loss does not waive failed controls.

For every probe, compare actual and saved linear-predicted loss/control values and report the existing dimer CC/CP method-sensitivity diagnostics. Use the original derivative direction. No CC energy is substituted into the PBE0 objective. The selected point is a consumed development candidate, not a validated paired potential; reserved24 outcomes remain closed and no protein/bulk sampling is released.

Root alone submits after actual zero-call preparation and independent source review. Record source/input hashes, aggregate start/end child CPU, parent CPU, wall time, scheduler costs and exact new invocation counts. No additional probes, parameter steps, epsilon changes or automatic extension follow this scout.
