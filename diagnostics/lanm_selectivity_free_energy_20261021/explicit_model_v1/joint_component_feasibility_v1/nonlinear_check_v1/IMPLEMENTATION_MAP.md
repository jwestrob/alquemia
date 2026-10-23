# One finite native candidate batch

Source-only map, 22 October 2026. No adapter, evaluation, build or submission was executed. Root must first accept the saved prediction and freeze the origin, direction, fractions, exact vectors and limits. This map does not release them.

## Recommended containment

Use a small, separately pinned `FiniteBatchEvaluator` subclass of the qualified evaluator, with a new batch configuration/release and fresh ledger. Keep the original scientific evaluator configuration and prepared manifest as immutable parent authorities. Retain 12 cell workers; the shortest launch path is the existing one-task, 12-CPU Slurm wrapper. Do not invoke `REFINE.py`, change module globals, supply a fictitious fit release or forge Slurm environment variables.

The subclass is an explicit new execution policy, not an assertion that the old policy admits this batch. Its source, review, release and executor identity must accompany every new result/receipt. The inherited `implementation` used in scientific cache identities continues to identify the unchanged rendering/parsing/scoring implementation; it must not be presented as the complete new launcher identity.

## Exact reusable authorities

Paths below are relative to the repository. Let `D=diagnostics/lanm_selectivity_free_energy_20261021`, `W=workspaces/lanm_selectivity_free_energy_20261021`, and `J=explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1`.

| Authority | Path | SHA-256 prefix |
|---|---|---|
| Qualified evaluator | `D/J/evaluator_v1/EVALUATE.py` | `51e0ba135ae67f41` |
| Scientific configuration, already bound to actual binaries | `W/J/fit_v1/CONFIG_EVALUATOR.json` | `8a2863f09be68877` |
| Literal 46-cell manifest, baseline cells and 48 static terms | `W/J/evaluator_v1/preparation_v4/MANIFEST.json` | `8f29dde62c7462cf` |
| Parameter box / physical coefficient mapping | `D/J/PARAMETER_DOMAIN.json` | `a9e53eb59749b459` |
| Four-block objective | `D/J/objective_v1/OBJECTIVE.py` | `af7d26f6f6d7f0de` |
| Native source manifest | `W/J/native_v1/source_v1/MANIFEST.json` | `82886b766d192b1f` |
| Completed build | `W/J/native_v1/build_job_v1/execution/RESULT.json` | `95e1cea12ac779646` |
| Analyze / testgrad | `W/J/native_v1/build_job_v1/execution/build/{analyze,testgrad}` | `24a18bdda2998ed2` / `91b7bd20b8fa84ab` |
| Actual independent native qualification | `W/J/qualification_v1/independent_actual_review_v1/REVIEW.json` | `7be6786d93887537` |
| Original serial native invoke/parser | `D/explicit_model_v1/RUN_PANEL.py` | `4189602145853bc3` |
| Existing Python 3.12 / 12-CPU wrapper | `D/J/RUN_PY312.sbatch` | `05f6a45f31f9a8e8` |

The full parent configuration carries exact hashes for all helpers, targets, native assets and baseline admissions. Bind that configuration and verify its full records; prefixes here are navigation aids. The historical `fit_protocol` remains parent provenance, while a new batch protocol supplies current execution authority.

## Reuse and narrow policy changes

| Source location in `EVALUATE.py` | Reuse / required delta |
|---|---|
| `checked`, `effective`, `render` (59–100) | Reuse unchanged. Twelve box values map to physical a/b/c/d and R/epsilon; only the declared metal–O override and response record change. POLPAIR, all other parameters and literal geometries remain fixed. |
| `local_readback` (102–148) and `_parse` (283–312) | Reuse unchanged for new cells. Preserve dual-field moment, induction, spectrum, kappa, exported A, solve-residual and force admission. Four susceptibilities retain the admitted analytic linear-model evaluation convention and frozen quantum h=.001 targets. |
| `__init__` (212–228) | Retain source/manifest/helper/target checks and fresh-cache initialization. Its 673/24000/12 constants are the historical parent ceiling, not the batch budget. Store a separate batch fingerprint and bind the active smaller budget before any request. Never initialize against the old fit cache directory. |
| `_gate` (254–272) | Override explicitly. Copy the build/binary/qualification checks intact; replace only old fit status/24-step release/cap/cache checks with exact batch release, vector allow-list, new executor/configuration/protocol/review pins and actual launcher envelope. Preserve global stop and single-thread checks. |
| `_invoke` (274–281) | A narrow override is needed to enforce the smaller native-attempt ceiling under the same lock immediately before invocation. Preserve the bounded semaphore, failure stop, native helper, argv and receipt behavior. Keep request counts separate from invocation attempts and distinguish pre-call failures. Do not hold the ledger lock throughout the subprocess. |
| `_cell`, `_saved_inputs` (314–366) | Reuse existing behavior for new Slurm cells: analyze is admitted before testgrad starts; exact input bytes, two export names, receipts and no-retry failure semantics remain. Login use needs a versioned `_cell` provenance change because it currently indexes `os.environ['SLURM_JOB_ID']`. |
| `_assemble` (368–378) | Reuse unchanged: all 46 cells, 10 E/F fixtures, 34 pair moments, 4 C moments, 4 susceptibility aliases, 48 static terms and unchanged 531-entry score. Add paired P/B energy diagnostics afterwards using the admitted literal coefficient dictionaries; do not silently append them to the objective. |
| `evaluate`, `snapshot` (380–421) | Wrap evaluation with exact vector allow-list and batch-vector cap before delegating. Evaluate the frozen vectors in declared order, without choosing alternatives or generating probes. Snapshot must distinguish historical parent caps from active batch caps, current calls and inherited costs. |

Required native qualification remains `passed_local_response_native_qualification`, exact `source_manifest`, `build_receipt`, `executables`, verified reviewer implementation/result, and all four true flags: `zero_baseline_parity_passed`, `force_checks_passed`, `moment_response_checks_passed`, `export_checks_passed`. The completed build must retain its manifest/executable joins. A root assertion alone cannot replace this independent receipt.

## Cache reuse without rerunning the fit

The manifest already contains the zero-response baseline and all static water/shell/ion records. `_cell` automatically reuses a baseline cell only when that metal's complete effective six-value tuple equals its initial tuple. Static terms never require native calls. Baseline reuse is therefore automatic even inside a vector changing the other metal.

For a nonzero saved origin or unchanged-metal reuse, import only explicitly pinned, already admitted cell records into the new in-memory cache after exact `cache_identity` comparison. Identity includes effective values, geometry, parameter bytes, rendered-key hash, scientific evaluator/helper pins, source manifest and executables. Bind the frozen vector artifact, its original per-cell producer/result/input pins, completed fit admission and original cost. Label reuse explicitly; never rewrite a historical cell, attach a new release to its old receipt, or increment new-call cost for it. A complete saved origin can be reported directly rather than evaluated again.

Do not point the new evaluator at `fit_v1/run_v2/execution/cache`: its constructor can load that ledger and subsequent requests mutate it. Do not copy its ledger into a new scope. A mismatched/failed cache record stops this batch; it must not silently trigger a duplicate native call.

There is one known replay trap: the original disk-cache branch compares all parsed fields with exact Python equality. Previously admitted export diagnostics and scalar susceptibility can differ in host-side LAPACK rounding. The completed fit review retained exact raw arrays, inputs, energies, moments, forces and stored score/chi bytes while allowing only the already qualified derived arithmetic. Prefer reuse through that completed immutable admission. If an explicit raw replay is required, use the narrowly reviewed `continuation_v2/CONTINUE.py` comparison helpers and unchanged scout arithmetic bounds; retain both numerical gate sets and every difference. Do not relax `_parse`, replace stored objective values or re-audit the full fit. Relevant authorities are `W/J/fit_v1/run_v2/ROOT_FINAL_REVIEW.json` (`1fe38e62`) and the completed continuation review (`8658d5af`); its source is `D/J/fit_v1/independent_actual_review_v1/continuation_v2/CONTINUE.py` (`4506a258`).

## Finite inventory and launcher

Each completely changed two-metal vector requires at most **46 analyze EM + 10 analytic testgrad = 56 native programs**, with four analyze cells producing two matrix snapshots each. An unchanged-metal cache match reduces this to at most 28. All 46 objective rows still appear in every complete result. If root freezes N new vectors/fractions, the safe ceiling is 56N, lowered only by an explicit exact cache inventory. No extra derivative probes, field calls, six omitted regression trimers or quantum references belong to this batch. No requalification build is needed when the exact native binaries and kernels remain unchanged.

**Slurm path:** reuse the existing wrapper interface with the new runner, `run --config <new batch config> --out <fresh directory>`, one task and 12 CPUs, standard_expand, default memory/time, no GPU. Bind the actual wrapper/interpreter and execution allocation. Scheduler availability was not queried in this source task.

**Login alternative:** requires a separately declared and reviewed login execution envelope. The native `invoke` helper itself has no Slurm requirement; its standalone panel main does, and is not used here. Preserve the helper's 300-second per-program timeout. Replace only the subclass's resource gate and the `_cell` scheduler-only metadata with truthful host/PID/run-ID/`job_id:null` records. Never fabricate `SLURM_*` values. Retaining 12 workers avoids changing inherited executor/semaphore behavior and leaves room within the authorized 16 aggregate login-thread budget; root must reserve the controller's numerical work and other concurrent agents explicitly. Sixteen workers would require further executor/semaphore changes and offers no containment advantage. Each native/numerical worker keeps all five inherited thread variables at 1 (also BLIS=1 if present).

## Required batch binding and small preparation checks

The new release must bind parent configuration/manifest, new executor and source review, accepted prediction/result, exact origin artifact, explicit fractions and full vector bytes/hashes, expected row IDs, domain/targets, native qualification/build/binaries, launcher/interpreter, active vector/program/worker caps, fresh output/cache path and `reserved_outcomes_authorized=false`. Any failure stops the batch; retain complete/failed/not-attempted rows and actual costs. No 24-step search or automatic continuation is implied.

Only focused source/saved controls are needed before release: unchanged rendering/parsing functions; exact baseline/static and admitted-origin cache joins with zero calls; rejection of an unlisted/out-of-box vector, changed native qualification and wrong launcher; a mocked invocation ceiling and first-failure stop including queued cells; and exact scoring/decomposition of the saved origin. These controls do not rerun the old qualification suite or fit. Numerical success of the eventual new point, the 25 E/F gates, paired P/B changes and all old response blocks must be reported separately from any scientific decision.
