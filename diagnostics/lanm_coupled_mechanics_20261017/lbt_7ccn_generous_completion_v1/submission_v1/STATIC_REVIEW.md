# Independent static review

**PASS for the declared single run.** No material source, schema or budget blocker was found. This review permits the intended saved admission within the GPU allocation; it does not claim that the new adapter has executed successfully or that the model will converge.

## Scope and preservation

- The optimizer changes only its exact resume-admission block: checkpoint digest, old/new limits and explanatory messages. All text outside that block is byte-identical to optimizer `c5092f24`. Direction construction, rigid projection, history, trust, Armijo arithmetic, convergence and physical guards are unchanged.
- Admission is restricted to the exact step 172 checkpoint `14cd220d`, canonical digest `2d9a27465af93f883099a3e9d0562af1ffe98f41ad8d590587cb10dce3a06535`, and terminal POINT `67690a33`. The cumulative caps are 1,172 accepted steps / 2,172 evaluations: at most 1,000 additional accepted steps and 2,000 distinct evaluations from actual 172/172.
- The loader reuses pinned prior `prepare` and `run_gate`, including actual earlier preflights and decisions, with the original lazy evaluator. It preserves the original initial point/source, recipe, charge −2, physical La multiplicity 1, all 256 movable atoms, eight history pairs and saved trust 0.1 Å.
- All 193 archived point records / 189 coordinate keys and 172 prior requests are admitted and retained. The prior cache is a required subset; exact initial/old/current endpoint coordinate, gradient and native-component replay precedes search. No source rebasing, cache reset or repeated saved chemistry is introduced.
- Actual audit 15489 is hard-pinned and required complete, with 81/80 admissions, all 2,828 checks passing and zero calls/model loads/trials/optimizer steps. Its recipe, settings, initial/final checkpoint, source, helper and prior optimizer match the literal saved predecessor schema. `inherited['rationale']` is populated by the pinned previous loader and is valid.

## Execution and failure boundaries

Only the single 7CCN–La search is reachable. Saved admission, the zero-cost/no-model assertion and START receipt precede lazy native startup. The original evaluator and search enforce the finite limits; the result additionally checks each component's 2,000-call ceiling. No DFT, GFN, MM, second state or extra search path is added. The one-state denominator and unavailable failure record remain present.

The inspected launcher requests one qualified H200, 12 CPUs and the original Python environment, without memory/time overrides. OMP follows 12 CPUs; other BLAS thread settings remain 1. Root must use the declared `--qos=standard_expand` submission override: the preserved launcher header says `standard`. READY's `launcher_arguments` correctly omit the Python executable already supplied by the launcher. The output directory was absent at review.

Original force criteria and termination behavior remain binding. The new budget is a finite model-minimization attempt, with no automatic extension or intermediate reference. Audit success and earlier incremental DFT support do not qualify later force accuracy, electronic stability or the eventual physical minimum.

## Reviewed pins and method

- CONTINUE.py: `1fc6b168c0807a50223b2e5bb2b77837e19efd9d210d4c98f3b3a9e8785a400a`.
- Optimizer: `6635c0f062ec9f41cd5bc84e12606af2490bd34a5ed216843d1735ea7b017db2`.
- SPEC: `10dc171ca4c472079209fb4dfc0a59d860b1f1f2baff42d935155595438435ec`.
- CONFIG: `1bfaebc040879107ac72d5f7a3777591670c76857bd4bd9af57314fdd173b454`.
- READY: `83635e671ce50d260f2fa9c6f62c3df2aa68c03e5d989f8e3b8cb646b8a14693`.
- Science decision: `e8c4bda02f6714f452251b2c4edd83dbbdebfdbafa81b84cfdef75cf5574cc80`.
- Launcher: `e8a0f2b21922061d8c07f3333c3e57c741b469245ca2be9bb065c7fa2abb2bfd`.

All 16 direct unique CONFIG/READY pins were verified. Review used source reading, exact text comparison, AST parsing, JSON field/count/equality checks and file hashes only. No scientific imports, numerical replay, molecular calls, tests, submissions or production edits occurred. No separate duplicate preflight allocation is requested; the declared run must complete its saved admission before new calls.
