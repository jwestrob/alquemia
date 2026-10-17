# Coupled La/Dy mechanics is active

17 October 2026 UTC. Jacob activated “Implement the saved La/Dy mechanics plan. Work autonomously, validate the physics, and document the results.” PLAN.md is the finite protocol. Root owns all submissions and integration; PQQ production remains unchanged.

## Actual results

- Mex qualification12952 passes all34 composite evaluations and all component/total derivative checks, including backbone, boundary, water and spectator movement. Search preflight12955 passes.
- LBT partition12953 passes all four preparations, with native MM XML unchanged from failed12950. Baseline172/expanded210 electronic atoms each include two mapped caps.
- LBT original preflight12954 failed only old exact-transform replay. Saved diagnostic12956 and reviewed V2 preserve the archived transform and all physical/scientific gates. Preflight12957 and all32 derivative cells12958 pass.
- Old LBT44-cell pool replay12959 passes with zero calls and an explicit required-missing-cell test.
- LBT scout12960 passes:12 accepted steps/12 new composites; energy−63.9031 kcal/mol and max force83.2803→12.7867 kcal/mol/Å. Physical/state checks pass, but it is not converged. Read LBT_QUALIFICATION_AND_SCOUT.md.

## Owned execution

Root submitted the six planned searches after reviewing actual scout admission and cost. LBT7CCO/La12961 completed/admitted80 steps, maxforce2.14091/RMS0.921311 and energy−98.5724 kcal/mol from the prior restricted terminal; still nonconverged. The other five initial jobs12962–12966 failed before molecular calls because their node lacked a pinned historical compiler path. Configuration-identical host-only recoveries use the qualified H200 node:

| System/state | Job |
| --- | --- |
| LBT7CCO La, exact scout continuation | 12961 |
| LBT7CCO Dy | 12967 (pre-call failure12962 retained) |
| LBT7CCN La | 12968 (pre-call failure12963 retained) |
| LBT7CCN Dy | 12969 (pre-call failure12964 retained) |
| Mex La | 12970 (pre-call failure12965 retained) |
| Mex Dy | 12971 (pre-call failure12966 retained) |

Zero-call LBT hybrid preflight12972 and Mex common5 baseline replay12973 are also submitted with completion wakes. Both source adapters passed independent review. Queued jobs remain submitted under normal scheduler policy; no priority changes.

All searches use12CPUs/oneGPU and the frozen80-step/200-evaluation limits. Scout work counts toward its continuation. Actual submission commands, pins and detached completion watchers are under each workspaces/lanm_coupled_mechanics_20261017/{lbt,mex}_run_*_v1/SUBMISSION.json. Inspect actual terminal results before any pool evaluation; never duplicate submissions.

## Ownership and next commands

- Root: source/config pins, scientific submissions, result integration, shared status.
- Leibniz: LBT/Mex source-local paired pools, preserving complete checkpoint caches and failed cells.
- Harvey: independent pool/hybrid review and bounded embedded-DFT derivative feasibility.
- mechanics_mex_adapter: distinct LBT local-electronic/MM model under lbt_hybrid/; zero-call preflight and nine-cell scout precede its eight-origin matrix.

Next: collect the six actual terminal searches; bind their results and run reviewed common-pool cross-scoring (at most four LBT/two Mex missing opposite-metal cells). Continue independent local-electronic/MM preflight and embedded-reference preparation. No reference launch until exact inputs, force-admission policy and measured scout are declared.

Current observed standard policy is12CPU/job120CPU/user3days, with standard_expand24CPU/job500CPU/user5days available. All numerical work uses Slurm, allocation-derived resources, unchanged parameter files and actual completion wakes. No priority changes, memory/time overrides, source-position springs or old clamped-core diagnosis. Goal remains active and unmet.
