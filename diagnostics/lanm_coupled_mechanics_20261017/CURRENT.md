# Coupled La/Dy mechanics is active

17 October 2026 UTC. Jacob activated the goal: “Implement the saved La/Dy mechanics plan. Work autonomously, validate the physics, and document the results.” PLAN.md defines the finite stages; GOAL.md remains the expanded scope. The initial planning pause is over.

Root recovered all six exact admitted restricted-search terminals in START_BINDINGS.json. This metadata binding made zero numerical or molecular calls. Existing LBT common22 and Mex common5 pools remain pinned for reuse. No completed calculation is restarted.

## Current ownership

- Root: source/configuration pins, integration, scientific submissions, receipts and shared status.
- Leibniz (`lanm_science_next`): LBT full-coordinate evaluator/optimizer and finite qualification adapter in this directory.
- Harvey (`lanm_resume_inventory`): LBT electronic/MM partition and interaction ledger under `lbt_partition/`.
- `mechanics_mex_adapter`: Mex arbitrary-coordinate evaluator and newly free-coordinate finite-difference checks under `mex/`.

Mex34-composite qualification12952 PASSED all ten replay/directional checks in307scheduler seconds on12CPUs/oneGPU. Both step sizes and every component passed; warm1MACE load. Search guard preflight12955 also passes. Read mex/QUALIFICATION_REPORT.md and immutable result pins in INITIAL_EXECUTION_RECORD.json. No Mex optimization launched before the LBT scout gate.

LBT partition recovery12953 PASSED4/4 preparations with actual MM XML exactly unchanged from failed12950. Baseline172/expanded210 electronic atoms include two mapped caps; both sources/elements preserve charge−2. Read lbt_partition/PREPARATION_REPORT.md. A new local-electronic/MM adapter is being implemented by mechanics_mex_adapter under lbt_hybrid/; no molecular calls yet.

LBT zero-call preflight12954 failed an old exact transformation-matrix replay. Saved-only diagnostic12956 proves every other lineage conjunct passes; matrix differences are1.11–1.94e−16, and archived/recomputed whitening residuals remain around1e−15 under the original1e−10 gate. Harvey owns a narrow versioned replay repair retaining the archived matrix as authoritative; historical code/results remain intact. Leibniz owns LBT/Mex paired-pool adapters. No new LBT energy yet. All seven initial jobs are terminal and recorded; no duplicate submissions.

Current observed scheduler policy is standard12CPU/job120CPU/user3days; standard_expand24CPU/job500CPU/user5days. Initial GPU qualification uses normal standard12CPU and oneGPU; no memory/time overrides or QoS bypass.

Next precise action: review LBT_CARTESIAN.py and its frozen spec/config contract when ready, bind actual parent12767 plus source/spec hashes in a fresh workspace, run zero-call real-fixture preflight onSlurm, then its declared32-cell derivative qualification. Only passed actual qualification enables the12-step/40-evaluation LBT scout. Mex prepares a separate34-composite qualification; root prevents duplicate chemistry and resource oversubscription.

All science uses Slurm. Warm workers, exact cache reuse, finite budgets, corrected/qualified parameter files and completion wakes are mandatory. The stopped clamped-core/internal-strain branch remains stopped. PQQ production is unchanged; the new goal remains unmet.
