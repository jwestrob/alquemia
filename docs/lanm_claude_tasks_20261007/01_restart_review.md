# 01_restart_review: Independent restart review

Read COMMON_RULES.md in full and pass it to every subagent.

## Start condition

none; Claude authentication required; review session can run while Slurm is unavailable

## Inputs

Read recovery_20261007/{REPORT,CONTINUATION_PLAN}.md, scripts/metal_environment_interrupted_seed.py, scripts/metal_environment_lady_embedded_prepare.py, scripts/metal_environment_lady_embedded.py and their real-fixture tests. Canonical prepared manifest: workspaces/metal_environment_response_20260926/lady_FNR_interrupted_continuation_20261007_v1/manifest.json. Baseline implementation commit7a2b0ce; verify current state rather than resetting to it.

## Work

Review the existing continuation, not a rewrite. Verify unchanged target geometry/field/state/basis; cancellation and GBW provenance; explicit MORead/TRAH execution evidence; no PModel fallback; fresh-convergence/gradient admission; resource mapping and failure visibility. Delegate code and scientific-state review independently. Tests, including parser tests, execute via Slurm. Propose narrowly scoped fixes for concrete defects; do not modify the shared implementation until the integrator assigns a fix task.

## Completion criteria

Deliver REVIEW.md and FINDINGS.json with severity, exact file/line, real reproducer where available, existing-test coverage and acceptance recommendation. Complete when every critical defect is fixed in a reviewed follow-up or explicitly blocks execution. No molecular evaluation or fresh orbitals.

## Ownership

Own only this task’s worktree branch and diagnostics/claude_lanm_20261007/01_restart_review/; raw products under workspaces/claude_lanm_20261007/01_restart_review/. Request a scoped follow-up before editing shared executors/scoring. Root owns integration; this Claude session owns and may submit its assigned scientific jobs.
