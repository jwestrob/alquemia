# 02_endpoint_recovery: Recover the two missing Hans endpoints

Read COMMON_RULES.md in full and pass it to every subagent.

## Start condition

01 accepted; Slurm account access; verified session and Slurm completion wakes

## Inputs

Use the reviewed same-state continuation manifest and its pinned implementation. Reuse four accepted Hans8DQ2/Mex8FNS endpoints. Read recovery_20261007/RESULT.json and ORBITAL_INVENTORY.json; interrupted GBWs are guesses, never energies.

## Work

Prepare an exact submission request for only Hans8FNR La_A/Dy_A. This Claude session submits and reconciles its own jobs once scientific and scheduler gates pass. Retain the cancelled PModel and TRAH attempts. The old September-startup cutoff scripts must not be blindly resubmitted. Verify actual current allocation, executable, geometry/state/basis, MORead and convergence. Reuse repaired_exchange_20260928/COMPARE.py with a new explicitly named continuation config. No automatic retries, altered tolerances or iteration-limit increase. A failure generates a concrete diagnosis and returns control.

## Completion criteria

Deliver accepted endpoint receipts/gradients or explicit unavailable rows, complete cost accounting, failure diagnosis where applicable, and a six-cell comparison retaining both Hans sources. End waiting_slurm when jobs are pending/running; resume only on actual events. Task completion does not mean discrimination success.

## Ownership

Own only this task’s worktree branch and diagnostics/claude_lanm_20261007/02_endpoint_recovery/; raw products under workspaces/claude_lanm_20261007/02_endpoint_recovery/. Request a scoped follow-up before editing shared executors/scoring. Root owns integration; this Claude session owns and may submit its assigned scientific jobs.
