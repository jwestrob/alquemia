# Continue the planned checks without compute budgets or time limits

Jacob's explicit instruction in this conversation:

> ok. no finite compute budgets. no time limits. proceed.

This supersedes earlier development compute-budget stopping rules, including
the allocated core-second and endpoint budget checks. Costs remain measured
and reported, including failed attempts and unused allocated capacity. No
walltime or subprocess timeout is introduced. Actual available CPU/memory
determines concurrency; it does not stop the task list after a cost threshold.

Execute the existing 18-state physical-check schedule. Reuse all six validated
MBIS/ESP records and the one completed APBS state. The remaining 17 states
contain 102 independent charging solves; run these concurrently within one
existing-policy Slurm CPU allocation. Each solve retains its original molecule,
grid, dielectric and boundary settings; splitting the six independent blocks
into separate processes changes execution only. Identity and repeat checks
remain actual calculations, not substituted cached values.

No additional high-level endpoints are needed for this continuation. All
scientific settings, geometry, waters, charge representations, acceptance
tolerances, case selection and the baseline default remain unchanged. The
new preparation verifies every state against the hashes in the earlier frozen
numerical schedule. No additional biological validation or model variants are
silently added. Finish the specified checks and report their actual outcomes.

The original queued/failed/budget-stopped reports are immutable history; this
instruction authorizes continuing instead of treating their cost threshold as
a blocker. Scientific failures still remain failures and unavailable scores
remain unavailable.
