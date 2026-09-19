# Running checkpoint

2026-09-19T22:17Z: first six-endpoint pilot is **running**, Slurm1202429,
node-64-768g-15,64CPU/256GiB, four concurrent16-rank endpoints.
Elapsed1h14m. **No final CC energy or comparison is available yet.** All four
initial HF references converged. Alpha-Ca localMP2 guesses converged and are
performing the full CC integral transformation: prepared-alpha25 batches,
first two165.389/250.102s; original-alpha30 batches. GGR solvent containers
completed and PNO construction is underway. Two alpha-La tasks remain queued
inside the runner. Actual correlation output confirms CCSD, triplesON, PTES,
TightPNO1e-7 and correct core/electron accounting. GGR has158 correlated
electrons in both endpoints (Ca218−60; La208−50 plus46 ECP electrons).

Parent requested a clear unfinished checkpoint and yielding instead of frequent
polling during known multi-hour calculations. **The scientific assignment is
unfinished; recover final results or failures later.** Only the existing six
approved calculations continue autonomously. No additional model/expansion.

**Update2026-09-19T22:25Z: autonomous final collection installed at parent's
explicit request.** Detached login-host PID1633721, PPID1, polls only job1202429
once per60s. A terminal accounting state, empty queue and free execution lock
are all required before finalization. It collects actual outputs/receipts/cost,
retains partial failures, writes a report and final vault note, emails Jacob
only final results or a material failure, and removes only validated-successful
own temporary matrices. No new scientific tasks, reruns or allocations.

Receipt: `workspaces/electronic_accuracy_20260919/autocollect_v1/installation_receipt.json`.
Pinned code/config and real running-output preview are beside it. Six tests
pass; preview returned0/6, empty correlated contrasts, no email and no cleanup.
The new final-CC parser must still meet its strict checks on actual completed
outputs; a parsing/receipt/method failure is reported explicitly, without
inventing an energy. No placeholder successful CC fixture.

Automatic final artifacts will be under `autocollect_v1/final/`: `collection.json`,
`accounting.json`, `REPORT.md`, `scratch_cleanup.json`, `email_receipt.json`,
`DONE.json`; collection failures additionally receive `collection_error.json`.
The report is also written to this diagnostic's `REPORT_AUTOCOLLECTED.md`.
Final vault note: `2026-09-19_laca-independent-electronic-reference-final.md`.
The monitor exits after finalization. Its own lock prevents duplicate monitors;
existing email receipts prevent sending the same final notification again.

Actual transient storage exceeded513GiB while four tasks were active; the
filesystem had130TiB available. Peak observed host RSS~147GiB. After successful
completion, retain outputs/receipts/inputs/GBW and useful scientific artifacts;
remove only this pilot's then-unneeded `endpoint.runtime.*.tmp[.rank]` matrices,
recording the cleanup. **Never delete active scratch.** Final wall/core costs
and high-water storage remain uncollected. This is a one-time reference
diagnostic, not an affordable routine scorer.

Manifest: `workspaces/electronic_accuracy_20260919/direct_v1/manifest.json`
SHA256 `bfa67a6c204d11f5309409b9dd292a27c66ab14c7bf6c30ad1346bf847c8a19f`.
Original producer preserved in `direct_v1/implementation/`.
Six real-archive/source/algebra/incomplete-collection tests pass. Completed new scientific tests
remain unavailable until endpoints finish; no synthetic success fixture.

Parent owns no electronic files; this agent owns only new
`scripts/electronic_accuracy_*.py`, `tests/test_electronic_accuracy_pilot.py`
and this experiment's diagnostics/workspace. Prior response-probe files stay
untouched. No commit/default change. Parent commits response work separately.

Next, verify actual completed CC component parsing/core accounting, collect
all six, report contrasts and measured allocation. PQQ1H4I/4MAE same-method
guardrail is declared but **not submitted** pending measured full throughput.
GGR2FW0's narrow native prepared-alpha margin is a possible subsequent
independently versioned method-uncertainty test, not another already-running
calculation.

Owned files: `scripts/electronic_accuracy_pilot.py`,
`scripts/electronic_accuracy_collect.py`, `scripts/electronic_accuracy_autocollect.py`,
`tests/test_electronic_accuracy_pilot.py`,
and `diagnostics/electronic_accuracy_20260919/{PLAN.md,AUDIT.md,CURRENT.md,
LEGACY_B97_AUDIT.json,run_pilot.sbatch}`. No commit made. Parent may reactivate
this agent after job1202429 finishes for the final collection/report/vault update.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/electronic_accuracy_collect.py --manifest workspaces/electronic_accuracy_20260919/direct_v1/manifest.json --output workspaces/electronic_accuracy_20260919/direct_v1/collection_1202429.json
```

Use a fresh collection filename if collecting an intermediate checkpoint;
the writer never overwrites. The collector distinguishes incomplete/failed
tasks and has no baseline fallback. Freeze results only after all receipts are
available and inspect actual triples/CC energy labels before trusting the new
parser. Current parser already rejects a real native DFT output as a CC result.
