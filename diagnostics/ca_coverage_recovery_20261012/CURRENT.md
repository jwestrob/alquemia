# Ca reference coverage recovery — complete, 12 October 2026

Read [DELIVERY.md](DELIVERY.md), [RESULT.json](RESULT.json) and
[SOURCES.tsv](SOURCES.tsv). All 18 originally unavailable accommodated sources
now have separate recovery results. Final collector 3049 completed and root
acknowledged the actual event.

| Same 125 sources | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Static origin | 117 | 1 | 7 | 0 |
| Accommodated operational/math | 119 | 3 | 3 | 0 |

Original primary125 remains 105/1/1/18. No new fit, source selection, ensemble
policy or production promotion. These are correlated model sources.

## Completed work

- MxaF interrupted samples 1/4: both correct Ca; jobs 2869/2873.
- MMOL_2048: one existing Cys128 hydrogen mapping repair, no heavy-atom change;
  jobs 2864/2870/2875/2880. Static La-correct, accommodated inconclusive.
- Fifteen alternative two-Ca AF3 sources all have unique PQQ Ca and retained
  remote structural Ca: folds 2881/2912/2913, audits 2884/2925/2931.
- Spectator preparation 2969/2985: all 15 pass, source heavy coordinates retained,
  old single-Ca regression exact. No global spectator-response model is claimed.
- Q9/ExaA scoring 2996–2999: 5 correct Ca. C5/ExaF scoring 3018–3022: 4 correct,
  1 inconclusive; Q88/PedE: 3 correct, 2 wrong. Merger 3049 combines both batches
  and the prior three repairs.

PedE wrong samples 1/2 were inconclusive at origin. Numerical, state and selection
checks pass; boundary contact does not distinguish failures. Physical cause
unresolved. Read two_ca_scoring/PED_E_RECOVERY_FINDINGS.md. Do not retune bands,
choose favorable sources or invent a five-fold vote from these inspected results.

## Resume without rerunning

Authoritative workspace: workspaces/ca_coverage_recovery_20261012/
two_ca_scoring_final_v1/{RESULT.json,SOURCES.tsv,REPORT.md,EXECUTION_3049.json}.
ROOT_ACK.json acknowledges the actual completion. Compact results and exact pins
are in this diagnostic directory; raw workspaces and failed 2877/2878 remain.
Historical absolute-v1 parity remains failed. All 28 recovery scheduler receipts
total 74,388 allocated CPU-seconds / 387 allocated GPU-seconds, including failures
and input support.

Recovery is finished. Root owns integration; all preparation, fold and scoring
agents completed scoped reports and vault notes. No new molecular submission
follows automatically. The next question is the PedE accommodation error, not
more coverage repair. Production and original records remain unchanged.

At final review, LanM 2665/2666 still run and collector 2669 waits; old reports
2832/2833 still appear RUNNING on node-48-256g-18. Preserve under Jacob's explicit
no-cancellation instruction. Inspect live scheduler/receipts before acting.
Completion watchers remain armed; readiness probes stay off. All new science
uses Slurm under current policies. Commit only owned paths; unrelated dirty
work remains untouched.
