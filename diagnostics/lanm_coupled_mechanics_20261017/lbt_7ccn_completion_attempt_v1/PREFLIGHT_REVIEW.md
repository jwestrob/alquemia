# Independent actual preflight review — 15457

**Verdict: pass; no blocker to root binding the declared run.** Actual preflight reports saved-input admission with no model initialization, new request or molecular call. This is input/execution readiness, not a convergence or accuracy result.

## Exact pins

| Artifact | SHA-256 |
| --- | --- |
| Actual preflight `RESULT.json` | `2fde1d8dbcc0ef7d5c1e7d9ca268c824253d107f37d060e7b673d625869b8e52` |
| `CONTINUE.py` | `551d964badc41141cfc53a85d704135ec74c21f185e12b6cc50d578dc47910d2` |
| Preflight configuration | `664ffb95e45d297904364e888d12dac63f0895b848f79e9d64de39afef52c710` |
| Scoped optimizer | `c5092f2480fbf393d6ca0d67733ebf2131580c922a62a67925dea8945eb5acb9` |
| Exact restart checkpoint | `317c2b24eff75a2af4a6fb7066f1754d1df1b9a284405a8aad2ea3b20d5c74cd` |

Actual result: `workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_completion_attempt_v1/preflight/execution/RESULT.json`.

## Verified saved evidence

- Protocol/stage/status are the expected completion-attempt preflight and `preflight_passed`; the sole row is 7CCN–La, `saved_inputs_admitted`. Scientific binding exactly matches the frozen configuration, whose future preflight/run-decision fields remain null.
- The admitted checkpoint retains92 accepted steps and92 distinct evaluations, eight history pairs and trust0.1Å. Original initial POINT `50a597c8`, terminal `eb96df29`, prepared source `4d064ff4`, numerical settings and original-source guard are preserved. The reported guard is an exact saved match to scout15254, passes physical/state admission and covers all256 movable atoms.
- The archived union is exactly113 records/109 unique keys. Every cache row and selected archived authority matches that union. The checkpoint's92 historical requests exactly equal original80 plus scout12; new requests are empty. Fresh loader aliases have not replaced the archived authorities.
- Declared prior limits92/120, cumulative172/292 and additional80/200 match configuration and evidence. Whole charge−2, La multiplicity1, target255 and the original model/solvent recipe match the exact terminal.
- Actual execution used the qualified `node-224-2t-8gpu-1`, one CPU and no visible GPU. Every native-call, new-evaluation, model-load/native-time and cache-hit counter is zero; `model_initialized=false`. Producer elapsed time is12.308673535 seconds; scheduler accounting remains authoritative.

The successful actual preflight exercised the reviewed saved component/gradient/cache/guard admission on Slurm. This independent review checked reported evidence, literal JSON identities and158 explicit file/hash pairs only; it performed no scientific imports, numerical replay, submissions or molecular work.

## Remaining gate and limits

Run still requires root's exact actual-preflight/run-decision binding and declared12CPU/oneH200 allocation. No source change or extra test is requested. The attempt retains original force criteria and finite stopping rules. `physical_state_qualified`, `force_accuracy_qualified` and `scientific_accuracy_claim` remain false; `no_automatic_followup=true`. Prior reference support for one incremental segment is not a guarantee of subsequent force accuracy, convergence, stable electronic state or preference.
