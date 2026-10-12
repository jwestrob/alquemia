# Source-boundary drain supervisor

Root launched `DRAIN.py` as PID1952935 after inspecting it. The launched file is
SHA256 `175670b50febf41bf5eb5bfcfb1efe2cb7324a4780e59582c304bbfe2b48c017`.
It remains unchanged. This agent ran only the read-only preflight and self-tests;
all real signals are owned by root's launched supervisor.

## Checks before launch

`python DRAIN.py selftest` passed six no-signal checks: safe-boundary decision,
waiting for an active call, flushing a dispatched request, preserving a raced
optimizer, rejecting a different process start identity, and emergency thaw
ordering. No model, optimizer, signal or scientific result was generated.

The actual `preflight` passed all five PID/start-tick/command/UID/parent-child
identities. Its preserved receipt is `DRAIN_PREFLIGHT.json`. At that observation,
49 sources were terminal and two were active. The process identities were:
coordinator1620438, parents1622701/1622702, models1622703/1622704.

## Operation and receipt

The supervisor pauses the coordinator, waits for each active source's actual
terminal receipt, then pauses its parent. A raced optimizer is allowed to finish
its full source. Any dispatched next-origin/cross call must reach a real terminal
receipt. Both lanes must be safe before either is terminated. Process CPU ticks,
clock frequency and peak RSS are captured before model termination. Finally the
coordinator resumes and collects the partial batch. No frozen science files,
parameters, optimizer settings or source results are changed.

Root's command:

```bash
python diagnostics/orca_parameter_audit_20261012/slurm_resumption/DRAIN.py drain --batch workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1 --output workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/DRAIN.json
```

Do not run a duplicate. `DRAIN.json` must report `status: drained`, the original
manifest pin, and empty `live_owned_processes`, `stopped_owned_processes` and
`incomplete_optimizer_task_ids` before the Slurm wrapper accepts it. The original
batch coordinator's nonzero exit is an expected operational stop. Actual complete
receipts remain reusable; the drain does not invent terminal scientific results.

## Known cleanup edge in the launched version

Independent review confirmed the source/call boundary logic. It also found that
`send()` records a stopped PID after SIGSTOP and the event-log fsync. An exception
inside that narrow interval can omit a stopped target from normal thaw cleanup.
Root was told immediately after launch; this file is preserved rather than
altered while running. If any supervisor exception occurs, root must inspect all
five exact recorded process identities for stopped state, including identities
not listed in the internal stopped map. Future versions should record STOP
ownership before sending the signal. This report does not claim that edge is fixed.
