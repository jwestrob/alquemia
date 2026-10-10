# Slurm readiness monitor — 9 October 2026

Running detached on the login host, PID **669874**. Jacob requested a
probe every 30 minutes and automatic resumption after readiness, without another
human check-in. First attempt was rejected for invalid account/partition
association; no job ID allocated. Next scheduled attempt: 19:21:39 Pacific.

## Behavior

- One tiny test-partition probe at a time: hostname, date and a unique marker.
- Submission attempts at least 1800 seconds apart. Pending/running probes are
  checked every 30 seconds; no replacement is submitted while one is outstanding.
- Readiness requires top-level COMPLETED, exit 0:0 and the exact stdout marker.
  Then `codex queue` wakes this thread with [resumption instructions](ON_READY.md).
- Inherited SBATCH overrides are removed; the command pins test/one node/one task.
  No hardcoded memory, time or CPU-per-task request. No scientific compute here.
- File lock and persisted submission/notification intents prevent blind duplicate
  retries. An uncertain submission wakes root for reconciliation instead of
  resubmitting. A queue failure is preserved visibly for recovery.
- Claude sessions remain on the login host; scientific work goes through Slurm.
  Root checks current task ownership and scientific gates before dispatch.

## Verification

Five administrative tests passed (0.016 seconds), using explicitly mocked
scheduler responses: top-level record selection, readiness requirements,
submission ambiguity, queue intent persistence, and rejected → pending → success
flow with 30-minute spacing and no duplicate wake on restart. No mock is
scientific evidence. Run `python3 diagnostics/slurm_readiness_20261009/test_monitor.py`.

Actual first sbatch rejection and live process were verified. Current Codex CLI
queue interface accepted a separately labelled installation self-test; actual
root delivery is confirmed only by its separate acknowledgement receipt.

## Receipts and restart

- Live state/attempts: `workspaces/slurm_readiness_20261009/state.json`.
- Exact argv, PID, code hash and rearm command: same directory, `LAUNCH.json`.
- Wake installation test: `wake_interface_test.json`; root delivery, when seen:
  `wake_interface_ROOT_ACK.json`.
- Monitor stdout/stderr: `monitor.log`.

The process is detached and survives this turn ending; a login-host reboot kills
it. Inspect its saved PID, state and existing Slurm job before using the exact
`rearm_command` in LAUNCH.json. A queued/uncertain receipt requires reconciliation,
not deleting state. A closed Codex client, usage/authentication limits or a host
outage may delay processing a queued wake. Test readiness alone does not establish
that every scientific partition accepts jobs; verify actual production resources
on resumption. No molecular jobs, scientific results or production changes here.
