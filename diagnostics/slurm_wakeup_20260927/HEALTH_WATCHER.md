# Live native-SCF health wake — frozen v1

The new `scripts/slurm_native_health_watch.py` watches explicit owned job/manifest pairs, reports per-endpoint completion and queues root **while the overall job is still running**. It never submits/cancels jobs or changes inputs. This is diagnosis/notification, not a new scientific acceptance policy. It complements the existing completion watcher.

## Alert rules declared before new recovery calculations

- Last8TRAH macro error norms all>1 and less than twofold improvement from first to last: persistent poor convergence to review. The quantity is ORCA's printed `||Error||_2`, not a kcal score.
- Trial step `|dE|>1e6Eh` in at least3distinct TRAH macro iterations: repeated extreme local-solver steps. These are optimizer diagnostics, never physical accepted energies.
- Explicit native SCF/error termination or nonzero execution receipt: failure notification.
- Output unchanged for900s while unfinished: inspect process health; **not proof of a stall or permission to cancel**. No runtime ceiling or CPU budget is imposed.
- Each endpoint terminal and partial completion: report allocation CPUs minus maximum remaining endpoint-rank demand. This is structural allocation accounting, not measured utilization. Completed endpoints may leave CPUs without work despite remaining SCFs.

Thresholds deliberately diagnose gross numerical behavior rather than desired metal signs. They were developed using already-consumed failed Dy logs and are not independent validation. Root decides whether to preserve, cancel or launch a separately justified recovery. Queued events are uniquely deduplicated by job/endpoint/reason; no repeated alert flood.

## Actual verification

Four tests pass. Both saved Dy outputs trigger severe-micro-step alerts; DyB also triggers persistent macro residuals. A real DyB prefix before macro20 already triggers an alert, showing it does not wait for job completion. Both completed La outputs report convergence and no pathology. Actual four-cell/224CPU manifest gives2completed endpoints and maximum112remaining endpoint ranks, leaving112CPUs without remaining endpoint work. All actual output hashes/evidence are in HEALTH_POLICY_AND_TESTS.json.

No molecular runs or new test jobs were launched. The supported `codex queue` delivery path was already proven to resume an idle root; this watcher uses exactly that interface but has not yet been armed on a new recovery job.

## Arming / reboot recovery

Supply explicit `--manifest`, `--job`, `--thread`, `--codex`, `--receipt` and optionally `--interval 30`. Use the same matched binary and current root UUID as the completion watcher. Root must record the exact command/PID in the new job's SUBMISSION record. Launch detached using `subprocess.Popen(...,start_new_session=True,stdin=DEVNULL)` with explicit logs. End the current turn after the checkpoint so queued notifications start the next turn automatically.

The watcher verifies Slurm User ownership and pins the manifest hash. A stable receipt/file lock prevents duplicate watchers. A restart after ambiguous/failed enqueue stops visibly for queue inspection rather than duplicating it. Reboot kills the process; relaunch on the same owned job and receipt, never resubmit chemistry merely because the monitor disappeared.

Offline read-only reproduction of the retained scout (does not send a message):

```bash
python scripts/slurm_native_health_watch.py \
  --manifest workspaces/metal_environment_response_20260926/lanm_ef3_hans_scout_v1/manifest.json \
  --job 1219501 --thread 01a0a63a-ee36-7483-b726-1f7f5b7f75cc \
  --codex /home/jwestrob/.local/codex/lib/node_modules/@openai/codex/node_modules/@openai/codex-linux-x64/vendor/x86_64-unknown-linux-musl/bin/codex \
  --receipt workspaces/slurm_wakeup_20260927/old_scout_health_inspect.json --inspect-only
```

Tests: `python -m unittest discover -s tests -p test_slurm_native_health_watch.py -v`.
