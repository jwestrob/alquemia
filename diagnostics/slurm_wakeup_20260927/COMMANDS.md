# Completion notifications that resume the root conversation

The matching running Codex CLI supports `queue`. Its pending message is consumed
when the current root turn finishes, starting a new turn automatically. A root
that stays forever inside a wait tool will delay this queued event; use a short
final progress reply after arming this watcher. This differs from agent mailbox
messages and emails, neither of which demonstrated idle-thread wake here.

## Reattach the tested monitor without resubmitting jobs

Run on the login host as the current user. This command checks the stable receipt
and exits without duplicate enqueue if it already completed:

```bash
python diagnostics/slurm_wakeup_20260927/watch_and_queue.py \
  --jobs 1219311 1219312 \
  --thread 01a0a63a-ee36-7483-b726-1f7f5b7f75cc \
  --codex /home/jwestrob/.local/codex/lib/node_modules/@openai/codex/node_modules/@openai/codex-linux-x64/vendor/x86_64-unknown-linux-musl/bin/codex \
  --receipt workspaces/slurm_wakeup_20260927/terminal_event.json \
  --event slurm-terminal-wake-proof-01 \
  --interval 2
```

For production monitors, use their explicit owned job IDs, a unique stable receipt
and event name, and a15-second interval. Monitor the collector's job as well as
workers when the actionable next step requires collected results. Launch with
`subprocess.Popen(..., start_new_session=True, stdin=DEVNULL)` and explicit log
files; archive itsPID/command. This survives a chat turn ending, not a host reboot.
After reboot inspect the receipt and relaunch this monitor only. Never resubmit
scientific jobs merely because a notification process disappeared.

The watcher takes a file lock; a second live copy fails instead of doubling
notifications. A recorded successful enqueue is never repeated. If the host dies
between enqueue and acknowledgement recording, the receipt remains
`enqueue_pending`; recovery stops visibly to avoid an ambiguous duplicate.
Inspect the supported session queue before deciding whether to re-enqueue.

There is no claim of wake with the CLI closed, host down, expired authentication,
or exhausted usage. The event remains a durable queued item, but processing
requires the user's existing Codex session to operate. Root receipt and timestamp
are the evidence for actual wake; CLI stdout proves only enqueue.
