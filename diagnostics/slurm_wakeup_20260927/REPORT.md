# Completion wake repaired and demonstrated

The supported installed `codex queue` command resumes this root conversation after its current turn finishes. Both a direct interface probe and a real Slurm completion event started a new root turn without a human message. Root explicitly acknowledged both events.

## Actual tests

- Interface `queue-probe-01`: root acknowledged after a final reply at2026-09-27T10:52:55UTC. Active-turn/wait-agent did not consume the message; it remained durably queued until final.
- Slurm success1219311: COMPLETED, exit0,45seconds,1CPU,testpartition.
- Slurm intentional failure1219312: FAILED, exit7,45seconds,1CPU,testpartition.
- Detached watcherPID2696269 saw both terminal at10:54:52UTC and enqueued messageID01a0e280-e380-7133-a972-fed6c7fc4256 at10:54:53UTC. Root then automatically resumed on `slurm-terminal-wake-proof-01` and acknowledged actual receipt.
- Real replay of the same watcher command detected the successful receipt and declined duplicate notification. Python syntax compilation passed.
- Total90allocatedCPU-seconds, zeroGPU, no molecular calculation, no other job modified. Initial32MiB submission failed before allocation because test registers only1MiB; successful sleeps used site-default `--mem=0`. No large-memory workload ran.

## Why earlier notifications failed

A saved collector result or email does not start a Codex turn. Subagent mailbox delivery is also distinct from starting an idle root turn. The previously used machinery did not invoke the supported queue interface. No session restart, terminal injection or daemon change is needed.

The running root binary is0.153.4, whereas shell PATH resolves0.157.1. The repair explicitly uses the running binary. Its `queue --help` advertises the supported operation. General official app-server turn/event documentation: https://developers.openai.com/codex/app-server . No unsupported assumption from that document substitutes for the local end-to-end test.

## Operating rule

Arm `watch_and_queue.py` with explicit owned job IDs, the current root UUID, pinned executable, stable receipt and unique event. Record PID and command. After a short progress reply ends the turn, terminal status automatically starts the next turn; no manual poll or user ping is necessary. For short interactive monitoring, waiting in tools remains possible but is not the mechanism tested here.

A receipt covers all supplied jobs. Include collectors when their outputs are needed before continuing. File locking and successful-enqueue deduplication prevent duplicate monitors/events. A crash during enqueue leaves a visible ambiguous state, requiring queue inspection rather than automatic replay. Host reboot kills the watcher: rearm it on original IDs from the checkpoint, never resubmit science for a missing alert. CLIclosed/authusage failures are limitations, not silent success.

See COMMANDS.md for exact tested reattachment and RESULT.json for actual accounting and root acknowledgement. Science remained paused throughout repair; integration into shared instructions and next scientific watcher belongs to root.
