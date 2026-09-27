# Slurm-to-Codex completion wake repair — 27 September 2026

User explicitly pauses molecular work until completion notification works and authorizes tiny test-partition tests. No production/PQQ jobs will be changed.

Question: does the installed supported `codex queue` interface deliver an explicitly automated terminal-job event into the existing root thread, without a user nudge? A successful CLI acknowledgement alone is insufficient. Root must acknowledge the unique receipt. Active-wait delivery and delivery after a finalized root turn are separate tests.

Installed running root binary: `/home/jwestrob/.local/codex/lib/node_modules/@openai/codex/node_modules/@openai/codex-linux-x64/vendor/x86_64-unknown-linux-musl/bin/codex`, version0.153.4. PATH binary0.157.1 is newer and is not used for delivery. Root thread `01a0a63a-ee36-7483-b726-1f7f5b7f75cc`.

First interface probe queued successfully, messageID `01a0e27d-eedf-7f43-aa90-92e19eda1701`; root acknowledgement pending. No Slurm molecular work in this test.

Initial 32MiB submission was rejected before creating a job: test registers only 1MiB on its sole idle node. Actual tiny jobs request test partition, one task, one CPU, site-default whole-node memory (`--mem=0`) and a two-minute limit; only sleep and intentional exit0/exit7. The partition may allocate exclusively under site policy; scripts do no compute and no memory stress. Root receives terminal state and accounting, not email.

Official capability reference: https://developers.openai.com/codex/app-server ; installed `codex queue --help` establishes the local command. `codex-watch` only handles process/transport recovery and does not watch Slurm. No TTY injection, external resume of the running session, or invented human message is used.

## Actual interface result

Root explicitly acknowledged queue-probe-01 at2026-09-27T10:52:55UTC after its final reply. The official queue generated the next root turn without human input. An active root turn did not consume it. This resolves finalized-idle wake capability, not only mailbox delivery during a waiting turn.

## Real Slurm tests

Jobs1219311(exit0) and1219312(exit7),45-second sleeps, submitted ontest. Detached monitorPID2696269 waits bothterminal and queues `slurm-terminal-wake-proof-01`. Runtime receipts are in `workspaces/slurm_wakeup_20260927/`. Root acknowledgement remains required before declaring end-to-end success.
