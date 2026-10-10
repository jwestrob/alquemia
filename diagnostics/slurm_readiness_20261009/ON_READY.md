# Automatic LanM resumption after Slurm recovery

Jacob authorized this monitor on 9 October 2026: try every 30 minutes and wake
Codex as soon as a real test job succeeds, then resume submitting the work.
This is execution authorization, not a request for another permission check.

1. Read workspaces/slurm_readiness_20261009/state.json and the referenced
   submission/terminal/stdout records. A readiness event requires the top-level
   test job COMPLETED, exit 0:0, and its exact unique output marker. Write a
   separate ROOT_ACK.json recording actual receipt. A wake self-test or uncertain
   submission alert does not establish scheduler readiness.
2. Read current AGENTS.md, applicable CLAUDE.md, docs/AGENT_PIPELINE.md,
   diagnostics/metal_environment_response_20260926/CURRENT.md, and
   docs/lanm_claude_tasks_20261007/{README,COMMON_RULES}.md. Inspect actual jobs,
   task session logs and worktree ownership before dispatch; don't duplicate jobs.
   A test-partition success establishes that route only: inspect science partition
   access and current resources before the molecular submission.
3. Start task01 restart review in an isolated worktree using run_session.sh,
   with a recorded model, effort, session UUID, PID and fresh output directory.
   Preliminary task03 interpretation may run alongside it. Both Claude parents
   run on the LOGIN HOST, may delegate their own subagents, and may submit owned
   scientific jobs. ALL scientific computation and scientific tests go through
   Slurm; do not put the Claude sessions themselves in Slurm. Propagate current
   instructions to children. No recursive shell-Claude spawning.
4. Close any critical restart-review finding, then execute task02's already
   prepared two missing Hans8FNR endpoints, not a new campaign. Read
   diagnostics/metal_environment_response_20260926/recovery_20261007/REPORT.md
   and the pinned continuation manifest. Preserve all accepted other endpoints
   and interrupted outputs. Reconcile current allocation with native MPI and
   memory; don't reuse the obsolete September startup cutoff. Use current
   scheduler rules and normal priority. PQQ production belongs to another session.
5. Record owned submissions and arm the proven completion wake for worker and
   collector with the actual current Codex UUID/executable. Add a thin tested
   exit-and-notify supervisor for Claude tasks before leaving them detached; its
   integration remains to be qualified. Model sessions should return after
   submitting and resume on receipts, not consume tokens polling.
6. Integrate actual results and choose dependent task03/04/05 work by their
   scientific gates. No automatic promotion, hidden failures, reserved-label
   inspection, old explicit-f/vacuum restart, or blind whole-library launch.

The recovery goal is still a defensible La/Dy discriminator. Four accepted
repaired embedded origins and a small Glu91 response do not establish that goal.
Work forward from actual endpoints rather than repeating completed chemistry.
