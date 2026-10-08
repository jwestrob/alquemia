# LanM work delegated to headless Claude — 7 October 2026

These are task specifications and a launch contract, not a launched campaign.
Root integrates results. Claude sessions and subagents run headlessly on the login
host and independently submit their owned scientific jobs through Slurm. Claude
sessions themselves do NOT go in Slurm. Read COMMON_RULES.md first.

## Current prerequisites

- Installed `/home/jwestrob/.local/bin/claude`: version2.1.283. Local help verifies
  `-p`, `--agents`, `--tools`, `--permission-mode dontAsk`, JSON/stream output,
  explicit session IDs and resume. No upgrade or environment change performed.
- Claude authentication is now verified: loggedIn=true, claude.ai authentication
  in `/home/jwestrob/jwestrob/.claude`. No credentials copied or printed in artifacts.
- The user-authorized tiny test-partition probe was rejected: invalid account or
  account/partition combination. Account-association queries return no rows.
  Idle nodes are not proof jobs can run. No molecular job was submitted.
- Headless parent → reviewer → analyst delegation **passed on the login host**.
  Actual Agent calls at both levels are pinned in VERIFICATION.json; no tool
  denials. The selected sonnet alias resolved to claude-sonnet-5. CLI-reported
  cost estimate: $0.3180161. No scientific commands or submissions in this test.
- Claude-session completion wake and live Slurm job submission remain untested.

## Task graph

| Task | Can overlap | Result needed before next scientific execution |
|---|---|---|
| 01 restart review | Preliminary03, mapping/ledger portion04 | Reviewed continuation with critical defects resolved |
| 02 two endpoints | Archived analysis03 and preparation04 | Accepted pair or diagnosed failure |
| 03 interpretation | Preliminary work independent of02 | Root selects a finite discriminating experiment |
| 04 coupled response | Preparation may precede03 decision | Required physical checks actually pass |
| 05 discrimination | Follows selected scientific gates | Matched source-specific utility report |

The five numbered files define inputs, ownership and completion. A process exit
never launches the next scientific module automatically. Root reads the result,
checks readiness, and dispatches the next bounded task under existing authority.
Failures redirect dependent work; they do not redefine the project as completed.
Do not reimplement the already prepared interrupted continuation or rerun its
33 passing fixture checks without a changed implementation or review finding.

## Isolation and permissions

At dispatch, create a separate `git worktree` from a recorded integration commit
and an owned `claude/lanm-<task-id>` branch. Keep canonical scientific artifacts at
their existing absolute paths; worktrees must not rewrite archived manifests.
Session outputs belong in fresh `workspaces/claude_lanm_20261007/<task-id>/` paths.
The headless parent has Bash and Slurm submission access. Scientific compute
must be submitted through Slurm; no session-level Codex approval is required for
owned jobs in scope. The read-only reviewer/analyst definitions cannot submit,
but a task owner may define execution subagents with disjoint job ownership.
These are task/permission rules, not an OS sandbox guarantee.

`agents.json` explicitly exposes Agent to two custom agent types and includes
all shared rules in each prompt. The session configuration permits two child
agents at once and two delegation levels. Nested agents must complete before
parent finalization. Both inherit the session's permission constraints. No
bypass-permissions flag, recursive shell Claude launches or automatic model
fallback. Root records an explicit chosen model and effort at each launch.
The verification run used sonnet with low effort; task models remain dispatch choices.

## Minimal qualification before task01

With authentication verified:

1. Start the read-only headless/subagent smoke test on the login host; no Slurm dependency.
2. Ask the headless session’s
   reviewer child to read COMMON_RULES.md, delegate to analyst for one distinct
   read-only check, then return the two findings. No chemistry and no repository edits.
3. Inspect the stream for actual Agent calls/results at both depths, tool denials,
   model IDs, result status and session UUID. A verbal claim to have delegated is
   insufficient. Check that the child reports the Slurm-only rule correctly.
4. When scheduler account access returns, confirm a tiny terminal job event reaches root via the existing wake watcher.
   Use the current root UUID, not a historical hardcoded session. Do not test an
   expensive molecular run to debug headless tooling.
5. Launch task01 from a fresh worktree and explicit model/effort. Independent
   preliminary03 can be a second top-level session; at most two roots concurrently.

## Launch and completion contract

Use `run_session.sh` on the login host with explicit arguments (its header
documents them), launched as a detached process with a saved PID/session UUID
and fresh log directory. Never submit this launcher to Slurm. It does not wait
for scheduler access. The session's model/effort are explicit; artifacts retain
resolved model and usage. A thin supervisor must record exit and notify root on
completion/failure; the notification integration still needs an actual test.

For scientific jobs, the owning Claude session registers a watcher for both worker and collector. The dispatcher resumes the
responsible Claude session only after evidence is available. Use explicit
`--resume <UUID>` with a new phase prompt and new phase output directory; never
`--continue` across unrelated tasks. No minute-by-minute model polling loops.
The wrapper's exit receipt is process evidence; task acceptance requires the
specified artifacts and actual scientific checks. If auth, usage or scheduler
access fails, record the blocker and stop that phase without retry loops.

## Verification sources

- Local `claude --version`, `claude --help`, `claude auth status` and
  `claude agents --help`, read October7.
- [Official headless CLI documentation](https://code.claude.com/docs/en/headless)
- [Official subagent documentation](https://code.claude.com/docs/en/sub-agents),
  including Agent permissions, nesting/concurrency controls and inheritance.
- [Official permissions documentation](https://code.claude.com/docs/en/permissions)

Scientific execution is blocked by scheduler account access. Claude review and
planning can run now. The nested headless test passed; scientific execution
remains unavailable until scheduler account access is restored.
