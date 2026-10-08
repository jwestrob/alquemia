# Rules inherited by every Claude session and subagent

## Mandatory startup reading

Before doing task work, read the applicable CLAUDE.md and AGENTS.md files: global Claude instructions under the active CLAUDE_CONFIG_DIR (currently /home/jwestrob/jwestrob/.claude/CLAUDE.md), /home/jwestrob/.claude/CLAUDE.md when applicable, canonical repository and worktree instructions, and any ancestor/subdirectory instructions governing your files. Locate them rather than assuming a repository-root CLAUDE.md exists. Identical copies need not be read twice. Record the actual paths read in STATUS.json under instructions_read; report missing required instructions instead of claiming to have read them. Pass applicable instructions and this task's rules to every child, and require the child to read the files relevant to its work.

Jacob's current explicit instruction is authoritative: Claude API sessions run on the login host, but computation runs through Slurm. Historical instructions mentioning login-node analysis, bypassing a wedged queue, SSH execution or non-Slurm daemons do not authorize those actions for this task. Before submitting jobs, also read /home/jwestrob/jwestrob/obsidian-vault/biotite-cluster-ops.md; current user instructions and actual scheduler allocation take precedence over historical operational examples.

## Goal and scientific scope

Develop a defensible La/Dy discriminator across the consumed Hans8DQ2, Hans8FNR and Mex8FNS LanM sources. A completed review, converged endpoint, fast calculation or negative experiment alone does not complete the project. Preserve PQQ production and the four accepted embedded endpoints. Keep source water inventories, occupancy, spectator metals, physical/effective spin and frozen-f assumptions explicit. Use no PQQ bands, invented biological labels, favorable-source selection or raw cross-element totals as affinities. Reserved SpyCI-LAMBS outcomes remain unopened.

## Claude on the login host; scientific compute through Slurm

Run Claude sessions and their subagents headlessly on the login host. They are API clients/orchestrators and MUST NOT be submitted as Slurm jobs or kept waiting for a compute allocation. Login-host file review, editing, git work, lightweight orchestration and scheduler commands are allowed. Molecular preparation that evaluates a potential, ORCA and integral-producing utilities (including orca_2json), MACE/OpenMM/GFN2, tests/validation that perform calculations, scientific data processing, and benchmarks MUST run through Slurm. Treat model-loading probes and supposedly quick calculations the same way; do not try them on the login host first. Never run those directly on the login host or by SSH to an unallocated compute node.

Each top-level Claude session is explicitly authorized to prepare, submit, monitor, collect and recover its OWN finite Slurm jobs within its task specification. It does not wait for Codex or Jacob to approve each submission. Record the manifest, question, resource layout and job IDs first; obey scientific prerequisite gates. Subagents can also own disjoint jobs when their parent explicitly assigns that ownership; the parent reconciles the complete job inventory and prevents duplicate chemistry. Different sessions must not operate the same manifest, workspace or job.

No account switching, partition evasion, priority changes or modification of someone else's jobs. Test partition is only for tiny explicitly authorized infrastructure probes, never science. Exclusive molecular jobs use --mem=0 and allocation-derived MPI/MaxCore with headroom; shared jobs respect their actual allocation. Scheduler failure blocks only dependent compute: continue useful review/code work. No unrestricted retry cascades.

## Delegation and ownership

You are explicitly authorized to invoke your configured reviewer and analyst subagents. They may delegate within the configured depth/concurrency limits. Give every child these rules, its exact task, inputs, ownership and acceptance criteria. Use at most two concurrent child agents per session and two nesting levels; await foreground child results before finalizing. Do not start independent Claude CLI sessions to evade those limits. Do not use an unrestricted shell as a substitute for a denied tool.

Each top-level task gets a dedicated branch/worktree prepared by the integrator. Read canonical archived artifacts by explicit absolute path; do not copy scientific results into a new identity or mistake a worktree's missing workspaces for absent data. Edit only the assigned scope. Never modify production, shared status, environments, locks or archived results. Commit only your own task files; no push, merge, stash, reset, clean or broad git add. Return commit IDs to the integrator, who owns shared runners/scoring and commits/pushes integrated work. One task failure may redirect the project; it must not trigger an unbounded recovery cascade.

## Reporting and wake

Keep a short CHECKPOINT.md plus machine-readable STATUS.json in the assigned output directory. Status is one of prepared, running, waiting_slurm, blocked_auth, blocked_scheduler, needs_decision, failed, complete. Retain the task ID, session UUID, specification hash, worktree/commit, owned files, jobs, actual artifacts, checks, missing evidence, measured costs and one precise next action. Missing science remains null. Distinguish parser/unit tests from actually executed molecular qualification.

Do not poll Slurm through repeated model turns. Submit/record the finite job, install the existing completion watcher on worker plus collector, checkpoint and end the session phase with waiting_slurm. Root receives the event and resumes the responsible session. A Claude exit code alone does not prove its task or molecular job succeeded. Nested tasks must return their reports before the parent declares completion. Keep full transcripts on disk and the final summary under 200 words.
