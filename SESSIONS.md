# SESSIONS — cross-agent coordination log

**Purpose:** Shared knowledge base for the multi-agent workflow feeding the DFT Ca²⁺/Ln³⁺ discriminator. Each concurrent agent (or fresh model — Claude 4.7, GPT-5.6-pro, whatever comes next) writes here so successors know what's running, what's owned by whom, and where to look before starting new work.

**Not a chat log.** This is durable state. Update it when you finish a substantial chunk of work, hand off to another agent, pause a subsystem, or make a decision that another agent needs to know about. Prefer terse "what/where/why" over narrative.

## How to use this file

1. On session start, read this file top-to-bottom before touching state.
2. When you take ownership of something, add or update the corresponding subsection under **Ownership map**.
3. When you finish or pause work, append a dated entry under **Session log** with:
   - your agent identity (model + role tag if known)
   - what changed on disk (paths, TSVs, DB tables, running processes)
   - what state you left things in (running / paused / broken / open question)
   - what the next agent should do first
4. If you disagree with a prior decision, don't overwrite it — add a new entry that says "revising X because Y" and update the ownership map.

## Model / host notes for future readers

- The code paths in this repo are model-agnostic Python + shell. Any model with tool-use can execute them.
- Pyenv-managed Python (`~/.pyenv/shims/python3` on both login and GPU nodes) is the runtime. SSH sessions to GPU nodes need explicit `PYENV_ROOT` + `PATH` setup (see `docs/upstream_fold_daemon.md`).
- Vault notes at `~jwestrob/obsidian-vault/` contain historical rationale but should NOT be treated as the primary source of truth — the repo is. The vault is on the user's personal home dir; the repo is on group NFS. Repo survives if home dir is unavailable.

---

## Ownership map

### DFT discriminator core — `alchemical_bvs/`

**Owner:** DFT-discriminator agent (identity TBD).
**Scope:** everything in `alchemical_bvs/` proper — `scripts/carve_generic.py`, `scripts/process_inbox.sh`, the recarve queue, `results/all_results.jsonl`, `results/confirmed_ln_binders.tsv`, `results/discriminator_panel_LATEST.tsv`, `results/ln_class_hits.tsv`, `results/pqq_refold_queue.tsv`, `results/recarve_*` logs, and the `*_qm/` per-candidate ORCA workdirs.

**Downstream inputs it consumes:**
- `alchemical_bvs/inbox/*.cif` — new candidate structures (dropped by the fold_daemon forwarder, or by hand)
- `alchemical_bvs/inbox/processed/*.cif` — audit trail of previously-processed inputs

**Outputs it publishes:**
- `results/all_results.jsonl` (line-per-candidate ΔΔE + class)
- `results/confirmed_ln_binders.tsv` (curated positives)
- `results/pqq_refold_queue.tsv` (candidates asking upstream to refold with a PQQ ligand — the fold_daemon owner reads this)
- `results/recarve_queue.tsv` (candidates needing re-processing after a carve-code fix)

**Long-lived processes it manages:** `recarve_watcher.log` shows one. Others as they come.

### Fold pipeline — `on_density_scanner/fold_daemon/`

**Owner:** fold-daemon agent (currently: Claude 4.7 sessions initiated 2026-05-09 through 2026-07-16; hands off cleanly to whichever model reads this next).

**Scope:** everything under `on_density_scanner/fold_daemon/`. Feeds this repo's `inbox/` via `forward_to_dft.py`.

Detailed handoff and resurrection procedure: **`docs/upstream_fold_daemon.md`**.

### LanM alchemical-BVS validation — see `HANDOFF.md`

Separate task-scope, presumably a different agent's remit. Not currently in flight (no active processes matching). If you're that agent, add an ownership entry here when you begin.

---

## Session log

### 2026-07-16 — fold-daemon agent (Claude 4.7)

**Model + session:** Claude 4.7 (1M context), same session continuous since 2026-05-09.
**Reason for log entry:** User asked for durability consolidation so successor sessions/models can take over.

**What I own on disk:**
- `on_density_scanner/fold_daemon/fold_daemon.py` — v4 parallel-Popen fold worker, opportunistic GPU use on `node-224-2t-8gpu-1` via `nvidia-smi` (bypasses SLURM GRES accounting on purpose)
- `on_density_scanner/fold_daemon/inbox_feeder.py` — login-node feeder from `results/fold_manifest.tsv`, excludes pH-1-3 clades (acidithiobacillus, acidimicrobiia)
- `on_density_scanner/fold_daemon/forward_to_dft.py` — bridge from `processed/` to this repo's `alchemical_bvs/inbox/`; runs the **geometry gate** (see below)
- `on_density_scanner/fold_daemon/stop_all.sh` — SSH-aware pidfile-based kill helper (reads `_pids/<name>.pid` files, SSHes to right host)
- `on_density_scanner/fold_daemon/sbatch_daemon.sh` — SLURM-launched alternative to SSH-launch
- `on_density_scanner/fold_daemon/msa_prep.sh` — sbatch template that invokes `foldit.py --no-submit`

**Current state (as of this log entry):**
- `fold_daemon` — **paused** (user request during giant-protein folding on the same GPU node)
- `inbox_feeder` — **paused**
- `forward_to_dft` — **ALIVE** on biotite login node, PID 23156, uptime ~2 months, will still drain any straggler in `processed/` if the fold daemon resumes

**Cumulative work since 2026-05-09:**
- 1029 protenix folds completed → `processed/`
- 200 folds errored (mostly "MSA COMPLETED but no manifest" from foldit.py's post-check; known false-positive failure mode)
- 586 folds forwarded to this repo's `alchemical_bvs/inbox/`
- 1029 entries in `_geometry_gate.tsv` (all decided by the current gate)

**Geometry gate (in `forward_to_dft.py`):**
- Counts all coordinating donors (O and N, protein and HETATM) within **3.0 Å** of La
- Passes if `n_total ≥ 4` AND `mean_d ≤ 2.7 Å`
- Calibrated 2026-05-18 against 11 La-verified + 14 Ca-verified PQQ controls AND 278 DFT-rated protenix folds
- Retroactive result on 616 pre-gate folds: 202 PASS, 414 FAIL (33% pass rate)
- Rationale for these numbers is inline in `forward_to_dft.py` (constants section) and expanded in `docs/upstream_fold_daemon.md`

**Known-open decisions I'm passing to future readers:**

1. **PQQ refold (A0A840IK71_pqq_la) result is inconclusive.** Protenix put PQQ 12–47 Å from La in all 3 ranks — the "co-fold" became an effectively-apo fold. The DFT discriminator scored it against a *different* pocket than the queue's hypothesized D547/D549/D554 site. If someone wants to actually test the PQQ hypothesis on this protein: try AF3 server, or superpose XoxF/MxaF template + manual PQQ transplant. Ticket lives in `results/pqq_refold_queue.tsv`.

2. **Manifest tracking is stale.** `results/fold_manifest.tsv` shows ~217/8929 done but `processed/` has 1029 dirs. `scripts/fold_status.py` isn't running in real-time. Cosmetic; fold_daemon doesn't consult the manifest for "what to fold next" — the feeder does that on its own.

3. **Disk usage** is ~1.5 GB in `processed/` + `errored/`. Will scale linearly with fold count; audit if it approaches 100 GB.

**Next agent picking this up should:**
- Read `docs/upstream_fold_daemon.md` — full resurrection procedure, including the pyenv-in-SSH gotcha that will otherwise cost you 30 minutes
- Not rebuild the gate. If you have reason to change the threshold, do it in one place: the constants block near the top of `forward_to_dft.py`. All rationale is in the docstring there.
- If biotite has restarted since this entry: `forward_to_dft` is dead, `_pids/forward_to_dft.pid` is stale, the recovery procedure is documented in `docs/upstream_fold_daemon.md`.

---

*Template for future entries:*
```
### YYYY-MM-DD — <agent identity>
**Model + session:** ...
**Reason for log entry:** ...
**What changed:** ...
**Current state:** ...
**Open questions / next actions:** ...
```
