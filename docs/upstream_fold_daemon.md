# Upstream fold_daemon — architecture, resurrection, gate

**Audience:** any agent (any model, any session) that needs to understand, restart, or modify the pipeline feeding this repo's `inbox/`. Written to be read cold — no assumed context from a prior conversation.

**Prerequisite reading:** `SESSIONS.md` for who owns what right now.

**Live code location:** `/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/fold_daemon/` (outside this repo — the fold_daemon is upstream infrastructure feeding the DFT verifier). Its own README is at `<that dir>/README.md` and covers the drop protocol for one-off custom folds.

## What feeds `alchemical_bvs/inbox/`

Three long-lived processes:

```
     manifest (results/fold_manifest.tsv, 8929 candidates)
             │
             ▼
  ┌─ LOGIN NODE (biotite) ────────────────────────────────┐
  │  inbox_feeder.py                                       │
  │  drops <acc>.faa symlinks → fold_daemon/inbox/         │
  └────────────────────────────────────────────────────────┘
             │
             ▼
  ┌─ node-224-2t-8gpu-1 (under --gres=gpu:0, or plain ssh) ┐
  │  fold_daemon.py                                        │
  │  1. sbatch msa_prep.sh on `standard` partition        │
  │  2. picks nvidia-smi-idle GPUs opportunistically       │
  │  3. runs `protenix pred` with CUDA_VISIBLE_DEVICES     │
  │  4. writes rank_0..2 CIFs to processed/<acc>/          │
  └────────────────────────────────────────────────────────┘
             │
             ▼
  ┌─ LOGIN NODE (biotite) ────────────────────────────────┐
  │  forward_to_dft.py                                     │
  │  reads processed/<acc>/<acc>_rank_0.cif                │
  │  runs GEOMETRY GATE (see below)                        │
  │  PASS → symlink to alchemical_bvs/inbox/               │
  │  FAIL → record in _geometry_gate.tsv, do not forward   │
  └────────────────────────────────────────────────────────┘
```

Everything downstream of `alchemical_bvs/inbox/` is this repo's own DFT discriminator pipeline (owned by a different agent — see `SESSIONS.md`).

## The geometry gate

Filter running inside `forward_to_dft.py`. Full rationale is in the docstring of that file at the constants block; short version:

- **Cutoff:** 3.0 Å (La-O bond ~2.6 Å; 3.0 is generous, 3.5 was too permissive)
- **Count:** every O and N atom (protein-side and HETATM/ligand-side) within cutoff of the La in the predicted CIF
- **Pass criteria:** `n_total ≥ 4` AND `mean_d ≤ 2.7 Å`
- **Calibrated 2026-05-18** against 11 La-verified + 14 Ca-verified PQQ controls AND 278 DFT-rated protenix folds

Empirical properties on the calibration data:
- All 25 PQQ controls pass (they're all real coordination pockets; discrimination is DFT's job, not the gate's)
- 75% (12/16) of DFT-Ln-preferring protenix folds preserved
- 95% (95/100) of DFT-EXCLUDE-under-carved cases rejected before wasting DFT compute
- Rejects 34% of DFT-Ca-evolved (fine — those go to DFT to confirm anyway)

**If you want to change the threshold:** edit the constants `GATE_CUTOFF`, `GATE_MIN_DONORS`, `GATE_MAX_MEAN_D` at the top of `forward_to_dft.py`. Then re-run gate over `processed/` retrospectively (see "Retroactive rescoring" below) and restart the forwarder.

State files:
- `fold_daemon/_geometry_gate.tsv` — every gate decision (columns: accession, n_prot_O, n_prot_N, n_het_O, n_het_N, n_total, mean_d, decision, ts)
- `fold_daemon/_forwarded.txt` — accessions that PASSed and got symlinked (legacy state, kept for backward compat)

Both files are append-only. Idempotent — a fold already in either is skipped on next scan.

## Where state lives (survives biotite restart)

Everything the pipeline needs to resume is on group NFS (`/groups/banfield/...`), not on user home dirs. Concretely:

| File / dir | Purpose |
|---|---|
| `fold_daemon/inbox/` | Pending FASTA drops (symlinks from feeder) |
| `fold_daemon/claimed/<acc>/` | In-flight workdirs (recovered on daemon restart) |
| `fold_daemon/processed/<acc>/` | Completed folds (rank_0..2 CIFs + confidence JSONs + _spec.json) |
| `fold_daemon/errored/<acc>/` | Failed folds with debug info |
| `fold_daemon/_pids/*.pid` | `<pid> <hostname>` per running daemon (used by `stop_all.sh`) |
| `fold_daemon/_geometry_gate.tsv` | Gate decisions per accession |
| `fold_daemon/_forwarded.txt` | Successfully-forwarded accessions |
| `fold_daemon/logs/` | Daily rotated logs per process |

If biotite restarts, the running processes die. The files above persist. Follow the resurrection procedure below.

## Resurrection after biotite restart (or any cold start)

Follow these steps in order. Total time: ~2 minutes.

### 0. Assess state

From the login node (biotite):

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/fold_daemon
bash stop_all.sh --status
```

Expected after a biotite restart: all three daemons show `NO PIDFILE` (files stale from before restart), OR `dead pid=X` (pidfile there but process gone).

### 1. Clean stale pidfiles

```bash
rm -f _pids/fold_daemon.pid _pids/inbox_feeder.pid _pids/forward_to_dft.pid
```

The daemons refuse to start with a live-looking pidfile from the same host. Explicit `rm` avoids the "was it just paused?" ambiguity.

### 2. Restart `fold_daemon` on the GPU node

This is the tricky one. The daemon needs to run on `node-224-2t-8gpu-1` with an interactive-shell PATH (pyenv). SSH non-interactive doesn't pick up `.bashrc`, so PATH must be set explicitly in the SSH command.

**Correct one-liner (verified working):**

```bash
DAEMON_DIR=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/fold_daemon
ssh -o BatchMode=yes node-224-2t-8gpu-1 \
  "export PYENV_ROOT=\$HOME/.pyenv; \
   export PATH=\$PYENV_ROOT/bin:\$PYENV_ROOT/shims:\$PATH; \
   cd $DAEMON_DIR && \
   nohup setsid python3 fold_daemon.py > logs/fold_daemon_stdout.log 2>&1 < /dev/null &"
```

Wait ~5 seconds and verify:

```bash
ssh -o BatchMode=yes node-224-2t-8gpu-1 'pgrep -af "python3 fold_daemon.py" | grep -v bash'
cat $DAEMON_DIR/_pids/fold_daemon.pid
tail -5 $DAEMON_DIR/logs/fold_daemon_$(date +%Y-%m-%d).log
```

Expected log content: `[INFO] recovered N in-flight workitems` (from `claimed/` dir), then GPU-idle-detection cycles.

**Common failures (in order of frequency):**

1. **`FileNotFoundError: 'protenix'`** — PATH didn't inherit. You forgot the `export PYENV_ROOT=...` in the SSH command. Kill the bad daemon (`ssh node-... 'kill <pid>'`), clean pidfile, retry with the full command above.
2. **Two daemons running** — you launched twice without checking. Kill both, clean pidfile, launch once. Confirm with the pgrep call above.
3. **SSH prompts for password** — should not happen; the user has passwordless key. If it does, the identity is `~/.ssh/biotite_internal`, config at `~/.ssh/config` under `Host node-*`.

Alternative: **`sbatch sbatch_daemon.sh`** launches the daemon under a SLURM `--gres=gpu:0` allocation. Cleaner accounting, but subject to `QOSMaxJobsPerUserLimit=10` on the `standard` queue. Use SSH launch if you're already at the QOS cap.

### 3. Restart `inbox_feeder` on login node

```bash
cd $DAEMON_DIR
setsid nohup python3 inbox_feeder.py > logs/inbox_feeder_stdout.log 2>&1 < /dev/null &
```

`setsid nohup` + input redirection is required so the process survives shell disconnect. Historical bug: on 2026-05-10 the feeder was launched without `setsid`, the user's zellij pane closed 6 minutes later, and the pipeline went idle for 14 hours before we noticed.

### 4. Restart `forward_to_dft` on login node

```bash
setsid nohup python3 forward_to_dft.py > logs/forward_to_dft_stdout.log 2>&1 < /dev/null &
```

### 5. Verify

```bash
bash stop_all.sh --status
```

Expected: all three daemons `ALIVE` with fresh PIDs; hostnames `node-224-2t-8gpu-1`, `biotite`, `biotite` respectively.

Tail each log for ~30 seconds to confirm activity:
```bash
tail -F logs/fold_daemon_$(date +%Y-%m-%d).log
```

fold_daemon should show either "started fold X on GPU Y" (if inbox has work) or "GPU X idle, but nothing fold-ready" (if inbox drained).

## Pausing the pipeline (safe operation)

To free GPUs for other work without losing state:

```bash
bash stop_all.sh fold_daemon
bash stop_all.sh inbox_feeder
# leave forward_to_dft alive; it just watches processed/ and is idle when there's nothing new
```

To resume: repeat steps 2 and 3 of the resurrection procedure.

## Retroactive rescoring (when you change the gate)

If you edit the gate constants and want to re-decide the historical `processed/` folds:

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/fold_daemon
# 1) back up the current TSV
cp _geometry_gate.tsv _geometry_gate.tsv.bak_$(date +%Y%m%d_%H%M%S)
# 2) delete + regenerate — the gate code will re-evaluate every processed/<acc>/*_rank_0.cif
rm _geometry_gate.tsv
# 3) run forward_to_dft in --once mode (add flag if not present, else briefly start+stop)
python3 forward_to_dft.py --once
```

The `--once` flag runs a single scan cycle then exits — perfect for retroactive re-scoring. All `processed/` accessions not in `_forwarded.txt` get gate-evaluated and written to `_geometry_gate.tsv`.

Note: accessions already in `_forwarded.txt` are skipped (they were already sent to DFT under the old gate; re-sending would be a duplicate to the DFT verifier). If you want to re-evaluate everything, back up and delete `_forwarded.txt` too, but this will cause the forwarder to re-symlink everything into `alchemical_bvs/inbox/` — check first that the DFT verifier is idempotent on re-drops.

## Common cross-boundary issues with the DFT verifier

1. **DFT verifier expects specific file naming.** `forward_to_dft.py` produces names like `<accession>_<clade>.cif` (clade looked up from `results/fold_manifest.tsv`) or bare `<accession>.cif` if no clade. The DFT verifier accepts either.

2. **PQQ co-fold candidates.** If a candidate needs re-folding *with* PQQ as a co-ligand (per `results/pqq_refold_queue.tsv`), drop a JSON manifest into `fold_daemon/inbox/`:

    ```json
    {
      "name": "<accession>_pqq_la",
      "fasta_path": "/abs/path/to/protein.faa",
      "ligands": ["LA:1", "PQQ:1"],
      "samples": 3, "seeds": "101"
    }
    ```

   Both LA and PQQ are valid CCD codes; `foldit.py` auto-prepends `CCD_`. Known limitation: protenix may not place PQQ near La — verified 2026-05-18 on A0A840IK71 (PQQ landed 12-47 Å away in all 3 ranks). Track the outcome of any PQQ refold before scaling.

3. **carve_generic.py bugs downstream.** The DFT verifier has had "empty-carve" incidents. If you see clusters of identical bit-for-bit empty ORCA outputs across unrelated proteins, that's not our (upstream) problem — flag to the DFT-discriminator agent via `SESSIONS.md`.

## Vault cross-references (rationale, not authoritative)

These live on the user's personal home dir and are historical documentation. Repo docs above are authoritative.

- `~/jwestrob/obsidian-vault/agent-captures/2026-05-11_fold-daemon-opportunistic-gpu-pipeline.md` — original design capture, ~250 lines
- `~/jwestrob/obsidian-vault/biotite-cluster-ops.md` §10 — cluster-ops perspective
- `~/jwestrob/obsidian-vault/projects/lanthanide-binding/CONTEXT.md` — project-level context including the active pipeline pointer
