# `scripts/` — per-script reference

Compact reference for every script in this directory. For deeper write-ups
(failure modes, idempotency, data flow), see
[PIPELINE.md](../PIPELINE.md) at the repo root.

Conventions:

- All scripts assume the workspace root is the parent directory of `scripts/`
  (i.e. `alchemical_bvs/`).
- Production interpreters: `lanm_qmmm` for almost everything,
  `fep` (PDBFixer + OpenMM) only for `protonate_cif.py`. See
  [environment.yml](../environment.yml).
- Listing order follows logical pipeline order rather than alphabetical.

---

## Pipeline core (production)

### `process_inbox.sh`

Top-level driver. Walks `inbox/*.cif`, detects PQQ via composition, routes
through the appropriate sub-pipeline, writes a SLURM submit script, and
`sbatch`s it. Idempotent; flock-guarded.

- **Invoke:** `bash scripts/process_inbox.sh` (env: `MAX_SUBMIT=N` to cap)
- **Inputs:** any CIFs in `inbox/`
- **Outputs:** `<stem>_qm/` per CIF; consumed CIFs moved to `inbox/processed/`
- **State:** `/tmp/process_inbox.lock` (flock mutex)

### `normalize_af3_cif.py`

Rename Protenix/AF3 `LIG_*` residues. La-bearing → `LA`; PQQ-like
composition (≥14C + ≥6O + ≥1N) → `PQQ`. Output: PDB.

- **Invoke:** `python scripts/normalize_af3_cif.py <src.cif> <out.pdb>`
- **Only used:** on the PQQ branch of `process_inbox.sh`

### `protonate_cif.py`

PDBFixer wrapper. Adds hydrogens at pH 7. Reads CIF or PDB; writes PDB.

- **Invoke:** `python scripts/protonate_cif.py <src> <out.pdb>`
- **Requires:** the `fep` env (PDBFixer + OpenMM)
- **Caveat:** non-deterministic protonation (15–66 kcal/mol absolute-energy
  drift between runs), but the noise cancels in ΔΔE.

### `carve_generic.py`

Generic QM cluster carver. Picks the metal site (most carboxylate donors, or
honour `--site-chain` / `--site-resnum`), enumerates first-shell donors,
carves ASP/GLU/ASN/GLN/SER/TYR sidechains + backbone-O contacts + inner-shell
waters, writes three XYZ files + a SLURM submit script.

- **Invoke:** `python scripts/carve_generic.py <pdb> <out_dir> --stem <stem>`
- **Optional:** `--site-chain`, `--site-resnum`, `--first-shell-cut <Å>`
- **Outputs:** `{stem}_{La,Ca,apo}_qm.xyz`, `sp_{stem}_{La,Ca,apo}.inp`,
  `submit_{stem}.sh`
- **Key constants:** `FIRST_SHELL_CUT = 3.0` Å; `SIDECHAIN_QM_ATOMS` dict
  (see [METHODS.md §3](../METHODS.md))

### `carve_with_pqq.py`

PQQ-aware carver. Same logic as `carve_generic` plus the 24 PQQ heavy atoms
at formal charge −2.

- **Invoke:** `python scripts/carve_with_pqq.py <pdb> <out_dir> --stem <stem> --metal-chain <chain> --metal-resname <name>`
- **Optional:** `--replace-metal La` (used when an X-ray structure has Ce/Y in
  place of La)
- **Key constants:** `FIRST_SHELL_CUT = 3.2` Å (slightly looser);
  `PQQ_CHARGE = -2`

### `update_results_jsonl.py`

Aggregator + classifier. Walks `*_qm/` dirs, parses the four `sp_*.out`
files, computes ΔΔE, classifies, appends one JSON object per candidate to
`results/all_results.jsonl`, regenerates
`results/discriminator_panel_LATEST.tsv`, and auto-shells out to
`rebuild_ln_tracker.py`.

- **Invoke:** `python scripts/update_results_jsonl.py .` (workspace root)
- **Optional:** `--rebuild` (wipe and re-process)
- **Idempotent:** mtime-based skip; append-only; deduped on read

### `rebuild_ln_tracker.py`

Rebuild `results/ln_class_hits.tsv` from `all_results.jsonl`, enriched with
spicy_lams Sharur DB annotations (Pfam, KO, organism bin, novelty) and
`known_false_positives.tsv`.

- **Invoke:** `python scripts/rebuild_ln_tracker.py`
- **Auto-fires:** at the end of every `update_results_jsonl.py` run
- **Outputs:** `results/ln_class_hits.tsv` (17 columns including rank)

### `process_recarve_queue.py`

Re-carve and re-submit driver for historical empty-carve entries. Reads
`results/recarve_queue.tsv`, retries each entry at 3.0 Å (then 3.5 Å) with
PQQ-aware routing, and classifies still-failing entries as `SOLVENT_EXCLUDE`
or `CARVE_AMBIGUOUS`.

- **Invoke:** `python scripts/process_recarve_queue.py [--dry-run]`
- **Outputs:** `results/recarve_log_<UTC>.tsv` per-stem outcomes
- **Key thresholds:** `SOLVENT_EXCLUDE_THRESHOLD = 4.0` Å,
  `PQQ_DETECT_RADIUS = 4.0` Å, `RETRY_WIDE_CUTOFF = 3.5` Å

---

## Drivers and watchers

### `inbox_watcher.sh`

Long-running poll loop. Every 45 min, if `inbox/*.cif` is non-empty and
`squeue -u jwestrob | grep -v fold/taxon_scan` has free slots, runs
`process_inbox.sh`. PID-lockfile convention at `/tmp/inbox_watcher.lock`.

- **Start:** `bash scripts/inbox_watcher.sh > watcher.log 2>&1 &`
- **Stop:** `rm /tmp/inbox_watcher.lock`

### `cleanup_orca_scratch.sh`

Walks `*_qm/` directories and deletes ORCA scratch (`.gbw`, `.tmp`, etc.)
in completed candidates. Mostly obsolete now that new carves auto-clean on
success; useful for legacy dirs that pre-date the auto-cleanup tail.

- **Invoke:** `bash scripts/cleanup_orca_scratch.sh`

### `tier1_pqq_batch.sh`

Legacy direct-driver for Colin's 25 PQQ-MDH controls, pre-dating PQQ-aware
inbox routing. Kept for reproducibility of the `colinpqq_la_*` /
`colinpqq_ca_*` panel.

- **Invoke:** `bash scripts/tier1_pqq_batch.sh <source_dir> <out_dir>`

---

## Pre-filters and pre-checks

### `prefilter_spicylams.py`

Donor-count screen for bulk batches. Drops CIFs with fewer than N
ASP/GLU/ASN/GLN/SER/TYR donors within 3.2 Å of La (default N = 2). Use
before staging the spicy_lams CIFs to avoid filling the inbox with carve-fails.

- **Invoke:** `python scripts/prefilter_spicylams.py <cif_dir> [min_donors]`
- **Outputs:** writes a `_filter_passed.txt` / `_filter_rejected.txt`
  manifest in the source dir

### `analyze_panel.py`

Legacy panel-analysis utility, predates `update_results_jsonl.py`. Reads
energies from `<stem>_qm/sp_*.out` files and prints a TSV.

- **Invoke:** `python scripts/analyze_panel.py <workspace>`
- **Status:** superseded by `update_results_jsonl.py`; kept for
  reproducibility of the 2026-05-07 panel.

### `carve_size_variants.py`

One-off helper for the calexcitin cluster-size sensitivity panel
(16/32/50/100 atom carves). Hard-coded for the calexcitin structure.

- **Invoke:** `python scripts/carve_size_variants.py`
- **Result:** see [VALIDATION.md §8](../VALIDATION.md)

---

## Legacy batch drivers (pre-discriminator pivot)

These were the overnight panel runners that built the 2026-05-07
validation panel before the inbox / watcher / aggregator stack existed.
They are kept for reproducibility of the original 54-candidate panel.

### `batch_overnight.sh` / `batch_overnight_v2.sh` / `batch_overnight_v3.sh`

Successive iterations of the original panel driver. v3 includes the
alphaproteo + Marco β-roll family expansion (2026-05-07 PM).

### `setup_b97_3c.sh`

Set up the B97-3c functional-robustness panel for the top 4 (XoxF, MxaF,
tannase, calexcitin). Output: `b97_3c_panel/{xoxf,mxaf,tannase,calexcitin}/`.
Used once for the cross-functional check; see [VALIDATION.md §7](../VALIDATION.md).

---

## 8DQ2 BVS pre-check (legacy)

Predates the discriminator pivot. These were the Pass 2 BVS residual /
donor analysis tools used to motivate moving from BVS-from-density to a
vertical Ca/Ln swap.

- `precheck_8dq2.py` — main BVS pre-check driver
- `submit_precheck.sh` — SLURM wrapper for the above
- `recon_8dq2.py` — density reconstruction comparison

Outputs at `results/precheck_8dq2_*.tsv` and the legacy `REPORT.md` /
`HANDOFF.md` at the repo root.

---

## Quick map: which script does what?

| stage                       | script                          | env       |
|-----------------------------|---------------------------------|-----------|
| inbox arrival + routing     | `process_inbox.sh`              | bash      |
| AF3 / Protenix normalisation | `normalize_af3_cif.py`         | lanm_qmmm |
| protonation (pH 7)          | `protonate_cif.py`              | fep       |
| QM carve (no PQQ)           | `carve_generic.py`              | lanm_qmmm |
| QM carve (with PQQ)         | `carve_with_pqq.py`             | lanm_qmmm |
| SLURM submission            | (inline in `process_inbox.sh`)  | bash      |
| aggregation + classification | `update_results_jsonl.py`      | lanm_qmmm |
| Ln-class enrichment         | `rebuild_ln_tracker.py`         | lanm_qmmm |
| empty-carve recovery        | `process_recarve_queue.py`      | lanm_qmmm |
| watchdog daemon             | `inbox_watcher.sh`              | bash      |
| disk hygiene                | `cleanup_orca_scratch.sh`       | bash      |
| pre-filter (bulk batches)   | `prefilter_spicylams.py`        | lanm_qmmm |
