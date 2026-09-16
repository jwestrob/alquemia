# Pipeline — runtime architecture

> **Current agent guide — 2026-09-16:** [Alquemia operations and protocol status](docs/AGENT_PIPELINE.md).
> The automatic inbox retains the existing v2 routes. Fixed-core PQQ and peptide-amide v3 are explicit separate paths. Manifested ORCA jobs use the current MPI runner with one OpenMP thread per rank.
> This guide supersedes conflicting operational/status prose below; dated scientific records remain historical.

This document is the engineering reference for the discriminator pipeline:
what every script does, in what order, with what state, and how to debug it.
For the chemistry see [METHODS.md](METHODS.md); for class labels see
[TIERS.md](TIERS.md); for the user-level how-to see [HOWTO.md](HOWTO.md).

---

## 1. Data flow

```
inbox/<stem>.cif
        │
        ▼
[ process_inbox.sh ] ──┬── PQQ-aware routing ──┐
        │              │                       │
        │  no PQQ      │  has PQQ              │
        │  (LIG_* with │  (≥12C + ≥6O + ≥1N    │
        │   <12C or    │   in a LIG_* res)     │
        │   <6O or 0N) │                       │
        ▼              │                       ▼
[ protonate_cif.py ]   │           [ normalize_af3_cif.py ]
   PDBFixer pH 7       │              LIG_*/LA1 → LA/PQQ
        │              │                       │
        │              │                       ▼
        │              │             [ protonate_cif.py ]
        │              │                PDBFixer pH 7
        │              │                       │
        ▼              │                       ▼
[ carve_generic.py ]   │           [ carve_with_pqq.py ]
   ASP/GLU/ASN/GLN/    │              same dict + 24-atom PQQ
   SER/TYR sidechains  │              at formal charge −2
   + backbone-O carve  │
   first-shell 3.0 Å   │
        │              │                       │
        └──────────────┼──────── identical from here ────┘
                       │
                       ▼
       <stem>_qm/    {La,Ca,apo}_qm.xyz + bulk_water.xyz
                     submit_<stem>.sh, sp_<stem>_{La,Ca,apo,water}.inp
                       │
                       ▼
                  [ sbatch submit_<stem>.sh ]
                       │
                       ▼  4 ORCA r²SCAN-3c CPCM(Water) SPs, sequential
                  sp_<stem>_{La,Ca,apo,water}.out
                       │
                       ▼  auto-cleanup on success (~30 MB → ~200 KB)
                       │
                       ▼
       [ scripts/update_results_jsonl.py . ]
                       │
                       ▼
       results/all_results.jsonl                ← append-only master record
       results/discriminator_panel_LATEST.tsv   ← TSV view, sorted by ΔΔE desc
       results/ln_class_hits.tsv                ← Ln-class tracker, Sharur-enriched
                       │
                       ▼
       [ scripts/rebuild_ln_tracker.py ]        ← auto-fired from update_results_jsonl
```

A single CIF travels through this in 20–60 min wallclock (4 sequential SPs on
a 64-thread memory-partition node).

---

## 2. Per-script reference

The production scripts live in [`scripts/`](scripts/). Each entry below covers
inputs, outputs, key parameters, and known failure modes. Compact tables are
in [scripts/README.md](scripts/README.md).

### 2.1 `process_inbox.sh`

**Top-level driver.** Reads every `*.cif` in `inbox/`, decides PQQ-vs-generic,
runs the appropriate carve pipeline, writes a SLURM submit script, and submits
it. Idempotent: skips stems whose `sp_<stem>_La.out` already shows a converged
SP, and consumed CIFs move to `inbox/processed/`.

- **Inputs:** any number of CIFs in `inbox/`.
- **Outputs:** `<stem>_qm/` workspace per CIF; `[OK] <stem> -> job <jid>` lines
  to stdout; CIFs moved to `inbox/processed/`.
- **Key parameters:**
  - `MAX_SUBMIT=N` (default 100) — per-invocation submit cap; prevents
    QOSMaxSubmit hits on huge batches.
  - `INBOX`, `PROCESSED`, `FEP_PYTHON`, `SCAN_PYTHON` — hard-coded near the top.
- **Failure modes:**
  - `[SKIP] <stem>: no La in CIF` — the regex `(^| )LA |HETATM[[:space:]]+[0-9]+[[:space:]]+La` didn't match. Either a true apo CIF, or an unusual atom-naming convention. Audit with `grep -i la <cif>`.
  - `[FAIL] normalize` / `[FAIL] protonation` / `[FAIL] carve_with_pqq` / `[FAIL] carve` — one of the sub-steps printed a non-zero exit. The workspace dir is left in place; inspect the partial files.
  - `[FAIL] sbatch: ...` — usually QOS limit or scheduler hiccup. Run again.
- **Concurrency:** flock-guarded on `/tmp/process_inbox.lock`. Concurrent
  invocations no-op until the first finishes.

### 2.2 `carve_generic.py`

**Generic QM cluster carver.** Reads a protonated PDB (or CIF), picks the
La/Ca/Ce/Y site with the most carboxylate donors (or honours `--site-chain`
/ `--site-resnum`), enumerates first-shell donors, and writes three XYZ files
(`{stem}_{La,Ca,apo}_qm.xyz`) plus the SLURM submit script.

- **Inputs:** structure path; output dir; `--stem` (filename prefix); optional
  `--site-chain`, `--site-resnum`, `--first-shell-cut`.
- **Outputs:** `{stem}_{La,Ca,apo}_qm.xyz`, `sp_{stem}_{La,Ca,apo}.inp`,
  `submit_{stem}.sh` (executable). The submit script bakes in the 4-SP loop
  (La → Ca → apo → water) and the success-conditional auto-cleanup tail.
- **Key parameters (in source):**
  - `FIRST_SHELL_CUT = 3.0` Å
  - `SIDECHAIN_QM_ATOMS = {ASP, GLU, ASN, GLN, SER, TYR}` (see [METHODS.md §3](METHODS.md))
  - `ORCA_PATH` — `/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg`
- **Failure modes:**
  - `RuntimeError: No La/Y/Ca/Ce in <path>` — input has no recognised metal
    atom. Check the CIF was normalised first (see `normalize_af3_cif.py`).
  - **Empty carve (1 atom)** — closest donor outside 3.0 Å, or coordinated by
    an unsupported residue (CYS/HIS/MET/LYS/ARG). Recovered via the
    `process_recarve_queue.py` driver.
  - **WARN: residue X not in sidechain dict** — informational; the residue's
    sidechain is dropped, but its backbone O is still carved if within cutoff.

### 2.3 `carve_with_pqq.py`

**PQQ-aware variant.** Same logic as `carve_generic` but includes the
24-atom PQQ cofactor in the QM region at formal charge −2.

- **Inputs:** protonated PDB (with a residue named `PQQ` if PQQ is present);
  output dir; `--stem`, `--metal-chain`, `--metal-resname`. Optional
  `--replace-metal` for Ce-modelled X-ray structures.
- **Outputs:** same as `carve_generic` plus 24 PQQ atoms in each metal carve.
  Typical 55–58 total atoms.
- **Key parameters:**
  - `FIRST_SHELL_CUT = 3.2` Å in source (slightly looser than `carve_generic`'s
    3.0 Å — historical; PQQ-bidentate distances sit near the cutoff).
  - `PQQ_CHARGE = -2` — see [METHODS.md §4](METHODS.md).
  - Sidechain dict matches `carve_generic` (no TYR / backbone-O; PQQ-MDH/ADH
    sites haven't needed those yet).
- **Failure modes:**
  - `No metal <resname> in chain <chain>` — wrong `--metal-chain` or the
    normalisation step didn't rename `LIG_*` → `LA`.
  - **Severe asymmetric carve** — Ca and La carves end up with different
    PQQ inclusion. Diff atom counts via the recipe in [TIERS.md](TIERS.md#severe-asymmetric-carve-manual-flag-no-auto-label).

### 2.4 `normalize_af3_cif.py`

**AF3 / Protenix CIF renamer.** Walks the structure and renames every
`LIG_*` residue:

- residue contains a La atom → renamed to `LA`; atom name `LA1` → `LA`.
- residue has PQQ-like composition (≥14 C + ≥6 O + ≥1 N) → renamed to `PQQ`.

- **Inputs:** source CIF path; output PDB path.
- **Outputs:** PDB file with normalised residue names. Heavy atoms only —
  PDBFixer adds the protons in the next step.
- **When called:** only on the PQQ branch of `process_inbox.sh`. The generic
  branch skips this step because `carve_generic.py` reads CIFs directly.
- **Failure modes:** silent — if a `LIG_*` residue contains neither La nor a
  PQQ-like composition, it is left alone. Downstream carve will fail to find
  the metal and `process_inbox.sh` will print `[FAIL] carve_with_pqq`.

### 2.5 `protonate_cif.py`

**PDBFixer wrapper.** Reads a CIF or PDB, adds missing residues / atoms,
adds Hs at pH 7, writes a PDB.

- **Inputs:** source path; output PDB path.
- **Outputs:** PDB with explicit hydrogens; intermediate temp PDB is removed.
- **Env requirement:** PDBFixer + OpenMM. The pipeline uses
  `/home/jwestrob/miniconda3/envs/fep/bin/python` for this step (see
  [environment.yml](environment.yml) for the `lanm_qmmm` env that drives
  carve + aggregation, and the FEP-env note for PDBFixer).
- **Non-determinism:** PDBFixer chooses tautomer / protonation states using
  internal heuristics that are not perfectly reproducible across machine
  states. Absolute SP energies drift by 15–66 kcal/mol between protonation
  rounds on the same input. **This is fine** — the vertical metal swap
  cancels it. See [METHODS.md §1](METHODS.md).
- **Failure modes:** PDBFixer chokes on certain CIFs (oddly-numbered chains,
  unusual residue mappings). The workspace exists but with no `_protonated.pdb`.
  Inspect `slurm_*.err` (when run inline through the inbox) or rerun the
  command standalone with verbose Python tracing.

### 2.6 `update_results_jsonl.py`

**Aggregator + classifier.** Walks `*_qm/` directories under the workspace,
parses the four `sp_<stem>_*.out` files, computes ΔΔE, classifies, and
appends a JSON object per candidate to `results/all_results.jsonl`. After
walking, writes `results/discriminator_panel_LATEST.tsv` and shells out to
`rebuild_ln_tracker.py` for the Ln-class TSV.

- **Inputs:** `workspace` (positional); `--rebuild` (optional, wipe and
  re-process).
- **Outputs:**
  - `results/all_results.jsonl` — append-only, one JSON per candidate. Schema
    in [HOWTO.md §5](HOWTO.md).
  - `results/discriminator_panel_LATEST.tsv` — TSV view, sorted by ΔΔE desc.
  - `results/ln_class_hits.tsv` — auto-rebuilt via `rebuild_ln_tracker.py`.
- **Key constants:** `E_AQUO_CA`, `E_AQUO_LA`, `DIFF_AQUO`, `HA2KCAL` (see
  [METHODS.md §2](METHODS.md)). Class thresholds in `classify()` (see
  [METHODS.md §6](METHODS.md)).
- **Idempotency:** for each `<stem>_qm/`, the row in JSONL is considered
  fresh if `completed_at` ≥ the dir's mtime. Stale rows are re-extracted and
  re-appended (the *latest* row per stem wins on read).
- **Failure modes:** ORCA `.out` parse failures emit `[warn] <stem>: <exc>` to
  stderr and skip the row. Common cause: SP died mid-run and the file has no
  `FINAL SINGLE POINT ENERGY` line. Re-submit and re-run the aggregator.

### 2.7 `rebuild_ln_tracker.py`

**Ln-class tracker + Sharur DB enrichment.** Reads
`results/all_results.jsonl`, filters to `Ln-preferring` and `Ln-evolved`
classes, joins against the spicy_lams Sharur DB classifications (Pfam, KO,
organism bin, novelty, gene id, length) and the FP catalogue, and writes
`results/ln_class_hits.tsv`.

- **Inputs:** `results/all_results.jsonl`, the two Sharur classifications
  TSVs at `/groups/banfield/users/jwestrob/bin/Sharur/data/spicy_lams/...`,
  and `results/known_false_positives.tsv`.
- **Outputs:** `results/ln_class_hits.tsv` — 17 columns including `rank`,
  `stem`, `source`, `iptm`, `ddE_kcal`, `class`, plus Sharur metadata.
- **Auto-invocation:** `update_results_jsonl.py` shells out to this script at
  the end of every run.
- **Failure modes:** missing Sharur TSV emits a stderr warning and continues
  with `-` placeholders in the enrichment columns. Genuinely standalone — no
  network calls.

### 2.8 `process_recarve_queue.py`

**Recarve driver for historical empty-carves.** Reads
`results/recarve_queue.tsv`, for each entry verifies a protonated PDB exists,
moves stale carve outputs to an audit subdir, picks PQQ vs generic carver,
retries at 3.0 Å (and 3.5 Å if still 1-atom), and resubmits SLURM jobs for
successful carves.

- **Inputs:** `results/recarve_queue.tsv`; optional `--dry-run`.
- **Outputs:** per-stem outcomes logged to `results/recarve_log_<UTC>.tsv` and
  to stdout. Outcome codes:
  - `RECARVED_AND_SUBMITTED` — multi-atom carve produced; SLURM job in queue.
  - `CARVE_AMBIGUOUS` — closest donor is unsupported (CYS/HIS/MET/LYS/ARG).
  - `SOLVENT_EXCLUDE` — closest donor > 4 Å — pocket is genuinely empty.
- **Key parameters:**
  - `SOLVENT_EXCLUDE_THRESHOLD = 4.0` Å
  - `PQQ_DETECT_RADIUS = 4.0` Å — any PQQ atom within this routes to `carve_with_pqq`
  - `RETRY_WIDE_CUTOFF = 3.5` Å
- **Failure modes:** an entry without a `<stem>_protonated.pdb` is logged as
  `NO_PROTONATED_PDB` and skipped.

### 2.9 `inbox_watcher.sh`

**Long-running poll loop.** Every 45 min: if `inbox/*.cif` is non-empty and
the SLURM queue has free slots (target ≤ 195 of 200 user limit), runs
`process_inbox.sh`. Lifetime tied to `/tmp/inbox_watcher.lock`.

- **Start:** `bash scripts/inbox_watcher.sh > watcher.log 2>&1 &`; the PID is
  written to the lockfile so the next iteration can sanity-check it.
- **Stop:** `rm /tmp/inbox_watcher.lock` — the loop exits at its next iteration.
- **Failure modes:** none in practice. If `process_inbox.sh` errors, the watcher
  just logs and continues.

---

## 3. Idempotency contract

| script                            | safe to rerun?                                                                                             |
|-----------------------------------|------------------------------------------------------------------------------------------------------------|
| `process_inbox.sh`                | yes — flock-guarded, skips stems with a converged La SP                                                    |
| `carve_generic.py` / `_with_pqq` | yes — overwrites XYZ + inp + submit; safe to call again to refresh                                         |
| `normalize_af3_cif.py`            | yes — pure file rewrite                                                                                    |
| `protonate_cif.py`                | yes — but each run picks a different protonation (PDBFixer non-determinism)                                |
| `update_results_jsonl.py`         | yes — designed idempotent; mtime-based skip, append-only JSONL, deduped by stem on read                    |
| `rebuild_ln_tracker.py`           | yes — overwrites `ln_class_hits.tsv` from scratch                                                          |
| `process_recarve_queue.py`        | yes — won't double-submit if `squeue` shows the stem in flight, won't re-carve an already-completed La SP  |
| `inbox_watcher.sh`                | only one at a time — lockfile-PID convention                                                               |

**JSONL append semantics.** `all_results.jsonl` is append-only: the
aggregator never rewrites old rows. On read, the *latest* row per `stem` wins
(via `load_existing(jsonl)` which builds a `dict[stem] = row` from a
linear scan). This means stale rows pile up in the file, but consumers
always see the freshest classification. Compact the JSONL only via
`--rebuild` if it gets unwieldy (currently ~2 MB at 1,453 entries; not
worth compacting).

---

## 4. Failure modes catalogue

Cross-referenced with [HOWTO.md §11](HOWTO.md). The compact version:

| symptom                                                              | most likely cause                                                | fix                                                                                       |
|----------------------------------------------------------------------|------------------------------------------------------------------|-------------------------------------------------------------------------------------------|
| `E[La] = −31.158775581143` exact, `n_atoms = 1`                      | empty carve — no eligible donor inside 3.0 Å                     | Add stem to `recarve_queue.tsv`, run `process_recarve_queue.py`                           |
| FAILED job in < 10 s                                                 | duplicate submission (legacy) or bad submit script               | Now mostly prevented by flock; check `slurm_*.err`                                        |
| `post_scan_fold` / `*_la_fold` / array task FAILED                   | upstream AlphaFold/Protenix folding (different agent's pipeline) | Not the discriminator. Filter `squeue` output with `grep -v fold`                         |
| Protonation FAIL, no `_protonated.pdb` written                       | PDBFixer chokes on odd CIF                                       | Inspect `slurm_*.err`; rerun standalone with verbose tracing                              |
| `\|ΔΔE\| > 30` on a small acidic pocket                              | OUTLIER tier — Mn²⁺/Zn²⁺ pocket                                  | Read as biology signal, not error; cross-ref Pfam/KO                                      |
| `\|ΔΔE\| > 100`                                                      | SEVERE asymmetric carve — different ligand sets in Ca/La         | Diff `head -1 <stem>_qm/<stem>_{Ca,La}_qm.xyz`; manual investigation                      |
| SCF didn't converge (no `FINAL SINGLE POINT ENERGY` line)            | Bad initial guess on La / level mismatch                         | Retry with `TightSCF` + `%scf shift shift 0.4 erroff 0.1 end end` ([HOWTO.md §11](HOWTO.md)) |
| Workspace dir at ~30 MB after job completion                         | Auto-cleanup didn't fire — one SP did not return "TERMINATED NORMALLY" | `bash scripts/cleanup_orca_scratch.sh` then inspect the failed SP                          |
| Disk filling on shared volume                                        | Many in-flight jobs, scratch accumulating                        | Same cleanup script; consider lowering `MAX_SUBMIT`                                       |

---

## 5. Cluster integration

### SLURM submit template

Every `submit_<stem>.sh` written by the carvers follows the same skeleton:

```bash
#!/bin/bash
#SBATCH -p memory
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J orca_<stem>
#SBATCH -o <out_dir>/slurm_%j.out
#SBATCH -e <out_dir>/slurm_%j.err

ORCA_PATH=/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export PATH=$ORCA_PATH:$PATH
export LD_LIBRARY_PATH=$ORCA_PATH:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
export MKL_NUM_THREADS=$SLURM_CPUS_ON_NODE

cd <out_dir>
export ORCA_TMPDIR=$PWD; export TMPDIR=$PWD
for kind in La Ca apo water; do
    $ORCA_PATH/orca sp_<stem>_$kind.inp > sp_<stem>_$kind.out 2>&1
done
# success-conditional auto-cleanup tail (see carve_generic.py)
```

Key points:

- **`-p memory` partition.** Memory-partition nodes have 64 cores and enough
  RAM (`%maxcore 8000` per thread × 64) for the larger PQQ-aware carves.
- **`--exclusive` and `-N 1`.** Full-node booking; no node-sharing.
- **OpenMP threading via `OMP_NUM_THREADS`.** ORCA's r²SCAN-3c implementation
  uses MKL + OpenMP. **Never use MPI / PMIX** with this ORCA build — the
  shared-openmpi418 binary will segfault under mpirun on these clusters.
- **`ORCA_TMPDIR=$PWD`.** Keeps scratch on local disk under the workspace
  (rather than `/tmp`), which prevents I/O contention when many jobs run
  concurrently.

### ORCA environment

```
ORCA: /home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
```

ORCA 6.1.1 with the shared-OpenMPI 4.1.8 build. Required: `PATH`,
`LD_LIBRARY_PATH` both include the ORCA bin dir.

### Conda envs

Two envs, both already provisioned on the cluster:

- **`lanm_qmmm`** — `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/`.
  Used by every script except `protonate_cif.py`. Contains gemmi, numpy,
  biopython, pandas. Python 3.11+.
- **`fep`** — `/home/jwestrob/miniconda3/envs/fep`. Used **only** for
  `protonate_cif.py` because PDBFixer needs OpenMM and a Python 3.10 build.

Spec: [environment.yml](environment.yml) for `lanm_qmmm`, plus a
`requirements_fep.txt` note for PDBFixer.

---

## 6. State files

```
/tmp/process_inbox.lock        — flock-guarded mutex for process_inbox.sh
/tmp/inbox_watcher.lock        — PID file + sentinel for inbox_watcher.sh

inbox/                         — drop zone
inbox/processed/               — CIFs that got submitted
inbox/apo_skipped/             — CIFs that filtered out as "no La"
inbox/stale_pre_cn_filter/     — legacy pre-prefilter rejects (audit)

<stem>_qm/                     — per-candidate workspace
results/all_results.jsonl      — master append-only JSONL
results/discriminator_panel_LATEST.tsv  — TSV view
results/ln_class_hits.tsv      — Ln-class tracker
results/known_false_positives.tsv       — manual FP catalogue
results/confirmed_ln_binders.tsv         — manual confirmed-binder catalogue
results/recarve_queue.tsv      — workqueue for empty-carve recovery
results/recarve_log_<UTC>.tsv  — per-stem outcomes from process_recarve_queue
```

Lockfile-PID convention for `inbox_watcher`: the daemon writes `$$` (its
PID) to `/tmp/inbox_watcher.lock` at start and checks at every iteration
that the lockfile still exists *and* still contains its own PID. Either
deleting the lockfile or starting a second watcher (which overwrites the
PID) cleanly terminates the older instance at the next 45-min tick.

---

## 7. Pointers

- One-page user how-to: [HOWTO.md](HOWTO.md)
- Methodology rationale: [METHODS.md](METHODS.md)
- Class label semantics: [TIERS.md](TIERS.md)
- Validation evidence: [VALIDATION.md](VALIDATION.md)
- Project state: [CONTEXT.md](CONTEXT.md)
- Per-script compact reference: [scripts/README.md](scripts/README.md)
