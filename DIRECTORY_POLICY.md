# Directory Layout Policy

**Enforce this layout. The repo root is for project-level artifacts only.
Per-candidate compute output belongs under `workspaces/`, always.**

Established 2026-07-16 after cleanup of 1757 `*_qm/` directories that had
accumulated at the top level. Do not undo this.

## The rule

- **`<repo_root>/<stem>_qm/`** — FORBIDDEN. Never write per-candidate workspaces at the top level.
- **`<repo_root>/workspaces/<stem>_qm/`** — REQUIRED location for all per-candidate ORCA workspaces.

Any script that creates a new candidate workspace MUST write it under `workspaces/`. This includes but is not limited to:

- `scripts/process_inbox.sh` (line 38: `out_dir="$ALCH/workspaces/${stem}_qm"`)
- `scripts/carve_generic.py` (called with `out_dir` argument — caller's responsibility)
- `scripts/carve_with_pqq.py` (same)
- `scripts/process_recarve_queue.py`
- Any new pipeline entry point

Any script that scans for existing workspaces MUST use `workspaces/*_qm` (not `*_qm`). Currently correct:

- `scripts/analyze_panel.py`
- `scripts/rebuild_recarve_queue.py`
- `scripts/update_results_jsonl.py`
- `scripts/cleanup_orca_scratch.sh`

## Repo root layout (what belongs where)

```
alchemical_bvs/
├── README.md                # public entry point
├── CHANGELOG.md             # dated methodology / pipeline changes
├── CONTEXT.md               # active project state
├── HANDOFF.md               # legacy handoff (pre-2026-07-16)
├── HOWTO.md                 # operator manual
├── METHODS.md               # full methodology (paper 2 material)
├── PIPELINE.md              # data flow diagram
├── RESULTS.md               # publishable findings
├── SESSIONS.md              # cross-agent coordination log
├── TIERS.md                 # tier assignment definitions
├── VALIDATION.md            # benchmark validation
├── REPORT.md                # historical pre-pivot report
├── CONTRIBUTING.md          # contribution guidelines
├── CITATION.cff             # citation metadata
├── LICENSE                  # (TODO — placeholder)
├── DIRECTORY_POLICY.md      # ← this file
├── environment.yml          # conda env spec
├── requirements_fep.txt     # legacy FEP-track requirements
├── docs/                    # supplementary docs (e.g., upstream_fold_daemon.md)
├── docs_buildout_brief.md   # doc-build notes
├── examples/                # walked-through example
├── inbox/                   # incoming CIFs for processing
│   ├── processed/           # audit trail of consumed inputs
│   ├── apo_skipped/         # skipped apo-cofactor entries
│   └── stale_pre_cn_filter/ # pre-filter-era entries
├── inputs/                  # curated input CIFs (non-inbox)
├── params/                  # ORCA + protocol parameters
├── paper_methods/           # Colin's-paper handoff (2026-07-16)
│   ├── colin_methods_paragraph.md
│   ├── citations.bib
│   └── paper2_pitch.md
├── protenix_jobs/           # Protenix fold job definitions
├── results/                 # published results (in git)
│   ├── all_results.jsonl    # main manifest
│   ├── ln_class_hits.tsv
│   ├── confirmed_ln_binders.tsv
│   ├── discriminator_panel_LATEST.tsv
│   ├── pqq_refold_queue.tsv
│   ├── recarve_queue.tsv
│   └── recarve_log_*.tsv
├── scripts/                 # all pipeline code
├── colin_handoff/           # Colin-specific inputs/outputs
├── workspaces/              # ★ per-candidate ORCA workspaces (gitignored)
│   └── <stem>_qm/           # one dir per candidate; created by process_inbox.sh
└── legacy/                  # ★ pre-pivot compute trees + completed validation panels (gitignored)
    ├── qmmm/                # pre-pivot QM/MM-MD track (~36 GB)
    ├── qmmm_lc/             # (~2.3 GB)
    ├── qmmm_md/             # smaller
    ├── b97_3c_panel/        # functional-robustness validation (~116 MB, results in VALIDATION.md)
    ├── calexcitin_size_panel/  # cluster-size sensitivity (~130 MB, results in VALIDATION.md)
    └── spicy_lams_low_cn_staging/  # old staging dir
```

## Rationale

Before the 2026-07-16 cleanup, 1757 `*_qm/` workspaces sat at the top level alongside project docs, scripts, and results. Consequences:

- `ls alchemical_bvs/` returned 1800+ entries, hiding legitimate project files
- New agents opening the repo couldn't find their bearings
- `find` and `git status` were slow to scan
- Auto-completion in shells was useless
- Even though workspaces were gitignored, they polluted the mental model of the repo

After cleanup, top-level `ls` returns ~35 entries, all of which are meaningful project artifacts.

## If you break this policy

If someone (agent or human) writes a `<stem>_qm/` workspace at the top level, they should:

1. Move it: `mv <stem>_qm workspaces/`
2. Update the `cluster_panel_path` field in `results/all_results.jsonl` for that stem to point to `workspaces/<stem>_qm`
3. Refresh derived tables: `python scripts/update_results_jsonl.py .`
4. Add a note in `SESSIONS.md` explaining what tool violated the policy so it can be fixed

## Where this policy lives

- **This file** (`DIRECTORY_POLICY.md`) at the repo root — canonical.
- **`.gitignore`** — enforces via `workspaces/` blanket ignore + safety-net entries for `legacy/` and its sub-dirs.
- **`CONTEXT.md`** — cross-linked in project state.
- **`SESSIONS.md`** — cross-linked in cross-agent coordination.
- **`README.md`** — should link to this file in the "developer notes" section.
