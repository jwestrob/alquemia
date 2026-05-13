# Changelog

All notable methodology and pipeline changes to the DFT vertical-swap Ca²⁺/Ln³⁺
discriminator (`alchemical_bvs`) are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/) with semantic-style
date headers; the project is methodology, so changes are described in terms
of what they alter about the numbers produced, not API surface.

For project state (panel counts, hit families, open work) see [CONTEXT.md](CONTEXT.md).
For how to run the pipeline see [HOWTO.md](HOWTO.md).

---

## [Unreleased]

### Added
- Public-repo documentation buildout (`README.md`, `METHODS.md`, `PIPELINE.md`,
  `TIERS.md`, `VALIDATION.md`, `RESULTS.md`, `CONTRIBUTING.md`, `CITATION.cff`,
  `environment.yml`, `scripts/README.md`, `examples/walkthrough.md`).

### Pending
- Choose and apply an OSS license (`LICENSE` is currently a placeholder).
- Fold-daemon upstream `carve_generic.py` path — still emitting empty-carves on
  some new CIFs because it references a stale carver. Diagnostic queue accumulates
  in `results/recarve_queue.tsv`.

---

## 2026-05-12 — Aggregator auto-rebuild chain validated; Eukaryotic PQQ-ADH discovery

### Validated
- `update_results_jsonl.py` now auto-calls `rebuild_ln_tracker.py` on every run.
  `results/ln_class_hits.tsv` is regenerated from the JSONL plus the Sharur DB
  classifications and `known_false_positives.tsv`. End-to-end chain
  (drop CIF → 4 ORCA SPs → JSONL → Ln tracker TSV) confirmed at scale (1,453
  candidates).

### Discovered
- Of the 50 eukaryotic PQQ-ADH structures from Colin's `03_Euk` batch processed
  to date, 21 score Ln-preferring or Ln-evolved (≥+5 kcal/mol). Top hits
  A0A2V0PRL1 domain3 (+23.45) and A0A2V0NW76 domain1 (+22.35) sit in the
  Ln-evolved tier. Four multi-domain proteins (A0A1R3GAW4, A0A5N6PFI5,
  A0A6J5WV68, A0A7J7PMX4) display intra-protein cross-domain
  differentiation — paralogous PQQ-binding domains of the same chain score in
  different tiers. Recorded in `results/colin_euk_pqq_adh_results.tsv`.

### Numbers (live)
- 1,453 candidates scored / 17 Ln-evolved / 278 Ln-preferring / 96 OUTLIER /
  90 EXCLUDE (deduped by stem) / 137 pending.

---

## 2026-05-11 PM — PQQ-aware inbox routing

### Added
- `scripts/process_inbox.sh` now inspects each CIF for a PQQ-like ligand:
  a residue named `LIG_*` with composition ≥12 C + ≥6 O + ≥1 N
  (canonical PQQ = 14 C + 8 O + 2 N, 24 heavy atoms). If detected, the CIF is
  routed through `normalize_af3_cif.py` → `protonate_cif.py` → `carve_with_pqq.py`;
  otherwise through the legacy `protonate_cif.py` → `carve_generic.py` path.
- `normalize_af3_cif.py`: renames Protenix/AF3 `LIG_*` residues to `LA` (if a La
  atom is present) or `PQQ` (if the composition matches), so downstream
  carvers can pick them up.

### Changed
- Relaxed the `process_inbox.sh` "has-La" filter from a bare " LA " grep to
  also match `HETATM   ...   La   ...` records and `LA1` residue names — Colin's
  240-CIF drop (`00_Lav`, `01_Cav`, `03_Euk`) was filter-rejected as "no La in CIF"
  under the old regex because Protenix writes atom name `LA1` instead of `LA`.

### Why it matters
- PQQ is the primary Ln coordinator in PQQ-MDH/ADH active sites (PQQ-O5 +
  PQQ-O4 are bidentate to the metal). A protein-only carve of a PQQ-binder
  systematically underestimates ΔΔE by ~10 kcal/mol; the routing fix
  reproduces Colin's 00_Lav set within 0.2–3 kcal/mol of the legacy direct-driver
  runs.

---

## 2026-05-11 (14:25) — Recarve queue

### Added
- `scripts/process_recarve_queue.py`: re-carve and re-submit driver for
  historical empty-carve entries (those with the `E[La] = −31.158775581143`
  signature, 1 atom in the QM file). PQQ-aware retry, 3.5 Å widening fallback,
  classification of still-failing entries as `SOLVENT_EXCLUDE` (closest donor
  > 4 Å — genuinely empty pocket) or `CARVE_AMBIGUOUS` (closest donor is an
  unsupported residue, e.g. CYS/HIS/MET).
- `results/recarve_queue.tsv` workqueue and `results/recarve_log_<UTC>.tsv`
  per-stem outcome logs.

### Reprocessed
- 19 historical empty-carve entries replayed through the new carve dict; the
  majority now produce multi-atom carves.

---

## 2026-05-11 (01:10) — Carve-dict tightening; TYR and backbone-O carving

### Changed
- `scripts/carve_generic.py`: `FIRST_SHELL_CUT` tightened from **3.2 Å → 3.0 Å**.
  Tighter cutoff aligned with physically meaningful first-shell distances on
  Ln/Ca clusters and dropped a number of false-near-shell donors.
- `SIDECHAIN_QM_ATOMS`: **TYR** added (full phenol ring + OH). Captures Tyr-OH
  coordination, which was previously a silent empty-carve failure mode.
- Backbone-O carving: when a residue's backbone carbonyl O is within
  `FIRST_SHELL_CUT` of the metal, the C=O group is carved into the QM region
  with two link-H caps (toward Cα and toward the next N). Captures GLY/PRO and
  other non-sidechain-coordinated sites.
- Donor scan now considers S (Cys SG, Met SD) — was O/N only, which silently
  dropped Cys-coordinated sites to bare-La carves.

### Policy
- An agent attempted to add THR/CYS/HIS/MET/LYS/ARG to `SIDECHAIN_QM_ATOMS`;
  those additions were **reverted out** per Jacob's explicit direction. Only
  TYR was authorized in addition to the original ASP/GLU/ASN/GLN/SER set.
  If a re-carve still produces a 1-atom cluster after the TYR + backbone-O
  fix, the entry is classified as `SOLVENT_EXCLUDE` or `CARVE_AMBIGUOUS` —
  the dict is not expanded. This restriction is intentional: the discriminator
  is calibrated on hard-O carboxylate/amide/hydroxyl donors. Adding soft
  S-donors or N-aromatic donors would change what is being measured.

---

## 2026-05-10 — Carve auto-cleanup

### Added
- `scripts/carve_generic.py` now appends an auto-cleanup tail to every SLURM
  submit script. When all 4 SPs return "ORCA TERMINATED NORMALLY", the script
  deletes ORCA scratch files (`.gbw`, `.tmp`, `.bas?`, etc.) and slurm logs in
  place. Workspace dirs drop from ~30 MB to ~200 KB on success. The bulk-water
  SP is also now baked into the default sequential loop instead of needing a
  post-hoc `sed` patch.

---

## 2026-05-07 — Production discriminator panel v1 lock

### Locked
- ΔΔE class thresholds: ≥+20 (Ln-evolved), +5..+20 (Ln-preferring),
  0..+5 (marginal), −5..0 (ambiguous), <−5 (Ca-evolved),
  |ΔΔE| > 30 (OUTLIER — typically Mn²⁺/Zn²⁺ pocket signature).
- DIFF_AQUO (Ca·6H₂O − La·6H₂O at r²SCAN-3c/CPCM(Water)) =
  **−646.064555551745 Ha**. Frozen reference.
- Electronic structure: r²SCAN-3c composite functional, CPCM(Water),
  DefGrid3, NoAutostart. La: def2-ECP + def2-TZVP large-core (4f^n in core,
  closed-shell singlet). Ca: all-electron. Multiplicity 1 throughout.
- The first production panel covered 54 sites including the XoxF/MxaF
  gold-standard pair (+24.1 vs −6.5), LanM EF1 (+47.5), tannase (+15.7), and
  the Ca-class negative controls (calmodulin +18.2, parvalbumin, calbindin
  +3.7). Recorded in `results/discriminator_panel_2026-05-07.tsv` and
  `results/findings_summary_2026-05-07.md`.

### Validated
- B97-3c rerun on the top 4 (XoxF, MxaF, tannase, calexcitin) gives the same
  sign in every case; magnitudes scatter <8 kcal/mol — class assignment
  preserved. See `b97_3c_panel/`.
- Cluster-size sensitivity on calexcitin (16/32/50/100 atoms) — the carve
  composition matters more than total atom count; the carbox-cluster
  contribution converges by ~30 atoms when the dict captures the right donors.
  See `calexcitin_size_panel/`.

---

## Earlier (2026-04 / pre-discriminator)

The work prior to the 2026-05-07 panel lock is summarized in `HANDOFF.md`
and `REPORT.md` at the repo root, and in `~/jwestrob/obsidian-vault/projects/lanthanide-binding/CONTEXT.md`. Highlights:

- 2026-03..04: ORCA pipeline diagnosis (serial-mode test), aquo Ln panel,
  EF1 cluster optimisations, AmberTools build, QM/MM single-point ORCA
  inputs, P3 vs P3v2 large-core ECP rerun.
- 2026-04..05-03: 8DQ2 X-ray pre-check (alchemical BVS Pass 2 residuals)
  and BVS-from-density experiments that motivated pivoting to a vertical
  Ca/Ln swap discriminator instead of refining BVS scaling parameters.
