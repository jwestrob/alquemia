# DFT Ca²⁺/Ln³⁺ Discriminator — How-To

**Purpose:** Quantitatively decide whether a predicted/known La³⁺-binding pocket actually
prefers a lanthanide over Ca²⁺, using a vertical metal-swap DFT cycle.

**Method, in one line:** ΔΔE = (E[Ca-cluster] − E[La-cluster]) − DIFF_AQUO. Positive = pocket prefers Ln³⁺. Negative = pocket prefers Ca²⁺. The aquo reference cancels absolute-energy noise (PDBFixer non-determinism, basis-set artifacts, ECP-La vs all-electron-Ca offset).

**Scale (2026-05-12):** 1,453 candidates scored, 17 Ln-evolved + 277 Ln-preferring hits.

---

## 1. Quick start (for impatient agents/humans)

You have a CIF/PDB with **La** placed in a candidate metal-binding pocket. The pipeline auto-handles: protonation → carve → 4× ORCA single-points → ΔΔE → JSONL aggregation → class assignment.

```bash
# 1. Drop the CIF in the inbox
cp my_protein.cif /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/inbox/

# 2. Trigger the pipeline (or wait 45 min for the inbox watcher)
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
bash scripts/process_inbox.sh

# 3. Check results when SLURM finishes (~20–60 min per candidate)
python scripts/update_results_jsonl.py .
column -t -s $'\t' results/discriminator_panel_LATEST.tsv | tail -20
```

The inbox CIF must contain La. PDBFixer protonation and PQQ detection happen automatically.

---

## 2. Architecture

```
inbox/<stem>.cif
       │
       ▼
process_inbox.sh ─── PQQ-aware routing branch ─────────────────────────────┐
       │                                                                   │
       │ "no PQQ" branch                            "has PQQ" branch       │
       │ (LIG_* with <12C or <6O or 0N)             (≥12C + ≥6O + ≥1N)     │
       │                                                                   │
       ▼                                                                   ▼
  protonate_cif.py                                                normalize_af3_cif.py
       │                                                          (LIG_* → LA / PQQ)
       │                                                                   │
       │                                                                   ▼
       │                                                          protonate_cif.py
       │                                                                   │
       ▼                                                                   ▼
  carve_generic.py                                                  carve_with_pqq.py
       │ ASP/GLU/ASN/GLN/SER/TYR sidechains                          │ same +PQQ (24 atoms)
       │ backbone-O carving                                          │ (cofactor in QM region)
       │ r_first_shell = 3.0 Å                                       │
       │                                                              │
       ▼                                                              ▼
  La/Ca/apo XYZ + submit_<stem>.sh  ←─────── identical from here onward
       │
       ▼
  sbatch submit_<stem>.sh
       │
       ▼  (4 ORCA SPs: La, Ca, apo, water — sequential)
  sp_<stem>_{La,Ca,apo,water}.out
       │
       ▼  (auto-cleanup on success, ~30 MB → ~200 KB)
       │
       ▼
  scripts/update_results_jsonl.py
       │
       ▼
  results/all_results.jsonl          (append-only master record)
  results/discriminator_panel_LATEST.tsv    (auto-regenerated)
  results/ln_class_hits.tsv          (Ln-class tracker, Sharur-DB-enriched)
```

### Key constants (locked by validation)

```python
DIFF_AQUO   = -646.064555551745  # Ha,  E[Ca·6H₂O] − E[La·6H₂O] at r²SCAN-3c/CPCM(Water)
HA2KCAL     = 627.5095
FIRST_SHELL_CUT = 3.0            # Å — first-shell donor cutoff (carve_generic.py)
```

`DIFF_AQUO` is your reference. It came from CPCM single-points on `[Ca(H₂O)₆]²⁺` and `[La(H₂O)₆]³⁺` at the **same level of theory** as every workspace SP. Don't change without re-running aquo references.

---

## 3. Class thresholds

```
ΔΔE (kcal/mol)      class                       interpretation
≥ +20               Ln-evolved                  textbook XoxF / classic Ln-binder
+5  to +20          Ln-preferring               likely real Ln-binder, moderate
0   to +5           marginal                    weakly Ln-preferring, often a paralog
−5  to 0            ambiguous                   essentially neutral, in noise band
< −5                Ca-evolved                  pocket prefers Ca²⁺
|ΔΔE| > 30          OUTLIER (verify SCF)        Mn²⁺/Zn²⁺ small-acidic-pocket signature;
                                                 SCF blowup is biology signal, not noise
|ΔΔE| > 100         SEVERE asymmetric carve    Ca and La cluster have different ligand sets;
                                                 manual investigation needed (atom-count diff)
```

The OUTLIER tier is a **useful annotation, not an error** — it usually means the pocket "wants" a smaller harder cation than La (Mn²⁺/Zn²⁺/Fe³⁺) and the discriminator can't accommodate it geometrically. Treat as biology signal.

The SEVERE asymmetric carve tier is a **flag for manual investigation** — when ΔE(metal) drifts ~0.18–0.65 Ha off the canonical −646.06, the Ca and La carves typically have different ligand sets (often PQQ inclusion mismatch). Diff their atom counts:
```bash
head -1 <stem>_qm/<stem>_{Ca,La}_qm.xyz
```

---

## 4. PQQ-aware inbox routing (2026-05-11)

`process_inbox.sh` inspects each CIF for **PQQ-like ligands**: a residue with name starting `LIG` containing ≥12 C + ≥6 O + ≥1 N (canonical PQQ has 14C + 8O + 2N = 24 heavy atoms).

If detected:
- **tier1 pipeline**: `normalize_af3_cif.py` (LIG_*/LA1 → LA/PQQ) → `protonate_cif.py` → `carve_with_pqq.py` (PQQ included in QM region as 24 atoms with formal charge −2)
- Typical QM cluster: 55–58 atoms (metal + 3 carboxylates + 24 PQQ atoms + 4 link Hs)

If not detected:
- **generic pipeline**: `protonate_cif.py` → `carve_generic.py` (protein-only)
- Typical QM cluster: 25–35 atoms (metal + 2–3 carboxylates + link Hs)

This matters because **PQQ is the primary Ln-coordinator in PQQ-MDH/ADH active sites** (PQQ-O5 + PQQ-O4 are tight bidentate to the metal). A protein-only carve of a PQQ-binder will systematically underestimate ΔΔE by ~10 kcal/mol — see the false-bad result on Colin's 00_Lav set before this routing was added.

---

## 5. Files & scripts inventory

### `scripts/` — pipeline guts

| script | purpose |
|---|---|
| `process_inbox.sh` | Top-level driver. Reads `inbox/*.cif`. PQQ-aware routing. Flock-guarded. `MAX_SUBMIT=100` per call. |
| `carve_generic.py` | Generic CIF/PDB → 4 metal-swapped XYZ carves + ORCA submit script. Sidechain dict: **ASP, GLU, ASN, GLN, SER, TYR**. Cutoff 3.0 Å. Backbone-O carving. CLI: `--stem`, `--site-chain`, `--site-resnum`, `--first-shell-cut`. |
| `carve_with_pqq.py` | PQQ-aware carver. Includes 24-atom PQQ at formal charge −2. Auto-detects PQQ residue by composition (≥12C + ≥6O + ≥1N) once residue is named PQQ. |
| `normalize_af3_cif.py` | Renames Protenix/AF3 `LIG_*` residues to `LA` (if La-containing) or `PQQ` (if PQQ-like). Output: PDB. Required before `carve_with_pqq.py` on AF3-style CIFs. |
| `protonate_cif.py` | PDBFixer wrapper, pH 7. Accepts CIF or PDB. Writes PDB. |
| `update_results_jsonl.py` | Walk `*_qm/` dirs, parse SP energies, compute ΔΔE, classify, append to `results/all_results.jsonl`, regenerate `discriminator_panel_LATEST.tsv`, auto-rebuild `ln_class_hits.tsv`. Idempotent. **Run this manually after a batch lands** — does not auto-fire on every SP. |
| `rebuild_ln_tracker.py` | Rebuild `ln_class_hits.tsv` from JSONL, enrich with spicy_lams Sharur DB annotations (Pfam, KO, organism bin, novelty). Called automatically by `update_results_jsonl.py`. |
| `process_recarve_queue.py` | Re-carve and resubmit historical empty-carve entries. Reads `results/recarve_queue.tsv`. PQQ-aware retry; 3.5 Å widening fallback; classifies still-failing entries as SOLVENT_EXCLUDE or CARVE_AMBIGUOUS. |
| `inbox_watcher.sh` | Long-running poll loop. Every 45 min, if inbox has CIFs and queue has free slots, calls `process_inbox.sh`. PID lockfile at `/tmp/inbox_watcher.lock`. |
| `tier1_pqq_batch.sh` | Legacy direct driver for Colin's 25 PQQ-MDH controls (pre-PQQ-aware-routing). Kept for reproducibility. |
| `prefilter_spicylams.py` | Donor-count pre-filter. Drops CIFs with <2 ASP/GLU/ASN/GLN/SER/TYR within 3.2 Å of La. Use before staging large batches. |
| `cleanup_orca_scratch.sh` | Walk `*_qm/` dirs, delete ORCA scratch for completed dirs. Idempotent. New jobs auto-clean on success; this is for legacy dirs. |

### `results/` — outputs

| file | content |
|---|---|
| `all_results.jsonl` | Master record. One JSON per candidate. Schema: `stem, source, iptm, n_atoms, ddE_kcal, class, cif_path, workspace_path`, etc. |
| `discriminator_panel_LATEST.tsv` | TSV view of all 1,453+ rows. Auto-regenerated by `update_results_jsonl.py`. |
| `ln_class_hits.tsv` | Curated Ln-class hits (ΔΔE ≥ +5). Auto-enriched with spicy_lams Pfam/KO/bin annotations. |
| `known_false_positives.tsv` | Manually curated FP catalogue with reason notes. |
| `confirmed_ln_binders.tsv` | Manually curated confirmed-real-Ln-binder catalogue with literature refs. |
| `colin_euk_pqq_adh_results.tsv` | The Colin Euk discovery (2026-05-12) — 50+ scored eukaryotic PQQ-ADHs with iPTM annotations. |
| `recarve_queue.tsv` | Workqueue for re-processing empty-carve entries (E[La]=−31.16 signature). |
| `recarve_log_<UTC>.tsv` | Per-stem outcome from `process_recarve_queue.py` (RECARVED_AND_SUBMITTED, CARVE_AMBIGUOUS, SOLVENT_EXCLUDE). |

### `inbox/` — drop zone

```
inbox/
├── *.cif              ← drop new CIFs here
├── processed/         ← consumed CIFs move here after sbatch
├── apo_skipped/       ← CIFs with no La (or filter-rejected) go here
└── stale_pre_cn_filter/  ← legacy pre-filter rejects, kept for audit
```

### Workspace dirs (`<stem>_qm/`)

After cleanup, each contains:

```
<stem>_qm/
├── <stem>_normalized.pdb         ← (PQQ pipeline only) renamed-residue PDB
├── <stem>_protonated.pdb         ← PDBFixer output
├── <stem>_{La,Ca,apo}_qm.xyz     ← carved clusters
├── bulk_water.xyz                ← water aquo reference
├── sp_<stem>_{La,Ca,apo,water}.inp    ← ORCA inputs
├── sp_<stem>_{La,Ca,apo,water}.out    ← ORCA outputs (energies + diagnostics)
└── submit_<stem>.sh              ← SLURM submit script
```

Pre-cleanup ~30 MB (mostly .gbw + .tmp scratch); post-cleanup ~200 KB.

---

## 6. End-to-end walkthroughs

### 6a. Single CIF, foreground

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs

# Stage
cp /path/to/my_protein.cif inbox/

# Process (carves + submits 1 SLURM job per CIF)
bash scripts/process_inbox.sh
# → [OK] my_protein -> job 1234567   (or → "PQQ detected → tier1 pipeline" first)

# Wait for SLURM
squeue -u $USER -j 1234567

# Compute ΔΔE
python scripts/update_results_jsonl.py .

# Read result
grep my_protein results/discriminator_panel_LATEST.tsv
```

### 6b. Bulk batch (hundreds–thousands of CIFs)

```bash
# Pre-filter (drops CIFs with no carve-eligible donors)
python scripts/prefilter_spicylams.py /path/to/staging_dir 2

# Stage with symlinks to save disk
mkdir -p staging
ln -sf /path/to/staging_dir/*.cif staging/

# Move to inbox in waves of ~100 (MAX_SUBMIT cap)
mv staging/CIF_BATCH_1*.cif inbox/

# Submit
bash scripts/process_inbox.sh   # submits up to 100, leaves rest in inbox

# Watcher will pick up the next wave at the next 45-min tick
```

### 6c. Long-running unattended scan

```bash
bash scripts/inbox_watcher.sh > watcher.log 2>&1 &
echo $! > /tmp/inbox_watcher.lock

# Drop CIFs whenever; they get picked up at the next 45-min poll.
# Stop with: rm /tmp/inbox_watcher.lock
```

### 6d. Manifest refresh

After a batch lands, refresh the master record:

```bash
python scripts/update_results_jsonl.py .
# Prints class breakdown, regenerates discriminator_panel_LATEST.tsv,
# auto-calls rebuild_ln_tracker.py for ln_class_hits.tsv.
```

The aggregator is **idempotent** — safe to call repeatedly.

### 6e. Recarve queue (empty-carve recovery)

```bash
# Dry-run
python scripts/process_recarve_queue.py --dry-run

# Live
python scripts/process_recarve_queue.py
# Outcomes:
#   RECARVED_AND_SUBMITTED   new multi-atom carve, jobs running
#   CARVE_AMBIGUOUS          closest donor is unsupported (CYS/HIS/etc.) — manual review
#   SOLVENT_EXCLUDE          closest donor >4 Å — genuinely empty pocket
# Per-stem outcomes logged to results/recarve_log_<UTC>.tsv
```

---

## 7. Reading the panel

```bash
# Class breakdown
python scripts/update_results_jsonl.py . | tail -15

# Top hits
sort -t$'\t' -k6 -gr results/discriminator_panel_LATEST.tsv | head -20

# Specific candidate
grep "MY_PROTEIN_NAME" results/all_results.jsonl | jq .

# Non-XoxF Ln-class hits (Sharur-DB annotated)
awk -F'\t' '$7~/^Ln-/ && $8!~/xoxF/ {print $0}' results/ln_class_hits.tsv

# Top AFDB UniProt-stem hits
python3 -c "
import json
rows=[(r['ddE_kcal'],r['stem'],r['class'])
      for line in open('results/all_results.jsonl')
      for r in [json.loads(line)] if r['stem'].startswith('A') and r.get('ddE_kcal') is not None]
rows.sort(reverse=True)
for ddE,s,c in rows[:25]: print(f'{ddE:+7.2f}  {s:<55}  {c}')
"
```

---

## 8. Confirmed Ln-binder families (validated by discriminator + biology)

| family | example stem | typical ΔΔE | notes |
|---|---|---|---|
| **XoxF** (PQQ-MDH, La-evolved) | xoxf | +24.1 | Methodology gold standard |
| **Tannase** | tannasepara_NZ_JYMT, FSRD01000002_1_1028 | +20 to +27 | Two-domain Brady-rhizobium-style |
| **SilE** | SilE_La_sample_0, SilE_like_LanM_associated | +21 to +22 | Hand-curated; small Ln binder |
| **Thermolysin Ca₂ site** | A0A226QCS9 (Geobacillus) | +24.4 | Rediscovers Holmquist & Vallee 1974 |
| **Colin's PQQ-MDH positive controls** (00_Lav) | colinpqq_la_i0jwn7, a0acd6b9f2, c5atj3 | +11 to +23 | 11/11 confirm Ln-binding; reproducibility 0.2–3 kcal/mol vs ref runs |
| **Eukaryotic PQQ-ADHs** (NEW 2026-05-12) | a0a2v0prl1_domain3, a0a2v0nw76_domain1, a0a423vtm4_domain1 | +6 to +23 | 40+/100 score Ln-preferring or Ln-evolved |
| **AFDB Acidobacterium/Actinobacteriota** | A0A7V5CTA7 (porin), A0A226Q2Q9, A0A7V4XRG4 | +5 to +18 | Clade-level signal in spicy_lams scan |
| **small_unannotated_candidate_marco** | VAZP01000222_1_12, NZ_VITY01000007_1_28, … | +7 to +13 | 31-member novel family awaiting characterization |

---

## 9. Known false-positive families

| family | example stem | typical ΔΔE | why it scores positive |
|---|---|---|---|
| **M20_amidohydrolase** (Zn²⁺) | NZ_BPQI01000118_1_36 | +5 to +10 | XoxF-mimic 4-Asp/His cluster |
| **VIT1 Fe³⁺ uptake** | NZ_JACIDR010000005_1_9 | +30+ | 5-carboxylate cage; **status: hypothesized FP, possibly real Ln transporter** (Jacob 2026-05-09) |
| **Diiron rubrerythrin / bacterioferritin** | C5CIA3_bryobacter, C1F75 | ~+10 | 4-Glu binuclear pocket, La grabs the second Fe site |
| **CcmF heme biogenesis** | NC_010511_1_4382 | +6–8 | Cytochrome C assembly factor (Fe handling) |
| **Mn²⁺/Zn²⁺ small acidic pocket** | RXJC01000147_1_13, NZ_FNHS01000007_1_195 | < −30 (OUTLIER) | Pocket prefers smaller harder cation; SCF signature |

**Rule of thumb:** discriminator gives the chemistry; biology context (Pfam, KO, operon, genome neighborhood) gives the function. Cross-reference against your candidate's annotations before treating a positive ΔΔE as a real Ln-binder.

---

## 10. Inbox watcher details

```bash
# Manual single-run
bash scripts/process_inbox.sh
# Honors MAX_SUBMIT=N env var. Default 100. Flock-guarded; concurrent calls no-op.

# Daemon
bash scripts/inbox_watcher.sh > watcher.log 2>&1 &
# Pollings: every 45 min (sleep 2700)
# Conditions: inbox has CIFs AND queue total < 195 (5-slot safety margin)
# Stop: rm /tmp/inbox_watcher.lock
```

---

## 11. Troubleshooting

### Empty-carve signature (E[La] = −31.158775581143 bit-identical, n_atoms = 1)

This means the carve script found no eligible donors within 3.0 Å. Causes:
- TYR-OH or backbone-O coordinated site that wasn't being captured (fixed 2026-05-11 — verify your carve_generic.py has TYR in `SIDECHAIN_QM_ATOMS` and the backbone-O carving block).
- CYS/HIS/MET coordinated site — out of policy (Jacob's call). Classify as `CARVE_AMBIGUOUS`.
- Genuinely empty pocket (closest donor >4 Å). Classify as `SOLVENT_EXCLUDE`.

Fix: drop into `results/recarve_queue.tsv` and run `process_recarve_queue.py`. If the upstream fold daemon keeps emitting empty carves, its `carve_generic.py` path is stale — point it at the canonical script.

### "FAILED" jobs

- **First few seconds + grep showing SCF iters + no FINAL**: duplicate submission collision (legacy). Flock guard now prevents this.
- **`post_scan_fold` / `*_la_fold` / array tasks FAILED**: AlphaFold/Protenix folding from a different agent's pipeline. **Not the discriminator.** Filter with `grep -v fold`.
- **Protonation failure**: PDBFixer chokes on certain CIFs (oddly-numbered chains, unusual residue mappings). Workspace exists with no .xyz files. Check `slurm_*.err`.

### Outlier ΔΔE > 30 kcal/mol

Almost always a Mn²⁺/Zn²⁺ pocket the SCF can't accommodate La in. Sometimes a TightSCF + level-shift rerun helps:

```bash
# Edit sp_<stem>_La.inp:
#   ! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3
# → ! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 TightSCF
#   %scf shift shift 0.4 erroff 0.1 end end
# then rerun: $ORCA_PATH/orca sp_<stem>_La.inp > sp_<stem>_La_v2.out
```

### Severe asymmetric carve (\|ΔΔE\| > 100 kcal/mol)

Different ligand sets between Ca and La carves. Diff the XYZ:

```bash
head -1 <stem>_qm/<stem>_Ca_qm.xyz <stem>_qm/<stem>_La_qm.xyz
# If atom counts differ → known issue; re-carve with explicit --metal-resname
diff <(awk 'NR>2{print $1}' <stem>_qm/<stem>_Ca_qm.xyz | sort | uniq -c) \
     <(awk 'NR>2{print $1}' <stem>_qm/<stem>_La_qm.xyz | sort | uniq -c)
```

### Wrong ΔΔE sign (Ln-evolved hit you don't trust)

- Check carve composition: `head -1 *_qm.xyz` — atom counts La/Ca should match exactly (only the metal differs).
- Check charge sanity in `.inp`: La carve `* xyzfile -1 1` typically (3 carboxylates − 3 + La +3 = 0; wait — with PQQ −2 it's `-2 1`); Ca carve is `(La-charge − 1) 1`.
- Check element composition: `awk 'NR>2 {print $1}' x_qm.xyz | sort | uniq -c`. La/Ca carves should differ only in the metal.

### Disk filling up

```bash
bash scripts/cleanup_orca_scratch.sh
```

Frees ~30 MB per completed candidate. New jobs auto-clean on success; only legacy dirs need this.

### Pipeline reproducibility

- Cross-run agreement on identical structures: <0.05 kcal/mol on Marco family paralogs (validated).
- PDBFixer non-determinism causes 15–66 kcal/mol absolute energy variation between runs, but **it cancels in ΔΔE** because both metals see the same protonation. This is the vertical-swap design's whole point.
- Cross-pipeline reproducibility on Colin's 00_Lav set: new run vs older `colinpqq_la_*` runs match within **0.2–3 kcal/mol** on all 10 entries (validated 2026-05-12).

---

## 12. ORCA settings reference

```
! r²SCAN-3c NoAutostart CPCM(Water) DefGrid3
%maxcore 8000
%basis
  NewECP La "def2-ECP" end
  NewGTO La "def2-TZVP" end
end
* xyzfile {charge} {mult} {stem}_{kind}_qm.xyz
```

- **r²SCAN-3c**: composite functional with implicit basis correction. Standard for our panel.
- **NoAutostart**: forces fresh SCF, avoids accidentally reading a prior `.gbw` from a different system.
- **CPCM(Water)**: implicit solvent. Required — gas-phase ΔΔE on charged clusters is wildly overstabilized for the higher-charge state.
- **DefGrid3**: integration grid level 3. Tight enough for ΔΔE convergence; lower (1, 2) introduces ~1 kcal/mol noise.
- **La with def2-ECP + def2-TZVP**: large-core ECP keeps SCF stable. Don't switch to small-core without re-validating against Tier-1 controls.
- **Multiplicity**: La³⁺ is f⁰ → singlet (mult=1). Ca²⁺ is closed-shell → singlet. Both runs use mult=1.

Each SP ~5–15 min on a 64-thread memory-partition node; 4 SPs sequential is 20–60 min per candidate.

---

## 13. Provenance / further reading

- **Project state:** see [CONTEXT.md](CONTEXT.md) at the same level.
- **Vault state-of-pipeline (compaction-resilient bootstrap):** `~/jwestrob/obsidian-vault/agent-captures/2026-05-11_dft-discriminator-state-of-pipeline.md`
- **Project log:** `~/jwestrob/obsidian-vault/projects/lanthanide-binding/CONTEXT.md`
- **Methodology paper draft pointer:** `~/jwestrob/obsidian-vault/agent-captures/2026-05-07_ca-ln-dft-discriminator.md` and `2026-05-09_dft-discriminator-spicylams-scale.md`
- **Sharur DB annotations:** `/groups/banfield/users/jwestrob/bin/Sharur/data/spicy_lams/dft_la_candidates/{holo,holo_low_CN}/_classifications.tsv`
- **Colin's source CIFs:** `/groups/banfield/projects/multienv/corkscrew/supplementary_structures/structures_alone/{La-verified,Ca-verified,droideka,Euk}/`

---

## 14. Editing this doc

When you change pipeline behavior, update this file (and the changelog in `CONTEXT.md`). The agents/humans coming next will read it before they read the code.

**Last updated:** 2026-05-12 — PQQ-aware inbox routing, 3.0 Å + TYR + backbone-O carve dict, recarve queue + `process_recarve_queue.py`, severe-asymmetric-carve anomaly handling, scale to 1,453 candidates / 17 Ln-evolved + 277 Ln-preferring, Colin Euk PQQ-ADH discovery documented.
