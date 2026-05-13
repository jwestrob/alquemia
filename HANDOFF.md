# Handoff — Alchemical BVS validation on spici_lams LanMs

**Audience:** A fresh Claude session (or human collaborator) picking up this
task with no prior conversation context. Read this whole document before
running anything.

**One-line goal:** Validate Pass 2 of the augmented hard-cation scanner —
alchemical Ln³⁺ swap with local relaxation + Bond-Valence-Sum residual —
by reproducing **published LanM Ln-selectivity curves** on the user's
spici_lams verified-LanM set.

**Success criterion:** For at least one well-characterized LanM
(canonical *M. extorquens* LanM is the gold standard), the predicted
preferred Ln³⁺ — i.e., the cation in {La, Nd, Sm, Eu, Tb, Dy, Yb, Lu}
with the smallest post-relaxation BVS residual at the Ln-binding sites —
matches the published Kd-minimum **within ±1 position on the lanthanide
series**. Documented Kd minima for canonical LanM are at Sm/Eu/Nd
(mid-series). If we predict Eu, Sm, or Nd as the preferred Ln, that's a
pass. If we predict La or Lu, we have a methodology problem.

**Failure mode interpretation:** if Pass 2 reproduces the bell-curve
selectivity, it's validated as a quantitative Ln-selectivity predictor
and can be deployed on the novel candidates we surfaced earlier
(rifoxy, A0A9E3VGV4, fern peroxidases, A0A800K4K6 Gemmat calexcitin-fold).
If it fails, this run identifies *which* part of the pipeline (force
field choice, relaxation protocol, or BVS r₀ values) needs refinement
before we can trust it on novel data.

---

## Background — what's already done

Two prior vault notes establish the context:

- `~/jwestrob/obsidian-vault/agent-captures/2026-05-02_on-density-scanner-build.md`
  — scanner build, calibration on PDB anchors (LanM 8DQ2, calmodulin 1CLL,
  carbonic anhydrase 1CA2, zinc finger 1ZNF, etc.), TED scan, AFDB
  phylum scan (150k structures), rifoxy Chloroflexi deep-dive,
  A0A975IXS4 Verruco β-helix multi-Ca controls.

- `~/jwestrob/obsidian-vault/agent-captures/2026-05-02_ln-scanner-overnight-results.md`
  — overnight Phase 1 v2 / Phase 2 / fern fold pipeline, foldseek
  organization, **Pass 1 BVS scanner implementation + 82% validation
  on labeled positive controls**, reframing of all top candidates
  (mostly Ca²⁺_likely with permissive Ln substitution; rifoxy Site B
  reclassified as binuclear-Mn-likely; only A0A9E3VGV4 La_3 passes
  strict Ln-likely + CN ≥ 7 + iPTM ≥ 0.7 filter).

Read the **post-spec section** of the second note in particular — it
explains why Pass 2 (this experiment) is the methodologically critical
next step. Pass 1 implementation lives at:

- `on_scanner/bvs.py` — BVS computation, 10-cation panel, Brown 1985 r₀
- `on_scanner/augmented.py` — Pass 1 site-feature extraction +
  `classify_site()` rules-based filter
- `scripts/validate_pass1.py` — validation panel against PDB controls
- `scripts/rerank_pass1.py` — applies Pass 1 to all our existing folds

---

## Inputs

### spici_lams verified-LanM list

**Path:** `/groups/banfield/users/jwestrob/bin/Sharur/data/spicy_lams/exploration/expanded_lanm_sets.json`

Schema:
```json
{
  "verified":     [<list of 529 verified LanM IDs>],
  "expansion_ef5": [<list of 303 EF5-embedding-expanded LanM IDs>],
  "expanded_all": [<union, 832 IDs>]
}
```

LanM IDs are NCBI-style protein accessions (e.g., `NZ_AXAY01000027.1_101`,
`WP_244528404.1`, `JAKSYB010000063.1_159`). To get sequences, you'll
likely need to look these up in a DuckDB table the user maintains; the
relevant database is `/groups/banfield/users/jwestrob/bin/Sharur/data/spicy_lams/sharur.duckdb`
(the standard Sharur database for that dataset). Schema docs are in
`/groups/banfield/users/jwestrob/bin/Sharur/CLAUDE.md` and the global
config at `/home/jwestrob/.claude/CLAUDE.md`.

Quick-start: pick a small subset for the proof-of-concept run, NOT all
529. Recommended subset:

1. **Canonical M. extorquens LanM** — the original Cook & Cotruvo 2019
   target. Sequence likely under one of the verified IDs from
   *Methylobacterium extorquens* / *Methylorubrum extorquens* (search
   the IDs by joining with taxonomy in the Sharur DB).
2. 5-10 additional LanMs spanning genera (one *Bradyrhizobium*,
   one *Hansschlegelia*-clade if present, one each from a few other
   alphaproteobacterial methylotroph genera).

If the methodology works on those, scale up to the full 529 in a
follow-up SLURM batch.

### Published Ln-selectivity ground-truth data

**Primary reference:** Cook E.C., Cotruvo J.A. *et al.* "Lanmodulin: a
designed protein for selective lanthanide binding." Nature Chemistry,
2019. The supplementary information contains the Kd values across the
Ln series for canonical LanM. Bell curve peaks at Sm/Eu/Nd (~10⁻¹² M)
and falls off at La (~10⁻¹⁰ M) and Lu (~10⁻⁹ M).

**Secondary references (multiple LanM variants):**
- Mattocks J.A., Cotruvo J.A. et al. "Discovery and characterization of
  lanmodulin orthologs." JACS or Inorg. Chem. (~2020-2022)
- Hooff et al. on HansLanM (different selectivity profile, more
  selective for late Ln in some reports)

**Where to find values:** search Cotruvo group publications on PubMed
or Google Scholar; SI tables typically have Kd values per Ln per
LanM variant. If the user has these locally, the file will be under
`/groups/banfield/users/jwestrob/` somewhere — ask the user.

If you can't find published values for your subset, **fall back to**
the canonical LanM (well-documented bell-curve peak at Eu) as the
single validation point. One quantitative agreement is enough to claim
methodology works.

---

## Pipeline overview

### Phase A — fold each LanM with La³⁺

For each LanM in the chosen subset:

1. Extract protein sequence from spici_lams DB (or directly from a
   FASTA — if you can find one).
2. Determine number of EF-hand sites — typically 4 for canonical LanM,
   sometimes 3 in truncated paralogs. The Pro-locked-loop signature
   (regex `D.[DN].P[DN]G` or similar) helps distinguish active EF-hands
   from degraded ones. Conservative: use **N = 4** for most LanMs;
   override per-protein if known.
3. Run Protenix fold with `--ligands LA:N`. **Critical settings:**
   - Model: `protenix_base_20250630_v1.0.0` (the default `protenix-v2`
     is broken on this cluster — see global CLAUDE.md)
   - Partition: `gpu_h200` with `--gres=gpu:1`
   - **NO walltime** (`-t`) in any sbatch
   - Use foldit at `/groups/banfield/projects/environmental/sr/srvp2020/protenix/foldit.py`
   - Strip the `-t 4:00:00` from foldit's printed sbatch command before
     submitting (foldit prints with walltime; cluster doesn't want it)
4. Output CIF: `protenix_jobs/<lanm_name>/<lanm_name>/seed_101/predictions/<lanm_name>_sample_0.cif`
5. iPTM ≥ 0.85 expected for canonical LanMs with proper MSAs. If lower,
   recheck MSA quality (Phase 1 v1 had a no-MSA bug; see vault notes).

### Phase B — alchemical swap + local relaxation per Ln³⁺

For each LanM fold, for each metal site:

1. Identify each La³⁺ atom in the CIF.
2. For each target Ln³⁺ in {La, Nd, Sm, Eu, Tb, Dy, Yb, Lu}:
   a. Replace La³⁺ atom with target Ln³⁺ atom (just swap element
      identity in the CIF).
   b. Build OpenMM system with Merz 12-6-4 parameters for the target
      Ln³⁺. Reference parameters: **Li, P. & Merz, K.M. 2014, "Taking
      into account the ion-induced dipole interaction in the
      nonbonded model of ions" J. Chem. Theory Comput. 10(1):289-297**
      and Pengfei Li et al. updates (2017+).
      - Tabulated σ/ε/charge/scaling-factor values for Ln³⁺ are in the
        AmberTools params database under `frcmod.ions12_6_4_*` or
        Merz group's published tables.
      - For implicit solvent: GBSA (OpenMM `Implicit/OBC1` or `OBC2`).
   c. Local minimization:
      - **Backbone heavy-atom restraint**: position-restrain all CA, N,
        C, O atoms with a stiff harmonic spring (e.g., 50 kcal/mol/Å²).
      - **Free**: sidechain heavy atoms within ~6 Å of the target metal
        atom. Restrain everything else.
      - Convergence: gradient norm < 0.01 kcal/mol/Å, or 5000 iters
        whichever first. L-BFGS in OpenMM is fine.
   d. After minimization, compute BVS for target Ln³⁺ at the new metal
      position using `on_scanner/bvs.py:compute_bvs()`. The first-shell
      donor list is the same set of donors as the un-relaxed structure
      (since backbone is fixed); their distances will have shifted.
   e. Record:
      - `bvs_relaxed_<Ln>`
      - `bvs_residual_<Ln>` = |BVS - 3| (target valence)
      - `relaxation_rmsd_<Ln>` (heavy-atom RMSD of the released sidechains)
      - Per-donor distance change Δd

### Phase C — selectivity prediction + comparison

For each (LanM, site) pair:

1. Sort the Ln panel by `bvs_residual_<Ln>` ascending.
2. The cation with smallest residual is the predicted preferred Ln.
3. Plot the residual curve across the Ln series — should look like
   the inverse of the published Kd curve (low residual ≈ high affinity).
4. Compare to published Kd-minimum.

**Output table** (TSV, one row per LanM-site pair):
```
lanm_id  site_id  predicted_pref_Ln  res_La  res_Nd  res_Sm  res_Eu  res_Tb  res_Dy  res_Yb  res_Lu  rmsd_La  rmsd_Nd  ...  rmsd_Lu  iptm
```

---

## Force field & parameter sources

- **Ln³⁺ Merz 12-6-4 parameters**: Li & Merz 2014/2017. Available in
  AmberTools (`$AMBERHOME/dat/leap/parm/frcmod.ions12_6_4_*`) and via
  the Merz group's online repository.
- **Ca²⁺ for control runs** (recommended sanity check): Liu & Thiel 2007
  Ca²⁺-specific r₀ for proteins, or the Joung-Cheatham parameters if
  going pure-MM. Don't mix parameter families within a single MD setup.
- **Implicit solvent**: GBSA OBC2 in OpenMM. Cheap, sufficient for this
  question. Don't bother with explicit water for the proof-of-concept;
  it's 10× more expensive and unlikely to change conclusions.
- **Force field for protein**: Amber ff14SB or ff19SB. Either is fine.
  Make sure your Ln³⁺ parameters were derived in the matching
  framework — Merz 12-6-4 was developed in Amber.

**Software:**

- **Protenix** for folding (existing pipeline; see `protenix/foldit.py`)
- **OpenMM** for relaxation (Python interface, easier than GROMACS for
  scripting). `pip install openmm` in the user's `fep` conda env may
  already have it; check `/home/jwestrob/miniconda3/envs/fep/bin/python3
  -c "import openmm"` first.
- **gemmi** for CIF parsing (already in `fep` env).
- **AmberTools** for parameter loading (`frcmod` files). Likely
  installed; confirm with `which tleap`.

If any of these are missing, set them up via SLURM job in a fresh conda
env; do NOT pip install on the login node for anything that touches
heavy compute libs.

---

## Existing infrastructure to reuse

```
/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/
├── on_scanner/
│   ├── bvs.py                        # BVS computation, cation panel, call_cation()
│   ├── augmented.py                  # Pass 1 site features + classify_site()
│   ├── geometry.py                   # donor enumeration, KD-trees
│   ├── scoring.py                    # composite geometric score (legacy)
│   └── scanner.py                    # grid scan + NMS (legacy)
├── scripts/
│   ├── validate_pass1.py             # validation panel runner
│   ├── rerank_pass1.py               # apply Pass 1 to existing folds
│   ├── batch_scan.py                 # parallel batch scanner
│   └── run_phylum_scan.sh            # SLURM wrapper for AFDB phyla
├── results/overnight_v2/             # Phase 1 v2, Phase 2, fern fold + foldseek outputs
├── calibration/pdb_controls/         # 8DQ2 (LanM), 1CLL (calmodulin), etc.
└── alchemical_bvs/                   # YOU ARE HERE — write your scripts under this dir
```

**Reuse `bvs.py:compute_bvs()` directly** for post-relaxation BVS — it's
already validated and handles the cation panel correctly.

---

## Hard constraints

Read the global CLAUDE.md (`/home/jwestrob/.claude/CLAUDE.md`) carefully
before running anything. Most importantly:

1. **No login-node compute.** Anything touching MM minimization, MSA
   search, or batch parsing of >100 files goes through SLURM.
2. **Full nodes always.** Use `--exclusive` and `$SLURM_CPUS_ON_NODE`
   (don't hardcode `--cpus-per-task=N` to a literal number).
3. **No walltime.** Strip `-t` from any sbatch command. Cluster
   uses `#SBATCH -p standard` or `#SBATCH -p gpu_h200`, no `-t`.
4. **Use the right Protenix model**: `protenix_base_20250630_v1.0.0`.
   The `protenix-v2` default is broken on this cluster.
5. **Set up SLURM completion watchers** with
   `bash -c 'until ! squeue -j JOBID -h ... grep -q . ; do sleep 60; done; ...'`
   and `run_in_background=true` so completion fires a notification.
6. Use `/home/jwestrob/miniconda3/envs/fep/bin/python3` for scripts that
   need gemmi.
7. `/home/jwestrob` is on a 100 GB SSD — write outputs to
   `/groups/banfield/...`, never to `/home/jwestrob/`.

---

## Pre-checks before running

Before doing ANYTHING expensive:

1. **Verify the spici_lams data is reachable**:
   ```bash
   ls /groups/banfield/users/jwestrob/bin/Sharur/data/spicy_lams/exploration/expanded_lanm_sets.json
   /home/jwestrob/miniconda3/envs/fep/bin/python3 -c "import json; print(len(json.load(open('/groups/banfield/users/jwestrob/bin/Sharur/data/spicy_lams/exploration/expanded_lanm_sets.json'))['verified']))"
   ```
   Expect: 529.

2. **Check Sharur DuckDB** for sequence retrieval:
   ```bash
   ls /groups/banfield/users/jwestrob/bin/Sharur/data/spicy_lams/sharur.duckdb
   ```

3. **Verify OpenMM and Merz 12-6-4 parameter availability**:
   ```bash
   /home/jwestrob/miniconda3/envs/fep/bin/python3 -c "import openmm; print(openmm.__version__)"
   ls $AMBERHOME/dat/leap/parm/frcmod.ions12_6_4* 2>/dev/null || echo "AMBERHOME not set; check tleap"
   which tleap
   ```

4. **Pick canonical LanM** for the first run. The verified set has 529
   IDs. Find one from *Methylorubrum extorquens* or
   *Methylobacterium extorquens* by joining with the taxonomy info in
   `sharur.duckdb`. If you can't find one, ask the user — they have
   the curated list in their head.

5. **Don't fold all 529 at once.** Start with 1 (canonical LanM).
   If iPTM ≥ 0.85 and EF-hand sites all populated, expand to 5-10.
   Only after Phase B+C succeeds on those, scale to 529.

---

## Output schema

Final output goes to:

- `alchemical_bvs/results/lanm_selectivity_predictions.tsv` — per-site
  BVS residual across Ln panel + predicted preferred Ln
- `alchemical_bvs/results/comparison_to_published.tsv` — predicted vs
  published preferred Ln per LanM, with discrepancy notes
- `alchemical_bvs/results/per_site_relaxation_logs/<lanm>_<site>_<Ln>.log`
  — minimization energy trace per swap
- `alchemical_bvs/REPORT.md` — top-level summary: did we reproduce the
  bell curve? For which LanMs? What went wrong on failures?

If methodology validates, write a vault note at
`~/jwestrob/obsidian-vault/agent-captures/YYYY-MM-DD_alchemical-bvs-validation.md`
and append a "Pass 2 validated, deployable on novel candidates" section
to the existing `2026-05-02_ln-scanner-overnight-results.md`.

---

## Open questions you may face

1. **Where exactly are the LanM sequences?** Verified IDs are in the
   JSON; sequences need to be retrieved. Likely path: query the Sharur
   DuckDB at the proteins table by accession. If unclear, ask user.

2. **How many EF-hand sites per LanM?** Default to N=4. Use the regex
   `D.[DN].P[DN]G` or similar to count Pro-locked loops as a sanity
   check. If a LanM has only 3 hits, fold with `LA:3`.

3. **Force field choice for non-canonical Ln³⁺ in OpenMM**: Merz 12-6-4
   has parameters for La, Nd, Eu, Gd, Tb, Dy, Ho, Er, Yb, Lu. If a
   specific Ln in the panel doesn't have published 12-6-4 params,
   either skip that Ln or use the closest available (e.g., Sm ≈ Eu).

4. **Validation when published Kd is missing**: if the user's selected
   LanM doesn't have published Ln-selectivity data, fall back to the
   canonical LanM (Cook & Cotruvo 2019) as the single validation point.
   One quantitative agreement is enough.

5. **Relaxation convergence failure**: if local min doesn't converge,
   it's usually a clash from the alchemical swap. Try smaller initial
   step size, longer minimization, or switch from L-BFGS to steepest-
   descent for the first 100 steps. Don't blindly raise iteration count.

6. **What if Pass 2 says LanM prefers La or Lu (the wrong answer)?**
   Possible causes:
   - Wrong Merz 12-6-4 parameters loaded.
   - Backbone is too rigid; some flexibility needed in EF-hand loop.
   - Sidechain set too small (try 8 Å instead of 6 Å).
   - LanM-specific Pro-locked-loop torsions need explicit dihedral
     handling (probably not, but check).
   Document the failure mode in REPORT.md and stop. Do NOT chain
   workarounds; the methodology has a real problem to debug.

---

## Cross-references

- Vault notes (must read for context):
  - `~/jwestrob/obsidian-vault/agent-captures/2026-05-02_on-density-scanner-build.md`
  - `~/jwestrob/obsidian-vault/agent-captures/2026-05-02_ln-scanner-overnight-results.md`
- Web-Opus spec (in user's prior conversation, not in vault): the spec
  for "Augmented Hard-Cation Site Scanner." Pass 2 design comes from
  there; relaxation protocol details should follow the spec.
- Sharur project CLAUDE.md (project conventions, key databases, ELSA):
  `/groups/banfield/users/jwestrob/bin/Sharur/CLAUDE.md`
- Global CLAUDE.md (cluster discipline, SLURM, pyenv vs conda):
  `/home/jwestrob/.claude/CLAUDE.md`

---

## Definition of done

You're done when:
1. Canonical LanM (Cook & Cotruvo 2019 reference) has been folded with
   La and alchemically swapped through the Ln panel; predicted preferred
   Ln matches published Kd-minimum within ±1 position on the series.
2. `alchemical_bvs/REPORT.md` documents the result with sufficient
   detail that another agent could replicate.
3. If validation succeeds, the scanner project's vault note is updated
   to reflect "Pass 2 methodology validated; ready for novel-candidate
   deployment." If validation fails, the vault note documents the
   failure mode and what would need to change.

Stop after the canonical LanM validation. Do NOT scale to all 529 LanMs
in this session — that's a follow-up scope.
