# DFT Ca²⁺/Ln³⁺ Discriminator — Colin's Handoff Package

## Context (read this first)

This is a self-contained pipeline for scoring protein structures by their
intrinsic Ca²⁺ vs Ln³⁺ binding preference using DFT cluster calculations.

Jacob built and validated this on a 50+ protein panel. He's hitting per-user
SLURM job limits, so we're distributing compute by having Colin's account
submit some structures.

**The science:** for each candidate, we carve a QM cluster (metal + first-shell
residues + cofactor if present), run 4 single-point DFT calculations
(La form, Ca form, apo form, bulk water reference) at r²SCAN-3c CPCM(Water),
then compute ΔΔE(Ca − La) using a fixed aquo reference.

- **ΔΔE > +20 kcal/mol** → Ln-evolved (architecturally amplified Ln preference; LanM, XoxF tier)
- **+5 to +20** → Ln-preferring (intrinsic carboxylate-cluster preference)
- **0 to +5** → marginal
- **−5 to 0** → ambiguous
- **< −5** → Ca-evolved (active-site geometry tuned for Ca²⁺; MxaF tier)

The pipeline is deterministic at the ΔΔE level (we measured ~0.5 kcal/mol
replicate noise even when absolute energies vary by 15-25 kcal/mol due to
PDBFixer protonation choices).

## What you (Colin's Claude) need to do

You are running on the same SLURM cluster as Jacob (banfield group).
ORCA 6.1.1 lives at `/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/`.
The shared Python env is `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/`.
Both are accessible from your account.

**Do not modify the scripts** — they were validated against gold-standard reference
proteins (XoxF +24.1, MxaF −6.5, LanM +47.5). Changing functional, basis set,
charge bookkeeping, or aquo references will break comparability.

## Step-by-step

### Step 1: Pick a workspace directory

Create a fresh directory in your account's writable area for the run output:

```bash
mkdir -p /groups/banfield/users/colinr/discriminator_workspace
cd /groups/banfield/users/colinr/discriminator_workspace
```

(or wherever you have write access)

### Step 2: Pick what to run

The CIF source directory is:
```
/groups/banfield/projects/multienv/corkscrew/supplementary_structures/structures_alone/
├── La-verified/    (11 PQQ-MDH structures, La-bound — controls)
├── Ca-verified/    (14 PQQ-MDH structures, originally Ca but modeled with La — controls)
├── Euk/            (167 eukaryotic PQQ-domain structures)
└── droideka/       (1 Acidobacterium with 6 La and NO PQQ — special case)
```

Jacob is running La-verified, Ca-verified, and droideka. **You should run the
167 Euk structures** (or a subset).

### Step 3: Submit a batch

For PQQ-bound structures (everything except droideka):

```bash
HANDOFF=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/colin_handoff

# Run the full Euk batch (will submit ~167 SLURM jobs over time as your
# QOS allows — the script is idempotent so re-running picks up where it left off):
bash $HANDOFF/run_pqq_batch.sh \
  /groups/banfield/projects/multienv/corkscrew/supplementary_structures/structures_alone/Euk \
  /groups/banfield/users/colinr/discriminator_workspace/euk_run
```

This will:
- For each `*.cif`, normalize → protonate → carve → submit
- Skip any that already have a converged La SP (idempotent)
- Print a summary

### Step 4: Monitor

Each candidate runs 4 ORCA single-points sequentially, ~30-60 min wallclock
per SP. Total ~2-4 hr per candidate, but they run in parallel up to your QOS
limit (typically 8-10 concurrent jobs).

```bash
# Live job queue:
squeue -u $USER -o "%.10i %.40j %.10T %.10M"

# Per-output check:
ls /groups/banfield/users/colinr/discriminator_workspace/euk_run/*_qm/sp_*_La.out 2>/dev/null | wc -l
echo "candidates with La SP done"
```

### Step 5: Aggregate results

When all (or most) jobs are done:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  $HANDOFF/scripts/analyze_results.py \
  /groups/banfield/users/colinr/discriminator_workspace/euk_run
```

This writes `discriminator_results.tsv` in the workspace with one row per
candidate: ΔΔE in kcal/mol, atom count, classification.

### Step 6: Hand back to Jacob

```bash
# The TSV is what Jacob needs:
cp /groups/banfield/users/colinr/discriminator_workspace/euk_run/discriminator_results.tsv \
   /tmp/colin_euk_results.tsv

# Tell Jacob the path (in Slack or whatever)
```

## Special case: droideka

The Acidobacterium "droideka" structure has 6 La³⁺ bound and NO PQQ. Each
La site needs to be carved separately. **Use the dedicated script:**

```bash
DROIDEKA_CIF=/groups/banfield/projects/multienv/corkscrew/supplementary_structures/structures_alone/droideka/acidobacteria_hcr00247.1-la6_model.cif

bash $HANDOFF/run_droideka_sites.sh "$DROIDEKA_CIF" \
  /groups/banfield/users/colinr/discriminator_workspace/droideka_run
```

Outputs 6 subdirs (one per La site at chains B/C/D/E/F/G). Aggregate same way.

(Jacob is already running droideka — but if it fails his queue, this is the
fallback.)

## How to interpret your results

Once `analyze_results.py` finishes, you'll see something like:

```
stem                                          atoms   ddE_kcal  class
A0A0D2A1I3_domain1chaina_pqq_la_model           58       +18.5  Ln-preferring
A0A0D2MSM5_domain1_contigchaina_pqq_la_model    44        −7.2  Ca-evolved
...
```

**Interesting hits to flag for Jacob:**
- Any with **|ΔΔE| > 30 kcal/mol** → likely a geometry artifact (short M-O distances,
  Protenix prediction errors) — flag as outlier; do NOT include in the final panel
- **Ln-evolved (≥+20)** → high-priority novel candidates worth wet-lab follow-up
- **Ln-preferring (+5 to +20)** → solid carboxylate-cluster Ln preference,
  comparable to canonical EF-hand proteins
- **Ca-evolved (<−5)** → architecture tuned for Ca²⁺, e.g. MDH catalytic site

**Sanity check:** the La-verified controls should land Ln-class (positive ΔΔE),
the Ca-verified controls should land Ca-class (negative ΔΔE). Use these as
your validation that the pipeline is working.

## What to do if something fails

### A candidate's La SP doesn't converge

Look at `*_qm/sp_*_La.out` — if it shows TRAH or SOSCF iterations, it's
just slow (sometimes 1-2 hr). If it shows "ERROR" or "NOT CONVERGED",
restart with TightSCF. Email Jacob the protein ID.

### sbatch returns "QOSMaxJobsPerUserLimit"

Normal — your account caps concurrent jobs. The pending jobs will start
when others finish. Just wait.

### Carve produces 0 carboxylates / weird charge

Open `*_qm/<stem>_protonated.pdb` and check that the LA HETATM is present
and the surrounding residues look like a proper coordination site. If the
AF3 prediction is bad (e.g., La placed in a non-carboxylate pocket), the
carve will be junk — flag and skip.

### "FINAL SINGLE POINT ENERGY" not in output

Job died mid-SCF. Check `slurm_*.out` for the actual error. Most common
cause: ORCA hit a node memory issue or a duplicate keyword in the input.
Email Jacob.

## File inventory

```
colin_handoff/
├── README.md                    (this file)
├── run_pqq_batch.sh             (main entry point for PQQ-MDH structures)
├── run_droideka_sites.sh        (special 6-site handler for droideka)
└── scripts/
    ├── normalize_af3_cif.py     (renames AF3 LIG_* residues → standard names)
    ├── protonate_cif.py         (PDBFixer adds Hs at pH 7)
    ├── carve_with_pqq.py        (carves QM cluster including PQQ cofactor)
    ├── carve_generic.py         (carves QM cluster, protein-only, with site selector)
    └── analyze_results.py       (computes ΔΔE TSV from completed SPs)
```

## The aquo reference (DO NOT change)

These are r²SCAN-3c CPCM(Water) DefGrid3 single-point energies on idealized
[M(H₂O)₈]^q⁺ structures, hardcoded in `analyze_results.py`:

```
E_aquo_Ca²⁺ (CN8) = -1288.921041 Ha
E_aquo_La³⁺ (CN8) =  -642.856486 Ha
DIFF_AQUO         =   -646.064556 Ha
```

`ΔΔE(Ca - La) = (E_holo_Ca - E_holo_La) - DIFF_AQUO` in Hartree, then
multiply by 627.5095 for kcal/mol.

## Reach out if stuck

Jacob is on Slack. Send him:
1. The protein ID(s) and ΔΔE values
2. The path to your workspace
3. Any error messages from slurm_*.out files

— Jacob's Claude (built this handoff 2026-05-08)
