# Contributing

This document is for collaborators who want to submit candidates from a
different cluster account, extend the carve dictionary, or otherwise modify
the methodology. For day-to-day operation see [HOWTO.md](HOWTO.md).

If you are Colin (the original cross-account use case), the older self-contained
handoff package at [`colin_handoff/README.md`](colin_handoff/README.md) is
still current.

---

## Cluster prerequisites

The pipeline expects to run on the shared SLURM cluster used by the
Banfield group at LBNL / UC Berkeley.

**Group membership.** You need to be in the `banfield` Unix group so that
the shared workspace (`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/`) and the ORCA install (`/home/jwestrob/jwestrob/bin/ORCA/...`)
are readable from your account.

**ORCA path** (hard-coded in `scripts/carve_generic.py` and `carve_with_pqq.py`):
```
/home/jwestrob/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
```

**Conda envs** (both pre-provisioned, both readable group-wide):

- `lanm_qmmm` at `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/` —
  used by every script except `protonate_cif.py`. Spec: [environment.yml](environment.yml).
- `fep` at `/home/jwestrob/miniconda3/envs/fep` — PDBFixer + OpenMM, used only
  by `protonate_cif.py`. Spec: `requirements_fep.txt` placeholder note.

**SLURM partition.** `memory` partition, 64-core nodes, full-node `--exclusive`,
OpenMP threading (no MPI). See [PIPELINE.md §5](PIPELINE.md).

---

## Submitting from a different account

### Option A — symlinks into the shared inbox (preferred for one-offs)

If your work has read access to the shared `alchemical_bvs/inbox/` directory
(it should — `inbox/` is group-writable), drop CIFs as **symlinks** rather
than copies:

```bash
ALCH=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs

# Symlink your CIFs (cheap, no disk duplication, easy to revoke)
ln -s /path/to/your/candidate.cif "$ALCH/inbox/your_candidate.cif"

# The 45-min inbox watcher will pick it up. Or trigger immediately:
bash "$ALCH/scripts/process_inbox.sh"
```

The SLURM job runs under whichever user triggered `process_inbox.sh`; if you
ran it yourself, the job is in your queue (and counts against your QOS).
If you let the watcher pick it up, the job is in the watcher owner's queue.

Use this for ≤ a few dozen CIFs you want scored against the master panel.

### Option B — separate workspace (preferred for bulk batches)

If you have hundreds or thousands of CIFs and don't want to flood the shared
inbox, set up a parallel workspace using the handoff scripts:

```bash
HANDOFF=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/colin_handoff

mkdir -p /groups/banfield/users/$USER/discriminator_workspace
bash $HANDOFF/run_pqq_batch.sh \
  /path/to/source/dir_with_cifs \
  /groups/banfield/users/$USER/discriminator_workspace/my_run
```

This:

- creates a parallel workspace tree
- runs `normalize_af3_cif → protonate → carve_with_pqq → sbatch` per CIF
- is idempotent: re-running skips stems that already have a converged La SP
- aggregates results into `discriminator_results.tsv` in your workspace

When the batch is done, hand the resulting TSV back via Slack or copy it
into `results/` for inclusion in the master panel (with Jacob's review).
See [colin_handoff/README.md](colin_handoff/README.md) for the full procedure.

---

## Extending the carve dictionary

The carve dict in `scripts/carve_generic.py` is currently restricted to:

```
ASP, GLU, ASN, GLN, SER, TYR    (+ backbone-O for any residue)
```

This is a deliberate policy choice; see [METHODS.md §3 ("Policy: what is not in the dict")](METHODS.md#policy-what-is-not-in-the-dict).
**Adding a new residue is a methodology change**, not a bugfix. If you have
a use case, the required workflow is:

1. **Open a discussion** with Jacob first. Adding S-donors (CYS, MET) or
   N-aromatic donors (HIS) changes what the discriminator measures (the
   chemistry shifts from hard-O / hard-N to mixed-donor). The class
   thresholds in [TIERS.md](TIERS.md) were calibrated against the current
   dict; they may need to move.
2. **Edit `SIDECHAIN_QM_ATOMS`** in `scripts/carve_generic.py`. Include all
   sidechain heavy atoms from Cβ to the donor, plus the link-H cap atoms.
   Mirror the existing dict entries for style.
3. **If the new residue can coordinate via a backbone atom**, no change
   needed — backbone-O carving already covers that.
4. **Re-run the Tier-0 (negative-control) and Tier-1 (positive-control) panel**:
   ```bash
   # Tier-1 positive controls: XoxF, MxaF, LanM EF1, tannase, Colin's 11 La-verified
   # Tier-0 negative controls: calmodulin (1CLL), parvalbumin (1B9A), calbindin (1A75), fern peroxidase
   ```
   Confirm the gold-standard ΔΔE values (XoxF +24, MxaF −6.5, etc.) are
   preserved within ~1 kcal/mol. If they shift by more, the threshold
   calibration needs revisiting.
5. **Add a CHANGELOG entry** documenting the change, the rationale, and the
   re-validation results.
6. **Update [METHODS.md §3](METHODS.md)** and any in-tree references that
   list the dict.

The same workflow applies if you want to widen the geometric cutoff
(`FIRST_SHELL_CUT`), add ring waters past the first shell, or change the
link-H placement convention.

---

## Modifying ORCA settings

`r²SCAN-3c / CPCM(Water) / DefGrid3 / NoAutostart` is locked at the
discriminator level. Changing the functional or basis set requires:

1. Re-running the aquo references (`[Ca(H₂O)₆]²⁺` and `[La(H₂O)₆]³⁺` at the
   new level of theory) to get a new `DIFF_AQUO`.
2. Updating `E_AQUO_CA`, `E_AQUO_LA`, and `DIFF_AQUO` constants in
   `scripts/update_results_jsonl.py`.
3. Re-running the Tier-1 validation panel to confirm class assignments are
   preserved (see [VALIDATION.md §7](VALIDATION.md) for the B97-3c
   cross-check that has already been done).
4. CHANGELOG entry.

---

## Repository hygiene

- **Do not commit:** `*_qm/` workspace directories (they are large and
  reproducible from the inbox CIFs), `slurm_*.{out,err}` logs,
  `*.gbw / *.tmp / *.densities` ORCA scratch.
- **Do commit:** changes to `scripts/`, `*.md` docs, `results/*.tsv` /
  `results/*.jsonl` (these are the persistent panel state — append-only
  changes only).
- **Authoritative state files** (`CONTEXT.md`, `HOWTO.md`, and this doc) should
  be updated when behaviour changes, not after the fact.
- **Pre-flight a methodology change** by running it against a known control
  (e.g. `xoxf_qm`) before applying it to the whole panel.

---

## Reaching out

- **Methodology questions / new candidate batches:** Jacob West-Roberts (Banfield Lab).
- **Cross-account submission troubles:** the handoff package at
  [colin_handoff/](colin_handoff/) has working examples.
- **Bugs / surprises in the pipeline:** open an issue with the failing stem,
  the workspace path, and the `slurm_*.err` if any.
