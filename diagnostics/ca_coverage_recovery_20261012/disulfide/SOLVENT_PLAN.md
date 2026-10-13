# Corrected solvent selection for the one repaired MMOL_2048 source

Run only after native2875 terminates with an actual one-source collection. Use
its exact admitted original/prepared coordinates and candidate matrix. New
source preparation is proven by2870, including the exact old group chemical-state
key and two unchanged real controls. The original125-source row remains an
unavailable primary result; this is a separately named preparation recovery.

Use the same frozen runtime and scientific settings as the completed corrected
Ca-conditioned solvent collection: published-parameter standalone xTB and the
independent tblite check, unchanged component-v2 acceptance, original energy
algebra, candidate/tie policy and exact canonical-only bands. Keep historical
absolute-v1 qualification failures visible. Do not introduce native ORCA's old
La shell assignment or alter the MACE model.

The adapter reuses MxaF SOLVENT.stage() for actual12CPU allocation validation,
then the unchanged GFN2_STAGE.run, GFN2.py dry_run/collect/verified_cell and
their pinned raw executors. This source had no valid original envelope pool;
its six candidate cells require at most12 distinct solvent inputs per backend
(two media each),24 fresh backend calls total. Exact input identities are
deduplicated. No new energy call for either real preparation control. No new
optimization, DFT, folds, fit or production modification.

Original core/source identities, updated protonation/core-parent provenance,
prepared maps, state anchor, paired charge/inventory and actual MACE receipts
are checked before solvent work. A missing required cell leaves the accommodated
pool unavailable; any valid origin is explicitly separate. The existing executor
retains partial failures and forbids repeating started unfinished input points.

Write one-source result and a125-row overlay keyed by source_case_id (the old
excluded row has case_id=null). Preserve the entire original row. If the actual
MxaF recovery overlay exists, verify its primary and each actual recovery row
before including its two distinctly named technical recoveries. Otherwise mark
that optional integration unavailable; do not delay or invent it. Counts are
recovery annotations, not replacements of the original primary experiment.

```bash
sbatch --parsable --dependency=afterok:2875 diagnostics/ca_coverage_recovery_20261012/disulfide/solvent.sbatch "$PWD/workspaces/ca_coverage_recovery_20261012/disulfide/SOLVENT_CONFIG_v1.json"
```

Root owns submission and completion wake.12CPU shared Slurm allocation, existing
one-thread child affinity, no memory/time/account/QOS flags. Metadata preparation
and syntax checks only occurred on login; scientific qualification/inputs remain
unrun until the actual Slurm stage. Report actual receipts separately from all
earlier work and reused preparation; do not claim a new-site speed benchmark.
