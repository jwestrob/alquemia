# All curated PLM XoxF: frozen-core preparation adapter

Approved scope: the 176 distinct full PLM proteins underlying the user's curated
XoxF tree selection. The two completed AF3 fixed-core results are reused; up to
174 new selected models can produce 348 La/Ca endpoints. This adapter performs
selection and chemical preparation only. The PLM parent runner owns submissions,
ORCA execution, monitoring and final publication.

Authorization: `EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/xoxf_all/authorization.json`.
Frozen inventory SHA256:
`13d85f07c00e759ee0997d622da538d4569d282971771a2aa5c4282113002884`.

## Interface

Use the pinned `lanm_qmmm` Python environment in an allocated CPU Slurm job:

```text
prepare_batch.py --inventory ABS --inventory-sha256 SHA \
  --fold-input-manifest ABS --fold-input-sha256 SHA \
  --fold-results ABS --approval ABS --approval-sha256 SHA \
  --output ABS [--workers auto]
```

The parent runs the folding collector first, including after an upstream stage
failure. The runtime fold-results manifest must pin the frozen inventory and
fold-input manifest. The output directory must not exist. The agreed production
output is `workspaces/plm_xoxf_all_fixed_core_20260916/prepared_batch`.

Outputs at the output root:

- `prepared_pairs.json`: frozen-executor-compatible
  `plm.adh9.prepared_pairs.v1`; **new ready cases only** in `cases`.
- `target_outcomes.json`: `per_target` contains all 176 proteins once, with
  `ready`, `reused`, `unsupported` or `preparation_failed` status.
- `candidate_manifest.json`, candidate and final implementation pins,
  `geometry_validation.json` and `resources.json`.
- `selection/<protein>/selection.json`: all three model outcomes, selected
  metrics, mapped residues and candidate Arg/Lys partners.
- `worker_results/<protein>.json`: one durable result per preparation attempt;
  raw preparation, minimizer observations, and any recovery remain separate.
- `worker_dispatch/<protein>.dispatch.json`: isolated process exit, stdout/stderr
  and retained per-target outcome, including native crashes. No retries/timeouts.

`reused_results` and `unsupported_targets` never enter the quantum task list.
An empty ready list remains a valid completed preparation inventory; the parent
skips quantum execution and publishes the available reused results. Geometry
aggregate PASS refers only to admitted cases, explicitly retaining unsupported
counts and target outcomes.

## Selection and chemistry

For each new protein, inspect the three existing AF3 models independently. Rank
structurally valid models by typed CN descending, protein–La iPTM descending,
protein CN descending, then sample index ascending. Apply CN >=7 at 3.1 A,
protein–La iPTM >=0.9 and direct N <=2 **after selection**. No alternative sample
replaces a failed selected model.

Four positions come from the independently validated saved alignment mappings:
conserved E/N, catalytic WD and its D+2 homolog. Actual residue identities are
required. A unique Arg/Lys side-chain contact to the catalytic Asp within 3.5 A
defines the cationic partner; residue offsets, partner identity, atom count and
formula are not hardcoded. Unsupported or ambiguous cores remain unscored.

The original frozen `prepare_candidates.py` calls the original fixed-core carver
and standard-only protonator. Chemistry remains dry oxidized PQQ3− plus the
released homologous side chains, without a radius-expanded carve, heavy-atom
relaxation, new ligands or field embedding. D+2 identity controls whether its
side chain is included and thus the charge. The separate parent uses the same
native ORCA r2SCAN-3c/CPCM(Water)/DefGrid3 endpoints and released decision bands.
No MACE computations or calibration refitting are part of this work.

## Hydrogen preparation and limitations

`hydrogen_preparation.py` copies the observed-original preparation and conditional
recovery functions from the completed two-XoxF adapter, with new workspace and
approval injection. The original helpers remain unchanged and independently
pinned. Each target runs in a fresh process, with OpenMM/BLAS/OpenMP threads set
to one; up to 80% of available allocation memory at a conservative 8 GiB per
worker limits parallelism. A lightweight thread dispatcher launches one explicit
subprocess per target; native child crashes become recorded per-target failures
without losing results or blocking sibling targets. This replaces the initial
process-pool controller only; preparation and selection are unchanged.

The original 50-step minimization is observed without changing its objective,
seed or iteration ceiling. Its residual H forces are recorded. **Clean core
geometry does not establish force convergence**; the previously completed two
XoxFs had clean geometry but RMS residual forces around 40–42 kJ/mol/nm.

Only actual malformed H/core geometry triggers the authorized numerical recovery:
restart retained H coordinates with the exact captured OpenMM System, same CPU
platform/objective/tolerance, at most 1,000 steps. All heavy coordinates, protonation
identities, PQQ atoms/H, caps and core ordering must remain fixed. Recovery must
reach full-precision H force RMS <=1; PDB serialization forces are separately
measured. Any preparation deviation and lack of a calibration rerun are explicit.

Scores classify protocol-specific calculated contrasts, not measured affinity,
substrate specificity or physiological metal use. The selected proteins are
phylogenetically curated XoxFs, not experimentally characterized validation controls.

## Checks completed before launch

Fourteen lightweight tests pass. They inspect the six existing
AF3 model metrics and two exact role mappings, test selected-model rejection
without replacement, reject the retained original malformed ADH9 H geometry and
changed frozen atoms, and test a synthetic zero-ready batch with complete
176-target accounting, plus native-crash and recorded-worker-failure receipts
without retry. **No new production protonation, folding or ORCA was run.**

The parent `run_quantum.py` independently accepted the two existing immutable
prepared cores using a temporary compatible outcome fixture. A first direct call
against the historical two-target directory lacked the newly required
`target_outcomes.json`; the temporary fixture supplies only that new interface
metadata, without editing historical files or running endpoints.

Validation receipts are in `software_validation.json` and `executor_review.json`.
