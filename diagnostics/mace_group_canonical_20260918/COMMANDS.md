# Typed-group PQQ operations

Run from the repository root. These are opt-in research operations; the baseline
is unchanged. Read PLAN.md and SOURCE_CLARIFICATION.md first. Existing output
directories are exclusive and must never be overwritten.

## Inputs and preparation

The exact input configuration is
`workspaces/mace_group_canonical_20260918/config.json`. It names the 25 source
preparations, canonical label inventory, unchanged model parent and plan.
Completed group preparation is `prepared_v1/preparation.json` in that workspace.
To independently replay preparation, choose the fresh output shown below:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_group_canonical.py prepare-inputs --config workspaces/mace_group_canonical_20260918/config.json --output workspaces/mace_group_canonical_20260918/prepared_review_v1
```

Model V1 had a local preflight import failure (missing two source files in the
execution snapshot); no molecular call or Slurm submission occurred. Model V2
includes those dependencies. Preserve V1 and its failure log.

## Preflight and finite execution

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_group_canonical_20260918/model_v2/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_group_canonical_20260918/model_v2/manifest.json
```

Job **1201524 completed** all 50 new endpoints; do not submit a duplicate.
For execution use the existing `diagnostics/mace_hybrid_20260916/run_pilot.sbatch`
runner and the exact absolute arguments recorded in `model_v2/submission.json`.
Inspect that receipt and live jobs before submitting: never duplicate an active
manifest. The runner executes only its 50 listed endpoints, takes an exclusive
execution lock and skips accepted receipts during recovery. A failed or missing
endpoint remains unavailable; it is never filled with a baseline value.

## Collect and report

These also work during partial execution. Calibration stays unavailable until
all 25 pairs pass numerical checks and the frozen class gap exceeds 0.02.

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_group_canonical_20260918/model_v2/implementation/mace_hybrid.py collect --manifest workspaces/mace_group_canonical_20260918/model_v2/manifest.json --output workspaces/mace_group_canonical_20260918/collection_review_v1.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_group_canonical.py report --manifest workspaces/mace_group_canonical_20260918/model_v2/manifest.json --output workspaces/mace_group_canonical_20260918/report_review_v1
```

The report includes all 28 cases, baseline fields, exact endpoint energies,
154 calibration pair comparisons, three retrospective transfers (including
unavailable 1KB0), optional research bands, and the prior grouping failures.
No broad or production-qualified classification is emitted.

Completed results are in `report_v1/result.json`: 106/154 raw class orderings,
gap −213.316367 model kcal, no bands. All 54 charge checks pass; inherited
grouping failures remain. Eight real-artifact tests pass in 407.534 seconds,
zero skipped. The saved-readout diagnostic is complete in `readout_audit_v1/`;
final costs and inspected figures are in `cost_v1/` and `figure_v1/`.
See REPORT.md and VALIDATION.md. No new inference is needed to inspect them.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_group_canonical.py -v
```

Tests distinguish archived algebra and actual source/receipt checks from new
molecular integration. Historical baseline and failed OMOL energies exercise
calibration logic only; they are never reported as new typed-group predictions.

## Figure and actual cost

After the report is available, export the figure and terminal job accounting:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/mace_group_canonical_20260918/plot_results.py --report workspaces/mace_group_canonical_20260918/report_review_v1/result.json --output workspaces/mace_group_canonical_20260918/figure_review_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/mace_group_canonical_20260918/collect_cost.py --manifest workspaces/mace_group_canonical_20260918/model_v2/manifest.json --submission workspaces/mace_group_canonical_20260918/model_v2/submission.json --preparation workspaces/mace_group_canonical_20260918/prepared_v1/preparation.json --output workspaces/mace_group_canonical_20260918/cost_review_v1
```

The cost operation refuses a live job unless `--allow-running` is explicit;
then CPU accounting remains unavailable, allocation/inference figures are
marked partial, and every failed/incomplete attempt is retained. Figure axes
use raw protocol-specific values without an invented aquo or plotting offset.
