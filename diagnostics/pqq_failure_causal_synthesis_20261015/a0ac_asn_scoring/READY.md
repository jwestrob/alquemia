# A0AC Asn290 scorer readiness

The scorer is ready to bind the actual ten-source preparation. No scoring
configuration or molecular job has been created by this adapter yet.

`SCORE.py` SHA256:
`7682a6265084039a1acc663bf4812cfff3aebd7b3acbf454a0f6642ab1701388`.
AST and Bash syntax checks passed. Independent static review found no remaining
blocker. Scientific preflight and model execution remain **UNRUN**.

## Scope and actual archive compatibility

The frozen parent protocol is `../ASN_ROTAMER_PLAN.md`. Preserve all ten sources,
sixty original cells, and twenty declared flip cells. Only paired preparations
with actual geometry admission are dispatched. Each complete new pool contains
origin, adaptive Ca, adaptive La and the single amide flip under both metals.
Mathematical and strict −0.1 kcal/mol operational selections remain separate.
The sign-blind ≥1 kcal/mol stage gate never changes class bands or launches
solvent automatically. Missing new cells leave common4 unavailable, while the
complete common3 baseline remains separately reported.

Metadata inspection of the sixty actual receipts found 24 complete CPU results,
12 computed H200 results with manifest-task provenance, and 24 complete H200
results with request provenance. Their model dictionaries match the qualified
worker; the CPU runtime dictionaries also match. Historical GPU receipts retain
their actual provenance. No GPU record is relabeled as CPU validation.

Before any new model call, Slurm preflight validates all sixty real 204-atom
states, requests and finite 204×3 saved force arrays. New calls use the unchanged
seven-thread warm CPU worker and its strict state/version/forward checks. The
existing actual Ca/La qualification receipts are reused; no new qualification
calls or optimizer starts are added. Parent numerical threads remain one.

The actual preparation inventory path is retained and checked against fixed
content SHA `c91b152b1fb8b548395c41aac72fb2370b5ba73d015abbfce14767ac2787d4d8`.
The diagnostic inventory copy is independently pinned as the reference. The two
paths are not required to be identical.

## Metadata binding after actual preparation

From the repository root, once the named actual preparation exists:

```bash
python diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_scoring/SCORE.py prepare \
  --preparation workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_preparation_v1/execution/PREPARATION.json \
  --output workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_scoring_v1
```

This command only verifies/copies metadata and writes a new immutable configuration.
The resulting snapshot wrapper accepts `CONFIG.json`, its printed SHA256, and
`run` as three arguments. Root owns its standard-partition seven-CPU submission
and completion wake. The `run` action performs the full scientific preflight
before loading the model. A separate `preflight` action is available but is not
required as a duplicate test job.

Execution uses an exclusive workspace lock and refuses existing output paths.
Partial results precede cleanup; cleanup or collection failures remain recorded.
The execution receipt survives resource-reporting errors, and those errors cause
a failed process exit. Original archives and failed geometry rows are preserved.
