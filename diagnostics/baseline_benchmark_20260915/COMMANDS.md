# Runnable operations

**Current status, 2026-09-16:** job1199299 completed all26 endpoints with no
failures/retries. [RESULTS.md](RESULTS.md) and [RESULT.json](RESULT.json) are the
completed records. The submission below is an archived execution/recovery
recipe; it is not a pending job or a new-analysis instruction. Current
protocol/default guidance: [agent operating guide](../../docs/AGENT_PIPELINE.md).

Run from the repository root. Existing preparation is immutable; use a new workspace path only for an explicitly intended new preparation. These commands submit no environmental/global work.

```bash
PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# Existing prepared manifest: input hashes, agreement, executable and runner checks.
"$PYTHON" scripts/affordable_workflow.py dry-run --manifest workspaces/baseline_benchmark_20260915/run_v1/manifest.json

# Execution completed as 1199299. Resubmission is unnecessary for recollection.
# The batch file runs affordable_workflow execute, then affordable_benchmark collect.
sbatch diagnostics/baseline_benchmark_20260915/run.sbatch

# Collect with strict output/receipt validation; outputs must be new filenames.
"$PYTHON" scripts/affordable_benchmark.py collect --manifest workspaces/baseline_benchmark_20260915/run_v1/manifest.json --output workspaces/baseline_benchmark_20260915/run_v1/collection_manual.json
"$PYTHON" scripts/affordable_benchmark.py report --collection workspaces/baseline_benchmark_20260915/run_v1/collection_manual.json --output workspaces/baseline_benchmark_20260915/run_v1/report_manual.md

# Real archived-artifact and preparation checks; no new QM calculations.
"$PYTHON" -m unittest discover -s tests -p 'test_affordable_development.py' -v
"$PYTHON" -m unittest discover -s tests -p 'test_affordable_benchmark.py' -v
```

Preparation operation used:

```bash
"$PYTHON" scripts/affordable_benchmark.py prepare --root . --output workspaces/baseline_benchmark_20260915/run_v1
```

The six original and six repaired site preparations are copied byte for byte. Native r2SCAN-3c/CPCM(Water)/DefGrid3, multiplicities, water states and source charge policy are frozen. The separate fixed-core 1KB0 pair uses the frozen PQQ preparation. Aquo reference outputs and receipts are verified, including the exact reference hash in the original panel. The copied manifest retains source paths/hashes and descriptors; runtime inputs differ only by the existing runner's allocation controls. Failed/partial output cannot satisfy the completion cache. A same-input retry requires a fresh output directory; no silent overwrite.

Collection preserves unrounded endpoint Hartrees and converts the contrast once using 627.509474. Repaired v3 uses the chemically compatible CN8 gauge for display, with no inherited classification threshold. Its difference from v2 can also be read as delta R without any aquo reference. Environment and response corrections are unavailable/excluded, never substituted by zero.
