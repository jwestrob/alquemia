# Runnable operations

Run from the repository root on biotite. All input paths are explicit. Frozen output names cannot be overwritten; use a new output path for a further revision.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
BENCHMARK_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
```

## Inspect and verify the delivered set

```bash
cat diagnostics/benchmark_set_20260915/RESULTS.md
"$BENCHMARK_PY" -m unittest discover -s tests -p 'test_affordable_benchmark*.py' -v
"$BENCHMARK_PY" scripts/affordable_workflow.py dry-run --manifest workspaces/benchmark_set_20260915/ready_tasks_v4/manifest.json
```

## Rebuild the inventory or ledger into a fresh version

```bash
"$BENCHMARK_PY" scripts/affordable_benchmark_set.py inventory --manifest workspaces/benchmark_set_20260915/sources/structure_sources_release.json --output workspaces/benchmark_set_20260915/structure_inventory_rebuild_01.json
"$BENCHMARK_PY" scripts/affordable_benchmark_set.py assemble --evidence workspaces/benchmark_set_20260915/evidence --inventory workspaces/benchmark_set_20260915/structure_inventory_v4.json --completed-ledger diagnostics/baseline_benchmark_20260915/RESULT.json --prepared-tasks workspaces/benchmark_set_20260915/ready_tasks_v4/manifest.json --output workspaces/benchmark_set_20260915/release_rebuild_01
```

Assembly verifies pinned evidence, joins existing results/preparations and emits JSON/TSV. It does not calculate accuracy or new energies. The delivered source inventory has a sibling FASTA with deposited entity sequences. PQQ target sequences and the B2/C5 apo models are pinned in `evidence/pqq/` beneath the workspace.

## Reproduce a prepared core from the exact archived protonation

```bash
"$BENCHMARK_PY" scripts/affordable_benchmark_set.py prepare-generic --root . --config diagnostics/benchmark_set_20260915/preparation_configs/1F6S.json --protonation-report workspaces/benchmark_set_20260915/prepared/alacta_1f6s_v1/preparation_report.json --output workspaces/benchmark_set_20260915/prepared/alacta_1f6s_rebuild_01
```

Omitting `--protonation-report` performs the existing PDBFixer preparation anew; archived protonated coordinates are preferred for exact reproduction. No DFT job is launched by preparation. The analogous 6IP9 config is in the same directory.

## Rebuild a task bundle

```bash
"$BENCHMARK_PY" scripts/affordable_benchmark_set.py tasks --root . --reports workspaces/benchmark_set_20260915/prepared/hans_lanm_v1/preparation_report.json workspaces/benchmark_set_20260915/prepared/alacta_1f6s_v2/preparation_report.json workspaces/benchmark_set_20260915/prepared/alacta_6ip9_v2/preparation_report.json --agreement diagnostics/benchmark_set_20260915/AGREEMENT.md --evidence workspaces/benchmark_set_20260915/evidence --output workspaces/benchmark_set_20260915/ready_tasks_rebuild_01
```

## Subsequent scoring

The current construction phase performed **zero new endpoint evaluations**. The prepared finite campaign is ten native r2SCAN-3c/CPCM endpoints under the existing runner, without a compute budget or walltime stopping rule. Its interpretations/conditions are frozen in RESULTS.md and the task manifest; no generic threshold is assigned. This command launches that subsequent scoring campaign, so it is separate from the verification commands above:

```bash
sbatch diagnostics/benchmark_set_20260915/score_ready.sbatch
```

The script records execution and collection receipts under `ready_tasks_v4/`, using the existing restart/cache policy and allocation-aware endpoint concurrency. Collection retains raw contrasts/common-gauge values; it does not turn the ordered sites into independent labels or calibrate a new threshold. No watcher, production rescore or automatic follow-on analysis is installed.
