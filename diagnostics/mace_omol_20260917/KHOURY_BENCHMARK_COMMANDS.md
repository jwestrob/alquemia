# Khoury author-domain masked-MACE challenge commands

Read KHOURY_BENCHMARK_PLAN.md and KHOURY_NUMERICAL_ADDENDUM.md. These commands
preserve all source states and the production baseline. All output directories
are exclusive-create; use a new explicit directory for a new run.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
KHOURY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
KHOURY_WORK="$PWD/workspaces/mace_omol_20260917"
OPENBLAS_NUM_THREADS=1 "$KHOURY_PY" scripts/mace_khoury_benchmark.py prepare \
  --sources "$KHOURY_WORK/d5sc02315g_reading_v1" \
  --development-collection "$KHOURY_WORK/charge_ablation_development_v2/collection_job_1200828.json" \
  --factorization "$KHOURY_WORK/factorization_report_v2/result.json" \
  --agreement "$PWD/diagnostics/mace_omol_20260917/KHOURY_BENCHMARK_PLAN.md" \
  --output "$KHOURY_WORK/khoury_author_domains_v1"
```

The recovery operation preserves unsuccessful preparation records and first
validates existing prepared artifacts. It launches no inference:

```bash
OPENBLAS_NUM_THREADS=1 "$KHOURY_PY" scripts/mace_khoury_benchmark.py recover \
  --manifest "$KHOURY_WORK/khoury_author_domains_v1/manifest.json" \
  --output "$KHOURY_WORK/khoury_author_domains_recovery_v1"
```

Actual endpoint execution uses the existing run_pilot.sbatch executor for
each of22 explicit site manifests. The submission receipt records the copied
runner, its SHA, the finite list, Slurm command and allocation. Do not rerun
the original preparation command on an existing directory. Site manifests
support the usual `mace_hybrid.py dry-run`, `execute` and `collect` operations.

After the recorded44 endpoint tasks finish, collect all site vectors and the
predeclared nine domain-mean comparisons:

```bash
OPENBLAS_NUM_THREADS=1 "$KHOURY_PY" scripts/mace_khoury_benchmark.py collect \
  --manifest "$KHOURY_WORK/khoury_author_domains_recovery_v1/manifest.json" \
  --ggr-development "$KHOURY_WORK/charge_ablation_report_v1/result.json" \
  --ggr-2fw0 "$KHOURY_WORK/ggr_masked_2fw0_report_v1/result.json" \
  --ggr-2fvy "$KHOURY_WORK/ggr_masked_2fvy_report_v1/result.json" \
  --output "$KHOURY_WORK/khoury_author_domains_report_v1"
```

These descriptors are model units, not binding free energies. No PQQ band or
new absolute classification is supplied. Any unavailable site makes the
primary corresponding domain mean unavailable; failed entries remain visible.
