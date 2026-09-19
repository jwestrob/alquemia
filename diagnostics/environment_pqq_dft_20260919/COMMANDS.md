# Charged PQQ native-context check

Run from the repository root with the existing environment. The finite8endpoint
job1202478 is already submitted; preserve it and inspect its outputs before any
new execution. Original inputs and execution implementation are frozen in v2.

```bash
ENV_DFT_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ENV_DFT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ENV_DFT_ROOT"
"$ENV_DFT_PY" scripts/environment_pqq_dft.py validate --manifest "$ENV_DFT_ROOT/workspaces/environment_pqq_dft_20260919/prepared_v2/manifest.json"
"$ENV_DFT_PY" -m unittest discover -s tests -p test_environment_pqq_dft.py -v
squeue -j1202478 -o '%.18i %.10T %.12M %R'
```

The Slurm script automatically collects only after all requested tasks finish.
A read-only collection replay after completion uses a new result path:

```bash
"$ENV_DFT_PY" scripts/environment_pqq_dft.py collect --manifest "$ENV_DFT_ROOT/workspaces/environment_pqq_dft_20260919/prepared_v2/manifest.json" --output "$ENV_DFT_ROOT/workspaces/environment_pqq_dft_20260919/collection_replay_v1.json"
```

To reconstruct the preparation without an energy calculation:

```bash
"$ENV_DFT_PY" scripts/environment_pqq_dft.py prepare --source "$ENV_DFT_ROOT/workspaces/environment_pqq_20260919/prepared_v1/manifest.json" --agreement "$ENV_DFT_ROOT/diagnostics/environment_pqq_dft_20260919/PLAN.md" --output "$ENV_DFT_ROOT/workspaces/environment_pqq_dft_20260919/preparation_replay_v1"
```

The pinned execution command used for1202478 is recorded for reproducibility,
not an instruction to duplicate the current/completed calculation:

```bash
sbatch --output "$ENV_DFT_ROOT/diagnostics/environment_pqq_dft_20260919/dft_%j.out" --error "$ENV_DFT_ROOT/diagnostics/environment_pqq_dft_20260919/dft_%j.err" "$ENV_DFT_ROOT/diagnostics/environment_pqq_dft_20260919/run_dft.sbatch" "$ENV_DFT_ROOT/workspaces/environment_pqq_dft_20260919/prepared_v2/manifest.json"
```

Missing/nonconverged tasks remain explicit; the collector never substitutes an
old baseline value for a context endpoint. Partial jobs retain their receipts
and require a distinct, explicitly pinned retry directory under the existing
runner rules. No new absolute calibration is created by this pilot.
