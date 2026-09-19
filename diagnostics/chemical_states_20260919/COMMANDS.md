# Proton-location experiment operations

Run from repository root. Job 1202444 completed the fixed 24-task manifest;
do not resubmit it while live or rerun accepted outputs.

```bash
CHEMICAL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$CHEMICAL_PY" scripts/chemical_states_proton.py validate \
  --manifest workspaces/chemical_states_20260919/proton_paths_v1/manifest.json
"$CHEMICAL_PY" scripts/affordable_workflow.py dry-run \
  --manifest workspaces/chemical_states_20260919/proton_paths_v1/manifest.json
"$CHEMICAL_PY" -m unittest discover -s tests -p test_chemical_states_proton.py -v
```

Actual input preparation (choose a fresh output directory for a read-only remake):

```bash
"$CHEMICAL_PY" scripts/chemical_states_proton.py prepare \
  --collection workspaces/hydration_occupancy_20260918/dft_v1/collection_1201910.json \
  --agreement diagnostics/chemical_states_20260919/AGREEMENT.md \
  --output workspaces/chemical_states_20260919/reprepare_v1
```

Allocated execution uses the existing ORCA task runner through
`diagnostics/chemical_states_20260919/run_dft.sbatch`. The script requests four
simultaneous 16-rank endpoints on 64 CPUs. No GPU, numerical gradients, trajectory
or optimization. Actual job/manifest identity is in `SUBMISSION.json`.

The legacy copied analyzer closure lacked a standalone runner dependency. The
scientific runner and inputs were unaffected. A separate analyzer snapshot
`analysis_implementation_v4` validates geometry and hashes directly, avoiding that
unnecessary dependency. It can collect partial or final results without rerunning
chemistry; choose fresh write-once collection/report destinations:

```bash
"$CHEMICAL_PY" workspaces/chemical_states_20260919/analysis_implementation_v4/chemical_states_proton.py collect \
  --manifest workspaces/chemical_states_20260919/proton_paths_v1/manifest.json \
  --output workspaces/chemical_states_20260919/proton_paths_v1/recollection_v4.json
"$CHEMICAL_PY" workspaces/chemical_states_20260919/analysis_implementation_v4/chemical_states_proton.py report \
  --collection workspaces/chemical_states_20260919/proton_paths_v1/recollection_v4.json \
  --output diagnostics/chemical_states_20260919/export_replay_v1
```

All states stay in the denominator. Incomplete receipts retain unavailable values.
Raw paths and balanced metal differences are electronic diagnostics; no pH,
occupancy, entropy or released decision bands are assigned to them.

The v3 analyzer additionally reuses the already-fixed numeric energy parser from
`hydration_square.py`; the older copied parser mistook ORCA diagnostic ellipses for
numbers. A real partial collection now accepts four completed receipts and leaves
the other 20 unavailable. Both older collections/snapshots remain preserved.

Actual completed collection: `workspaces/chemical_states_20260919/proton_paths_v1/collection_1202444_recovered_v4.json`.
Actual report: `diagnostics/chemical_states_20260919/export_v1/result.json`.
