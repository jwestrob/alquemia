# Runnable replay

Run from the repository root; outputs are write-once. Completed chemistry should
not be rerun. The commands below validate and reanalyze existing actual results.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_accommodation_response.py -v
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/accommodation_response_20260920/analysis_implementation_v1/accommodation_response.py validate --manifest workspaces/accommodation_response_20260920/prepared_v2/followup/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/accommodation_response_20260920/analysis_implementation_v1/accommodation_response.py analyze --gate-collection workspaces/accommodation_response_20260920/prepared_v2/gate/collection_1203465.json --followup-collection workspaces/accommodation_response_20260920/prepared_v2/followup/collection_1203499.json --output workspaces/accommodation_response_20260920/result_replay_v3.json
```

Preparation and bounded execution are implemented via `prepare` and `execute`;
inspect `--help` for explicit required paths. The completed manifests already
contain all scientific inputs, software pins and the finite 60-task denominator.
Their original execution receipts are recorded in the two submission JSON files.
