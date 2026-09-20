# Completed reference-grid replay

Run from the repository root. No molecular calculation needs repeating.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_accommodation_reference_profile.py -v
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/accommodation_reference_profile_20260920/prepared_v1/implementation/accommodation_reference_profile.py analyze --manifest workspaces/accommodation_reference_profile_20260920/prepared_v1/manifest.json --mace-result workspaces/accommodation_reference_profile_20260920/prepared_v1/mace_1203879/result.json --low-collection workspaces/accommodation_reference_profile_20260920/prepared_v1/low/collection_1203880.json --output workspaces/accommodation_reference_profile_20260920/result_replay_v2.json
```

The `prepare --help` operation exposes all source, mapping, agreement, comparison
and q0 solvent manifest paths. The exact case/angles are intentionally fixed for
this experiment version. Actual submitted commands are in `SUBMISSIONS.json`.
Outputs are write-once; keep the original experiment and use fresh replay paths.
