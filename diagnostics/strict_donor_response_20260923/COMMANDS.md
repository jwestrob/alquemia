# Exact donor-response diagnostic

Job1211802 owns the finite32 fresh calls; do not duplicate it. Four strict4MAE origin cells are reused from the completed strict32 fresh branch. Production is unchanged.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RUN="$PWD/workspaces/strict_donor_response_20260923/run_v1"
"$PY" "$RUN/implementation/strict_donor_response.py" validate --manifest "$RUN/manifest.json"
"$PY" -m unittest discover -s tests -p test_strict_donor_response.py -v
```

The wrapper collects once to `COLLECTION.json`. Read its36 `rows`,18 `points`,12 `endpoint_works` and6 `differentials` without another molecular run. A technical recollection, if required, must use a new write-once output path and retain the first result. The original native DFT/MACE comparison remains pinned through the manifest.
