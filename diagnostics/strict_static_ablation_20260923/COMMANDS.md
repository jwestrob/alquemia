# Actual read-only operations

Run from the repository root. The first command has already completed; do not
overwrite its reference. Full transfer comparison awaits the owned strict225
collector. Each operation refuses to overwrite its output.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
scripts/strict_static_ablation.py calibrate \
--qualification workspaces/strict_native_pool_20260923/run_v1/COMPARISON.json \
--output workspaces/strict_static_ablation_20260923/run_v1/REFERENCE.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
scripts/strict_static_ablation.py compare \
--transfer workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json \
--reference workspaces/strict_static_ablation_20260923/run_v1/REFERENCE.json \
--output workspaces/strict_static_ablation_20260923/run_v1/COMPARISON.json

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
-m unittest discover -s tests -p test_strict_static_ablation.py -v
```

Reference SHA256:
`8f87bb14f5d138c7dfda4b3622ed055d019930323f3db60b3bcf79a6874e5741`.
Static bands: Ca_max `-405467.7497230548`, La_min `-405456.6749943662`
model kcal/mol. The25-member canonical gap is11.07472869; this is calibration,
not a finding of better transfer. No crystals, folds or PLM inputs were fitted.

Three actual-artifact tests pass in0.144s, zero skips. The initial test discovery
used the wrong module for the existing eV conversion constant; correcting that
import resolved the test failure. No scientific executable ran for this analysis.
