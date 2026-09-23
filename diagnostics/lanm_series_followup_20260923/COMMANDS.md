# LanM follow-up operations

Use the pinned adapter for actual execution. Paths below assume the repository
working directory; shell variables explicitly name the installation/workspace.
Do not rerun completed tasks: submissions and execution receipts are immutable.
The capability endpoint is excluded from the remaining tasks, preserving the
12MACE/48GFN2 total. See SUBMISSIONS_CAPABILITY.json,
SUBMISSION_CAPABILITY_CONTINUATION.json and SUBMISSIONS_REMAINING.json for the
exact actually submitted argument arrays.

```bash
lanm_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
lanm_python=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
lanm_work="$lanm_root/workspaces/lanm_series_followup_20260923/prepared_v2"
lanm_adapter="$lanm_work/implementation/lanm_series_followup.py"
lanm_authorization="$lanm_root/diagnostics/lanm_series_followup_20260923/EXECUTION_AGREEMENT.md"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

"$lanm_python" "$lanm_adapter" validate --manifest "$lanm_work/manifest.json"
```

`validate` performs no molecular evaluation. The source/paired-state preparation
is preserved in sources_v2 and prepared_v2; sources_v1 is the rejected alternate
conformer attempt, and prepared_v1 is an earlier unexecuted code snapshot.

## Native continuation

The adapter stages only actual matching successful initial GBW+xtbw. Its
`--task` option may be repeated for the explicit remaining23 IDs in
SUBMISSIONS_REMAINING.json. It must exclude `Hans_EF1__Dy__vacuum`, whose
continuation is already in native_capability_continuation. No source editing is
required to set paths or select the finite task subset.

For any newly declared, unrun prepared native manifest the bounded executor is:

```bash
# Run only inside its approved allocation, with that exact manifest argument.
"$lanm_python" "$lanm_adapter" execute-native \
  --manifest "$lanm_work/native_remaining_continuation/manifest.json" \
  --authorization "$lanm_authorization"
```

The supplied run_native.sbatch requests64CPUs/128GiB,8ranks per cell and noGPU.
Capability native jobs override to8CPUs/32GiB for their single cell. The warm
MACE wrapper requests exactly oneGPU, one task with32CPUs and200000MiB. All use
the previously supported explicit CPU-sharing GPU node and existing ORCA/Python.

## Collect the ordered comparison

After both native continuation manifests and MACE result files exist:

```bash
"$lanm_python" "$lanm_adapter" collect \
  --manifest "$lanm_work/manifest.json" \
  --mace-collection "$lanm_work/mace_capability/result.json" \
  --mace-collection "$lanm_work/mace_remaining/result.json" \
  --continuation-manifest "$lanm_work/native_capability_continuation/manifest.json" \
  --continuation-manifest "$lanm_work/native_remaining_continuation/manifest.json" \
  --output "$lanm_work/final_collection.json"
```

This is a write-once collection operation. Read final_collection.json once it
exists; do not overwrite it. Missing/failed cells remain null, and all12endpoint
and3ordered-site denominators remain present. No implicit aquo-reference value,
zero correction, primary-SCF fallback, fitted label or calibrated band is used.

## Real-artifact verification

```bash
"$lanm_python" -m unittest discover \
  -s "$lanm_root/tests" -p test_lanm_series_followup.py -v
```

Tests exercise actual prepared inputs and archived native results; they do not
generate synthetic molecular outputs. The original unrun-status test refers to
the earlier immutable prepared_v1, while executable preparation is prepared_v2.
