# Frozen pilot operations

From the repository root. Prepared inputs already exist; do not rerun preparation
or completed chemistry to replay a report. Runtime paths are recorded in manifests.

```bash
SS_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
SS_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
SS_PREP=$SS_ROOT/workspaces/second_shell_20260919/prepared_v2
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# Read-only preparation and manifested task checks.
"$SS_PY" "$SS_PREP/implementation/second_shell_context.py" validate --manifest "$SS_PREP/manifest.json"
"$SS_PY" "$SS_PREP/implementation/second_shell_mace.py" validate --manifest "$SS_PREP/mace_manifest.json"

# Replay collection into fresh files only, once actual execution is complete.
"$SS_PY" "$SS_PREP/implementation/second_shell_context.py" collect \
  --manifest "$SS_PREP/manifest.json" --output "$SS_PREP/dft_recollection_v1.json"
"$SS_PY" "$SS_PREP/implementation/second_shell_mace.py" collect \
  --manifest "$SS_PREP/mace_manifest.json" --output "$SS_PREP/mace_recollection_v1.json"

# Compare actual calculations, with separate DFT/OMOL energy scales.
"$SS_PY" scripts/second_shell_report.py \
  --dft "$SS_PREP/collection_1202082.json" \
  --mace "$SS_PREP/mace_collection_1202083.json" \
  --output "$SS_ROOT/diagnostics/second_shell_20260919/replay_v1"

"$SS_PY" -m unittest discover -s tests -p 'test_second_shell_context.py' -v
```

Jobs1202082/1202083 were submitted with the exact scripts and manifests pinned in
SUBMISSIONS.json. Their execution uses the existing ORCA manifested runner and
native OMOL worker/receipt checker. A small stage adapter keeps this experimental
manifest out of shared dispatch files. Source snapshots are immutable. Failed
prepared_v1 reached only JSON regeneration verification, not execution; no inputs
were changed to resolve its JSON integer/string key comparison.

No failure is replaced with baseline data. MACE partial attempts remain visible;
restarting reuses only fully accepted receipts. Failed tasks require an explicit
fresh unchanged attempt. This small pilot does not fit an absolute reference or
classification threshold. The output supplies matched contrasts and components.

## PQQ physical-H result replay

All four proposals and all four DFT checks are complete. The existing results
must not be overwritten. These commands create fresh read-only replays:

```bash
SS_H=$SS_ROOT/workspaces/second_shell_20260919/pqq_h_v1
SS_H_DFT=$SS_ROOT/workspaces/second_shell_20260919/pqq_h_dft_v1
"$SS_PY" "$SS_H/implementation/second_shell_hydrogen.py" validate \
  --manifest "$SS_H/manifest.json"
"$SS_PY" "$SS_H_DFT/implementation/second_shell_hydrogen.py" collect-dft \
  --manifest "$SS_H_DFT/manifest.json" --output "$SS_H_DFT/recollection_v1.json"
```

`PQQ_H_SUBMISSION.json` and `PQQ_H_DFT_SUBMISSION.json` pin the exact launch inputs.
The declared final-iterate rule includes constrained stationary proposals even
when SLSQP reaches its iteration limit. No third optimization, relaxed boundary
or score-dependent proposal selection was introduced. See `PQQ_H_REPORT.md`.
