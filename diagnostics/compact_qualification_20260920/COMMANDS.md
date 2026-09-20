# Read and replay the completed numerical checks

```bash
QUAL_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
QUAL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$QUAL_ROOT"
cat diagnostics/compact_qualification_20260920/REPORT.md
"$QUAL_PY" scripts/compact_solvation_qualification.py validate \
  --manifest "$QUAL_ROOT/workspaces/compact_qualification_20260920/prepared_v1/manifest.json"
"$QUAL_PY" -m unittest discover -s tests -p 'test_compact_solvation_qualification.py' -v
```

Read-only validation launches no scientific calculation. This explicit collection
command writes a new replay record from saved output, without rerunning chemistry:

```bash
"$QUAL_PY" scripts/compact_solvation_qualification.py collect \
  --manifest "$QUAL_ROOT/workspaces/compact_qualification_20260920/prepared_v1/manifest.json" \
  --output "$QUAL_ROOT/workspaces/compact_qualification_20260920/prepared_v1/replay_v3.json"
```

Existing destinations are immutable; the replay command refuses to overwrite a
previous replay. The successful-native and failed-ordinary records are retained
together, with separate qualification fields.

The original operations are recorded in SUBMISSION.json, RESTART_SUBMISSION.json
and EXPLICIT_SUBMISSION.json. The corresponding manifest stages are
`primary_checks`, `orbital_restart_fix` and `explicit_guess_fix`. The initial
32 calls plus two eight-call technical corrections have finished. There is no
pending numerical submission and no reason to resubmit those same manifests.

For usable fresh composite scoring, use the separate
[scanner commands](../compact_scanner_20260920/COMMANDS.md). That interface accepts
explicit supported prepared pairs; prior MACE scores are optional. It retains
the original frozen context calibration and does not change production defaults.
