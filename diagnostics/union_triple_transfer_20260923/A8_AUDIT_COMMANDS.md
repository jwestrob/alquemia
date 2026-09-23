# Reproduce the read-only A8 audit

From the repository root, use a new output path; existing artifacts are immutable:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/union_triple_transfer_20260923/audit_a8.py \
  --comparison workspaces/union_triple_transfer_20260923/COMPARISON_v1.json \
  --continuing-inventory workspaces/precision_pool_continuation_20260923/INVENTORY_v2.json \
  --output workspaces/union_triple_transfer_20260923/A8_AUDIT_replay.json
```

This reads all four declared A8 triples, existing raw energies, proposal angles,
SCF traces and retained seeds. It runs zero energies or forces. The saved v2
artifact is `workspaces/union_triple_transfer_20260923/A8_AUDIT_v2.json`.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  -m unittest discover -s tests -p test_union_triple_a8_audit.py -v
```

The source-distance audit likewise reads saved coordinates only:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/union_triple_transfer_20260923/audit_a8_geometry.py \
  --audit workspaces/union_triple_transfer_20260923/A8_AUDIT_v2.json \
  --output workspaces/union_triple_transfer_20260923/A8_GEOMETRY_replay.json
```

The saved artifact is `workspaces/union_triple_transfer_20260923/A8_GEOMETRY_v1.json`.
