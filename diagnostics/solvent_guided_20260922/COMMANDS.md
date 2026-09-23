# Solvent-guided archive replay and common8 pilot

Use the pinned driver from the repository root:

```bash
SOLVENT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$SOLVENT_PY" -m unittest discover -s tests -p test_solvent_guided_proposals.py -v
```

Completed immutable preparations and launches:

- `workspaces/solvent_guided_20260922/archive_rankings_v2.json`: archived255source replay.
- `probes_v1/specification.json`:32declared point slots,27admissible geometries,
  five geometry-only unsupported directions; all8sources prepared.
- `pool_v1/manifest.json`:54newnativeMACE and108nativeGFN2 calls,
  SHA `82c2ab5ca9361dae184e23acdca02ee0c0b2cd2b5a505b2306c089f483e4e340`.
- GPU job1210100 and solver job1210101 completed successfully. Do not resubmit them.
  `pool_v1/SUBMISSIONS.json` pins actual wrappers, allocations and dependency.

Each path after the first is under `workspaces/solvent_guided_20260922/`.
The GPU uses32CPU, oneH200,200000MiB host memory; solvent uses64CPU/128GiB,
eight concurrent eight-rank tasks on node128 with noGPU. Existing runners
provide individual receipts and collect solver failures even on nonzero exit.
No continuous optimization, newDFT, newfold, CPCM or full225rescore is included.

The actual collection and `comparison_v1.json` are complete. To reproduce only
the parser/algebra report, use a fresh output path (no molecular calculations):

```bash
"$SOLVENT_PY" diagnostics/solvent_guided_20260922/report.py \
  --collection workspaces/solvent_guided_20260922/pool_v1/after_solvent_0_1210101.json \
  --specification workspaces/solvent_guided_20260922/probes_v1/specification.json \
  --output workspaces/solvent_guided_20260922/pool_v1/comparison_replay.json
```

This reports frozen released/adaptive reference transfer, components and selected
source-pair spread. It fits no reference for the changed proposal protocol.
Required missing energies remain unavailable; no baseline is substituted.
See `REPORT.md` for the final negative utility result and `ARTIFACTS.json` for
exact output pins. The tested proposal rule is closed; no expansion is scheduled.
