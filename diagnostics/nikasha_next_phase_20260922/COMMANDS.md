# Reproduce the completed comparisons

The existing molecular jobs are finite and recorded in their workspace submission
receipts. Do not submit them again. These commands only collect/compare real saved
outputs. The output names are new because reports deliberately reject overwrites.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
NIKASHA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
NIKASHA_POOL=workspaces/adaptive_completion_20260922/primary225_pool_v1

"$NIKASHA_PY" "$NIKASHA_POOL/implementation/nikasha_pool.py" collect \
  --manifest "$NIKASHA_POOL/manifest.json" \
  --output "$NIKASHA_POOL/recollected_v1.json"

"$NIKASHA_PY" scripts/nikasha_pool_compare.py compare \
  --collections "$NIKASHA_POOL/recollected_v1.json" \
  --reference workspaces/adaptive_completion_20260922/original30_pool_v1/REFERENCE.json \
  --prior-comparison workspaces/nikasha_recovery_20260922/proposal_comparison.json \
  --output "$NIKASHA_POOL/recompared_v1.json"

"$NIKASHA_PY" scripts/nikasha_next_phase_compare.py \
  --context workspaces/consistent_context_20260922/primary225_v1/comparison_v1.json \
  --adaptive "$NIKASHA_POOL/recompared_v1.json" \
  --output workspaces/nikasha_next_phase_20260922/RECOMPARISON_v1.json
```

The final command writes both a full JSON ledger and a readable Markdown table.
It retains all225 sources,75 strict group summaries and100 correlated La triples.
Each candidate uses its own canonical-only reference; no threshold is fitted on
these fold outputs. Neither command changes the released scorer or historical
DFT access. The context and CPCM branches retain their own runnable operations in
their diagnostic directories.

## Implemented candidate preparation

`scripts/nikasha_pool.py prepare-expansion --proposal-family scaled-angular`
accepts the existing original30 or primary225 proposal collection. The exact
executed arguments, validated task counts, pinned implementation and scheduler
commands are retained in `primary225_pool_v1/manifest.json`, `PREFLIGHT.json`,
`SCHEDULER_PREFLIGHT.json` and the five `SUBMISSION_*.json` files. The same
preparation operation refuses incompatible chemical states, origin mappings or
incomplete candidate pairs. Execution uses the existing shared-pool MACE runner
and native solvent shards; there is no separate workflow system.

The released source-to-score command and explicit DFT mode remain documented in
[the production commands](../pqq_fast_release_20260920/COMMANDS.md).

## Cost interpretation

Five distinct admitted geometries require20 nativeGFN2 singlepoints per source
(five geometries × two metals × vacuum/ALPB), versus four for the released
single-geometry composite. Exact-coordinate deduplication can reduce this.
The completed base three-geometry pool is reused in this phase; the new two
adaptive candidates therefore need eight additionalGFN2 calls and two additional
cross-MACE calls per admitted pair, plus their recorded proposal searches.
Incremental development receipts are not a measurement of a cold, full
source-to-score production run. Reused work, failed searches and the local
collection/preparation overhead remain distinct in the final cost report.
