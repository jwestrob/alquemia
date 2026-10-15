# Canonical solvent v2: archived-task schema repair

Preflight5135 failed with `KeyError: 'reuse'` before molecular calls: SCC started/terminal0, native component calls0. Preserve its v1 CONFIG, implementation and `preflight_5135/EXECUTION.json`. Recorded cost was79.21594091400038s and554.5115863980027 allocated CPU-seconds.

Exactly12 cached endpoints—canonical C5AXV8 and ExaA, each three geometries × two metals—use the older `context_gfn2_component_v2/manifest.json` (`7d03f585…`). Its task schema omits `reuse`; the newer REEVALUATE raw reader indexes that field directly. The underlying results remain available. No missing scientific values are supplied.

`LegacyComponentFacade` calls that archive's pinned `COMPONENT_V2.raw_results` reader, verifies exact manifest/task/RESULT/path identities, and supplies a private raw resolver to the unchanged REEVALUATE component function. Shared globals, task dictionaries and archived files remain untouched. Strict raw-output verification remains active. The legacy archived-reference ≤1e−6Eh criterion is retained alongside the unchanged component/population criteria; reconstructed results must equal their actual archived component rows.

The root-submitted Slurm preflight must pass all12 real component replays, reject copied real tasks with a corrupted ID or XYZ pin, then complete all156 cached native/component checks and the original25 pool replay. No molecular calls occur in this preflight. The other144 cached descriptors retain the existing route. Failed checks stop execution.

The binder requires every preexisting scientific configuration field to equal v1. The original PLAN, full25 denominator, ≤144 new endpoints, ≤576 SCC/288 native calls, four-term expression, gates, calibration rule and GPU receipt admission are unchanged.

Metadata-only binding, after static review:

```bash
/usr/bin/python3 diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/canonical_solvent/SOLVENT.py bind \
  --score-config workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_search_GPU_v1/CONFIG.json \
  --score-collection workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_search_GPU_v1/execution/FINAL_COLLECTION.json \
  --output workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_solvent_v2
```

Root submits the copied `RUN.sbatch CONFIG.json <actual-config-sha> preflight`; execution remains a separate root action after the actual preflight passes. No numerical tests or chemistry run on login.
