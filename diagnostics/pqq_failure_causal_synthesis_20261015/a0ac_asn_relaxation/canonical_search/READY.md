# Canonical MACE/search stage: scout complete; CPU fallback bound

Root submitted the scout as5057, with completion wakePID4022277. Inspect
asn_canonical_search_scout_v1/SUBMISSION.json and current scheduler/receipts
before acting. Do not duplicate this scout; remaining23 still require its
actual completed result and root review.

Root's CPU scout5057 completed with both searches and common6 admitted; see [SCOUT_REPORT.md](SCOUT_REPORT.md). The immutable remaining23 CPU configuration is bound as an **unsubmitted fallback**, while root evaluates the separate GPU route. No jobs were submitted by this agent.

Executed scout command, retained for provenance; do not repeat:

```bash
sbatch workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_search_scout_v1/RUN.sbatch \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_search_scout_v1/CONFIG.json" \
  e77354e51ca5d742330475d8a4c713bdc1ec1f278d5b35a21a8c8a88864513cd
```

The wrapper requests seven CPUs on standard, with no host, memory, time, account or QOS override. Root owns submission and completion wake. The scout is A0A3F2YLY8 La-conditioned sample1,202 atoms,Ca−3/La−2. Its original selected four modes exclude anchor-Asn χ2; that remains unchanged. It dispatches two fixed-flip calls, two bounded searches and up to two cross scores. All25 rows and the exact already executed A0AC result remain in the collection.

| Artifact | SHA256 |
|---|---|
| Scout `CONFIG.json` | `e77354e51ca5d742330475d8a4c713bdc1ec1f278d5b35a21a8c8a88864513cd` |
| `RUN.py` | `9da0c95d4077b1c413f19bcf0d0307b6fa00fe5947de768f8329b75dde86a2bb` |
| `RUN.sbatch` | `77a53ad8dd70bf7d4bdfed7c0bf64a216aab16c8f963a64075944336d71865a4` |
| `PLAN.md` | `c85ae5730267b1648281ed19deed41df17647d278e84593e69fd863bb2682d58` |
| Actual5016 `PREPARATION.json` | `663ad699c2fae8baf6a42ee65078481a26c872d0fa6ad0d6cf17d3e8eeb0fbca` |
| Existing frozen canonical collection | `ee00fa03aba64c23e6a737fb2b3b4a9f40dc139e58bdb9225e8ed10c9af2b58d` |

## Remaining23 CPU fallback — bound, not submitted

Actual remaining config SHA: `59fab20f7f615249cd73b9c6d091e79d31b737ecc8ae03bc1c5dbcd3551fea67`. Worker manifest SHA: `525fc9ef9b9a490146565e9b9072d58329fc0b22fdf80aa00a0ac47de320920b`. The implementation and plan are byte-identical to5057. If root chooses this route, its exact command is:

```bash
sbatch workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_search_remaining_v1/RUN.sbatch \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_search_remaining_v1/CONFIG.json" \
  59fab20f7f615249cd73b9c6d091e79d31b737ecc8ae03bc1c5dbcd3551fea67
```

Do not submit concurrently with another route covering these sources. The following metadata binding command has already completed and must not be repeated into the same workspace:

Binding required the actual successful scout collection **and** execution receipt. It performed no scientific execution and submitted no job:

```bash
/usr/bin/python3 diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/canonical_search/RUN.py bind \
  --preparation workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_relaxation_v1/execution/PREPARATION.json \
  --stage remaining \
  --scout-collection workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_search_scout_v1/execution/FINAL_COLLECTION.json \
  --output workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_search_remaining_v1
```

The emitted config pin supplies the exact final argument to its copied `RUN.sbatch`. Remaining23 means46 fixed-flip calls,46 original searches and at most46 cross calls, with the scout and A0AC reused. Full new scope is48/48/48; line-search forwards vary under unchanged maxiter200. Rough3–4hour single-lane planning estimate is based on measured A0AC costs, not a cap or guaranteed runtime.

## Preserved evidence and limits

Slurm preflight verifies all150 archived cells and their actual force arrays/count/state; historical H200 records retain their original provenance. Actual qualified CPU fixtures/runtime and SciPy pins are checked before any call. Fresh flip forces qualify each new cached q0 against its rebased Jacobian. Full-source guards retain inherited short contacts and reject newly introduced severe contacts, using existing thresholds. Missing required cells block common6; no partial minimum or altered label-based gate is used.

Outputs are `execution/PREFLIGHT.json`, per-source fixed/proposal/cross receipts, `Q0_QUALIFICATION.json`, progress receipts, `FINAL_COLLECTION.json` and `EXECUTION.json`. Collection reports common3/4/6 mathematical and operational vacuum works only. No SCC/DFT/solvent, class decision or calibration is executed. The eventual exact25 calibration rule and separate old-band interpretation are frozen in PLAN.

A historical warm GPU route is being evaluated separately. It is not a direct transport replacement because request authorization, receipt schema and cleanup differ. This CPU package remains unchanged; no claimed GPU speedup or duplicate canonical chemistry is embedded here.
