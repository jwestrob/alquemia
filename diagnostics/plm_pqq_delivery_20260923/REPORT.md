# PLM PQQ operations and manuscript package

**The practical candidate now has a fresh-source preparation bridge, an SOP,
batch export and paper draft.** Jacob can start explicit PLM batches using the
[SOP](../../docs/PLM_PQQ_SOP.md). No cohort calculations were submitted here.
The scientific scorer, frozen reference, released static default and explicit
DFT access are unchanged.

## What was already complete

The three-source 4.3 Å envelope scorer, full reference transfer, two actual PLM
triples, source response figure and biological joins existed. Its execution API
accepted general complete prepared triples, but operating instructions started
from archived protonation. That was a real gap for newly folded PLM proteins.

## What changed

- `scripts/pqq_plm_prepare.py` validates explicit source requests, invokes the
  released seeded source preparer when fresh protonation is needed, and passes
  the resulting checkpoint to the existing envelope/scoring planners. Exact
  archived preparation remains an explicit option. It runs no endpoint energy.
- `params/pqq_plm_envelope_v1.json` pins the existing reference, release and
  numerical qualifications. It introduces an operations identity, not a new
  scientific score or calibration.
- The small batch wrapper prepares and invokes the existing executor under its
  established H200/32 CPU/200000 MiB profile. It does not reimplement optimization,
  scoring, numerical retries or collection.
- Every declared preparation group stays in `READY.json`. Complete triples can
  execute while failed triples remain unavailable; no reduced-member median.
- `scripts/pqq_plm_export.py` combines disjoint real result batches through the
  existing exact-sequence biological join. Original columns remain unchanged;
  failures and unscored inputs remain distinct. Duplicate protein results and
  mixed references are rejected.
- [SOP](../../docs/PLM_PQQ_SOP.md), [main-text draft](MANUSCRIPT.md) and
  [technical supplement](SUPPLEMENT.md) provide commands, interpretation,
  evidence denominators, costs, limitations and actual output paths.

## Actual verification

Fresh preparation was executed on the two previously consumed PLM triples,
six real source structures. All six succeed. Every final endpoint coordinate,
atom identity/order, charge, multiplicity, physical-motion map and state/protein
signature exactly matches the archived envelope preparation. This verifies the
fresh-source bridge without inventing energies or rerunning successful scoring.

Source/envelope preparation took **68.097288 seconds** measured local wall time;
the six existing source-preparer calls account for **50.748385 seconds** of that.
The CPU protonator is pinned to one thread. This was local preparation, with no
Slurm allocation or measured allocated-core receipt; CPU time and peak memory
were not separately measured. Scoring-plan creation and tests are additional.
No GPU, DFT, new fold, new solvent or new MACE energy evaluation ran in this task.

The resulting v2 execution plan passes the actual existing preflight with
12 MACE origins and 24 native GFN2 origin inputs; its prospective full batch is
12 searches, up to 12 cross-MACE calls and 72 scalar cells. None was launched.
Do not call this a fresh end-to-end molecular benchmark.

Five new real-artifact checks pass in 7.810 seconds, zero skips. They cover exact
fresh preparation, unexecuted complete preflight, original biological rows and
actual scores, duplicate-result rejection and explicit preparation failure
without a fabricated energy. Thirteen prior envelope API/integration regression checks also pass in
38.674 seconds, zero skips (`REGRESSION_TESTS.txt`). Python compilation and shell syntax pass.

The actual export CLI reproduces the previous 176-row, 53-column protein table
byte-for-byte (SHA256 below), retaining only the two actually scored candidate
groups; the other 174 remain unscored. It performs no expression normalization.
The SOP's combined submission wrapper is syntax-checked but was not submitted;
its preparation and lower scoring stages have separate actual evidence.

## Actual products and next operation

- Fresh preparation/request/preflight:
  `workspaces/plm_pqq_delivery_20260923/fresh_two_v1/`.
- Actual historical-score export:
  `workspaces/plm_pqq_delivery_20260923/export_v1/proteins.tsv`, SHA256
  `3a6eadbd22485c97e6f5c37bc9ff4bd2db9b7558f8422dc82fc21206a28a5cd1`.
- Original computed results:
  `workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json`.
- Editable scientific figure:
  `workspaces/plm_envelope_response_20260923/figures_v2/`.

Immediate read-only verification:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/pqq_three_source_envelope_execution.py dry-run \
  --plan workspaces/plm_pqq_delivery_20260923/fresh_two_v1/prepared/scoring/plan.json
```

For new proteins, start at SOP step 1. Its existing-cohort source builder requires
the real AF3 input contract and verified role mapping. New PQQ families need
their actual source manifests and biological crosswalks; the old 176-row XoxF
table does not silently become an all-PQQ cohort. Unsupported chemistry remains
unavailable. No source selection, label or band is changed to obtain a desired call.

The benchmark claims remain those of the completed reports: practical primary
91/100 correct complete groups and nine unavailable, separate qualified recovery
94/100 correct and six unavailable; strongest tenfold result 207/0/1/17 among
225 structures is a different context policy. All comparisons are development
evidence. Unknown PLM preferences remain predictions. This closes packaging of
the current La/Ca round while LanM data organization proceeds separately.
