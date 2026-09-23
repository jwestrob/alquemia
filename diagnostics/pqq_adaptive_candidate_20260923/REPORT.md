# Minimal adaptive PQQ candidate is usable through an opt-in source command

**Delivered and exercised the complete source-to-score path.** Fresh 1H4I and
4MAE preparations, paired native forces, unchanged four-mode searches, complete
three-geometry scoring and frozen-reference decisions all succeeded. The command
makes the previously demonstrated fold-robustness gain usable on explicitly
selected compatible sources. It does not change the released default.

The broader existing result remains a development benchmark: the recovered
minimal candidate has204correct/1wrong/1inconclusive/19unavailable among225
correlated fold samples. On common206 sources, released static scoring has
202correct/2wrong/2inconclusive. This integration adds operational evidence;
its two consumed crystals are not independent biological validation.

## Actual fresh integration

Job1210488 completed normally, exit0. Fresh request/plan/implementation were pinned
before execution. Both preparations reproduced archived context atom ordering and
coordinates exactly, including H. Charges:1H4I Ca−2/La−1;4MAE Ca−3/La−2; singlets.
No chemistry, assembly, source selector, water inventory or reference was changed.

| Source | Original R | Candidate R | Candidate call | Difference from archived candidate R |
|---|---:|---:|---|---:|
| 1H4I | −405506.647982 | −405505.425780 | Ca-supported | −0.00572905 |
| 4MAE | −405441.666967 | −405444.304425 | La-supported | −0.00000000565 |

R is Ca-minus-La in model kcal/mol. Original and candidate use their own frozen
bands. The large common offset is not binding free energy; no aquo score or
probability is invented. Maximum origin composite discrepancy was3.14×10⁻⁹
kcal/mol. Final proposals reproduce archived coordinates within3.91×10⁻¹⁴A;
the small1H4I scalar difference is retained as native numerical variability,
not interpreted as a physical gain. All four searches converged without a final
boundary flag. No class-directed pose or energy-based source selection occurred.

Actual calls: **72 native MACE** (4origin,64search,4cross), **24/24 native GFN2**,
zero failures and zeroDFT. Eight candidate-native energies reuse their own
appropriate proposal/origin calculations through pinned state/coordinate checks;
there is no hidden old-score fallback. No extra starts, retries or cohort rescore.

## Interface and checks

New command: `scripts/pqq_adaptive_candidate.py` with prepare, dry-run, execute,
collect and report. Actual requests and commands are in COMMANDS.md. The command
reuses released source preparation, source-connected Kinematics, force projection,
four-mode selector, scaledSLSQP200, warm OMOL worker, ORCA executor and pool algebra.
It leaves the old hardcoded recovery script and every running/historical input intact.

Fresh sources accept only existing explicit selected-crystal-chain or Protenix
PQQ normalization/configuration. Preparation defects and unsupported chemistry,
fewer than four independent modes, absent forces, failed proposals or incomplete
energy matrices remain unavailable. Raw Ca source conventions are not silently
normalized beyond the released policies. New source domains remain conditional;
physical/protocol compatibility does not create a biological label.

Archive replay is separate and currently accepts explicit completed origin-recovery
bundles. Both actual recovery cases reproduce their scores exactly with zero
molecular calls. This is intentionally a narrow reuse interface, not a generic
archive migration tool. Exact pinned artifacts and the existing executor's safe
partial-attempt behavior provide reuse; failed jobs are not automatically rescued.

**Eight focused tests pass, zero skipped.** Tests use actual recovery forces,
completed pools, released preparation and the fresh integration. They verify
selector replay, coordinates/state/recipes, frozen bands, no missing-cell fallback,
unsupported normalization, changed method invalidation and actual fresh outcomes.
Malformed-input tests explicitly corrupt copies of real records. Earlier test
interface failures were corrected before submission and remain in test logs.

## Measured cost

Slurm allocation: **220s elapsed ×32CPUs =7040 allocated core-seconds**, one H200
reserved for **220GPU-seconds**,200000MiB requested host memory. This includes
preparation, launching workers, native solvers, staging and in-job collection.
There is no matched-hardware speed claim.

Measured components (included above, not additional charges): source preparation
21.753s; origin/candidate ORCA stages44.943/89.028s; warm native workers4.699s
(origins),26.923s(searches),4.746s(crosses), including3.786s combined model loading.
Worker import/launch, validation and file/report overhead remain in total allocation.
Peak CUDA allocation6,055,696,384bytes; batch MaxRSS3,588,468KiB (not an aggregate
whole-job memory peak). Local preparation/tests/reporting outside Slurm are additional
unmetered CPU work. Reservation time is not GPU utilization.

## Status / next operation

The candidate is **opt-in**, scientifically unchanged from the supported minimal
pool. Fresh eligibility beyond these tested inputs remains conditional. Native
GFN2 convergence/force limitations remain documented by the earlier continuation
and force studies; this command adds no solvent-gradient optimization or new
numerical solver. No general affinity claim, default promotion or reference refit.

Read the actual result without rerunning chemistry:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/pqq_adaptive_candidate_20260923/two_crystals_v1/implementation/pqq_adaptive_candidate.py \
  collect --plan workspaces/pqq_adaptive_candidate_20260923/two_crystals_v1/plan.json \
  --output workspaces/pqq_adaptive_candidate_20260923/two_crystals_v1/review_collection.json
```

Authoritative result, immutable execution/solver/native receipts, preparation,
verification and scheduler accounting are under
`workspaces/pqq_adaptive_candidate_20260923/two_crystals_v1/`; compact pins are in
ARTIFACTS.json. Parent independently reviewed the adapter and accepted the
source/paired-state, common-pool and missing-result handling. No work remains
running in this delivery.
