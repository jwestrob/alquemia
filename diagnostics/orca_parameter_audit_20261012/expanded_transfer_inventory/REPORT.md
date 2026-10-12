# Expanded candidate: corrected reevaluation is ready

The existing 4.3 Å envelope candidate has enough real saved inputs for a finite
corrected calibration and structural-transfer test. No new geometry search is
needed. Root has launched the frozen canonical stage; completed scientific
results belong in its execution collection, not this preparation report.

## What is being tested

Each source retains its published three-candidate pool: origin, Ca proposal,
La proposal. Both metals see that same pool. All 810 archived MACE cells and
all exact coordinate/state/map references were checked. The published strict
−0.1 kcal/mol origin-retention rule and canonical-only extrema calibration
remain unchanged. Replaying all 134 complete historical pools gives exactly
zero score residual. The one historical missing native cell remains recorded.

The new variable is the corrected GFN2 solvent component, evaluated with pinned
standalone Linux xTB 6.7.1 and independently checked by isolated tblite 0.7.0.
The earlier absolute-energy parity gate remains failed; its geometry-only ATM
convention difference cancels in the matched water-minus-vacuum component.
No fitted correction is subtracted. Component and atom-charge gates remain
1e-6 Eh and 1e-5 e. Source charges vary and are preserved: Ca/La −4/−3, −3/−2,
and −2/−1, all supported singlets.

## Frozen population and staging

| Stratum | Pools | Interpretation |
|---|---:|---|
| Designated canonical calibration | 25 | 14 Ca / 11 La consumed references |
| Crystals | 3 | 1H4I, 4MAE, 1KB0; transfer only |
| Prescribed development probes | 6 | Known-reference folds, separate stratum |
| Additional fold/context pools | 101 | Complete declared transfer inventory |

The latter 101 plus three reused probes supply 104 source/context records from
98 physical structures, grouped into the original 100 triples. All 100 remain
in the denominator: 94 prepared, six historical preparation exclusions. Pools
are correlated within proteins. The full 135-pool inventory spans 26 biological
groups because 1KB0 contributes Q46444 beyond the 25 training groups.

All six probes have existing reference-class labels; no unknown PLM or reserved
labels were opened. Coverage is strongly La-conditioned: 129 La, three Ca,
three crystal pools. This does **not** establish full Ca-conditioned transfer.
The three Ca probes are A0A3F2YLY8 samples 1/3 and A0ACD6B9F2 sample 4.

24 corrected real input pairs (one xTB and one tblite evaluation per input)
are reusable from the completed context factorial, covering the two canonical
edge proteins and their own three candidates. Other corrected archive summaries
lack raw receipts and are not reused. The separate two-anchor common-five
experiment is not mixed into this candidate's calibration.

- First: 276 new xTB + 276 new tblite calls for the 25 canonical pools.
- Only if all 25 qualify and the published operational extrema gap exceeds
  0.02 kcal/mol: 1,320 additional calls per backend for the remaining strata.
- Maximum declared new work: 1,596 calls per backend; no new MACE, DFT, folding,
  optimization, water inventories, or threshold tuning.

Failure of calibration stops dependent transfer. It is not a reason to tune
bands on the alternate structures. Origin-only diagnostics retain their own
availability and calibration; missing proposal cells still invalidate the
primary pooled result.

## Execution, integrity and recovery

Final workspace is `workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3`.
Manifest SHA256:
`4f842a2d92364e1559dd5d266ae86e6b682a8ed81181a9e58f0a17cbd59ccfe1`.
Root launched 14 one-thread workers and an actual completion-queue watcher;
see `CANONICAL_LAUNCH.json`. This uses the user's current explicit login-thread
authorization. It is not a scheduler allocation template.

The driver checks executable/library/parameter pins, source settings and state,
MACE receipt identity, exact XYZs, raw process commands/stdout/stderr, and reparsed
energies/populations. Python-version differences in summation produce about
1e-13 e differences in derived electron sums; only those derived sums allow
1e-10 e comparison tolerance. Scientific admission gates are unchanged.

A lock prevents duplicate concurrent stages. Numbered immutable attempt receipts
allow untouched inputs to resume; completed inputs are never recomputed.
Interrupted point directories remain explicitly unavailable, with no automatic
retry. A surviving child process blocks resumption. Read receipts before any
restart. Draft v1/v2 preparation artifacts were superseded before molecular
execution and remain preserved.

## Cost expectation, not a new measurement

The preceding 80-input paired-backend factorial used 2,446.44 CPU-seconds over
160 calls. Extrapolating that mean gives about 8,440 CPU-seconds for calibration
(10.0 minutes at ideal 14-worker occupancy) and 48,806 CPU-seconds for the full
conditional test (58.1 ideal minutes). Using the measured 192-atom rate instead
gives roughly 12.8 and 73.8 ideal minutes. New inputs span 169–218 atoms, beyond
the prior largest 192-atom measurement; SCF variation and workload imbalance
make these planning estimates only. Actual receipts measure wall/CPU/RSS.

## Remaining evidence gap

The supplied October snapshot excludes the full corrected current34/current39
raw archive. Its two edge summaries are insufficient to reconstruct that result.
This test therefore explicitly recovers the real September23 published candidate
population. It does not rerun the 2,796-input original panel or assume missing
corrected outputs existed. Production defaults and prior results stay intact.

Recommendation: finish the canonical stage, then apply the unchanged calibration
gate. Proceed to the declared structural transfer only if it passes.

## Preparation checks

Two real-artifact tests passed: full 135-pool scope/state/reuse accounting, and
exact replay of 134 completed historical pools plus deliberately unavailable
cell handling. Independent read-only review confirmed canonical identity,
14/11 labels, all 100 triple memberships, runtime/reuse gates and published
selection/calibration policy. Preparation/testing ran no molecular evaluations.
