# Relaxed alternative amide: complete matched composite result

**All10 sources and40 new solvent endpoints completed**, with80 original endpoints reused. The expanded common6 pool changes only Ca-conditioned sample4 and La-conditioned sample4; the other eight retain their original minima. Both mathematical and operational selectors give the same result. Here Ca4/La4 describe source conditioning, not biological class.

The alternative amide followed by the unchanged four donor modes accesses lower composite states in both sample4 sources. This supports a missed basin in the current model. It does not establish experimental accuracy, a uniquely causal H-contact energy, or a calibrated replacement classifier.

## Saved response and coverage

`R = E(Ca) − E(La)`. Positive changes mean relatively greater lowering of the La endpoint. All entries below are copied from the actual4971 collection, in kcal/mol.

|Source|Common3→4 ΔR|Common4→6 ΔR|Common3→6 ΔR|MACE contribution3→6|Frozen-solvent contribution3→6|
|---|---:|---:|---:|---:|---:|
|Ca4|+10.974999511|+9.908598533|+20.883598044|+17.988344487|+2.895253557|
|La4|+0.673057892|+18.617905858|+19.290963751|+15.992516580|+3.298447170|
|Other8|0|0|0|0|0|

Both final sample4 pools select `relaxed_flip_Ca` for Ca and `relaxed_flip_La` for La. Relative to common3, their selected total energies fall by7.583891025/28.467489069 for Ca4's Ca/La endpoints, and1.411029366/20.701993116 for La4's Ca/La endpoints. These are finite shared-pool minima, not a claim of a stationary composite optimum.

The fixed-flip solvent benefit of roughly9.5kcal/mol was shared by both sources. After donor relaxation, the *additional* common4→6 differential solvent terms oppose the vacuum contribution: −6.618764700 for Ca4 and −6.258081843 for La4. Thus the final response cannot be attributed simply to a uniquely favorable solvent correction in the previously inconclusive fold.

## Fold spread and explicitly old-band sensitivity

Saved-only job5030 reproduces each original common3 score against its exact frozen-distribution record, then summarizes all ten sources. No band is fitted.

|Source group|Common3 range|Common4 range|Common6 range|
|---|---:|---:|---:|
|Ca-conditioned5|23.325984036|12.350984525|5.677901201|
|La-conditioned5|20.459434524|19.786376632|2.222540672|
|All10|30.128351231|19.786376632|12.480268395|

Ranges are kcal/mol; mathematical and operational results agree. The remaining conditioning-associated separation is visible in the all10 range.

Using only the **old common3 canonical bands** (`Ca_max = −405461.9631022772`, `La_min = −405453.5127788069`) gives9 correct/1 inconclusive at common3 and10 correct at common4/common6. Ca4 is the sole changed call. This is an old-band transfer sensitivity on consumed sources. **Calibration of the alternative-amide protocol remains unavailable** until the same policy competes fairly across the canonical reference; no production promotion follows.

## Numerical admission and physical limits

Actual4971 preflight replayed all20 final geometry/full-source admissions and all reused common4 evidence. The saved5030 audit passes all120 component-v2 records, all four native water-diagonal checks per cell, informative vacuum-charge bounds, and identical water/vacuum cavity and Born-matrix hashes. Largest saved component-parity error is1.492362003e−7Eh; population discrepancy1.149115054e−7e; diagonal absolute residual3.409454341e−8Eh; vacuum-charge rounding bound2.687731597e−7Eh. Existing gates were unchanged. All240 historical absolute endpoint-parity flags remain failed: component-v2 admission does not erase those failures. Saved algebra closes within6.650324735e−11kcal/mol against the existing1e−8 tolerance.

The [named-contact analysis](../contact_review/FINAL_REPORT.md) shows HD22–Gly292H relief in **both** sample4 sources while retaining oxygen-facing amides. Included short H contacts remain, and the inherited ND2–H lengths remain approximately1.183–1.188Å. Ca4/La4 searches are not boundary-limited, but four other searches reach the original0.8Å bound. The result remains conditional on the fixed dry204-atom state, cofactor/proton preparation, frozen-charge ALPB/cavity model, and MACE-only constrained proposals. It is not a full-protein equilibrium or affinity result.

## Artifacts and cost

- [4971 collection](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_solvent_v1/execution/FINAL_COLLECTION.json), SHA256 `64668c7ec758da42d34b039855bdc370b39d06815db47a00dfaf69c41a69368d`.
- [5030 saved summary](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_solvent_summary_v1/execution/RESULT.json), SHA256 `2d5ab6278f90b0705e96727399a5fc60ff440d2c6077ef34213f6492c7d31649`; [all scores/choices/calls](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_solvent_summary_v1/execution/SCORES.tsv), [ranges](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_solvent_summary_v1/execution/RANGES.tsv), [editable SVG](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_solvent_summary_v1/execution/COMMON_POOLS.svg), [PNG](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_solvent_summary_v1/execution/COMMON_POOLS.png).
- Plot offset is the single fixed common3 canonical A0AC La3 score, −405429.1999088685kcal/mol, for every source and pool. It does not independently center each curve.
- Job4971: COMPLETED0:0,942s×7CPUs=6594 allocated CPU-s; actual160 SCC and80 native replays, zero new MACE/DFT/optimizations. Internal wall941.196263s. Job5030: COMPLETED0:0,3s×1CPU=3 allocated CPU-s, zero molecular calls; internal summary/plot wall2.202171s. Accounting is scheduler allocation cost, not measured CPU utilization.

This report reads saved results only. All arithmetic/plotting was performed by5030; no numerical analysis ran on login.
