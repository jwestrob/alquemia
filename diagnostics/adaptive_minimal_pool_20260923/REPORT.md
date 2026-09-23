# Preserve the adaptive benefit with fewer candidate geometries

**The older terminal-only proposals are unnecessary for every available class
decision in the completed reference-fold test.** A common pool containing only
the origin and the two four-mode adaptive proposals retains the earlier gain:
202 correct, one wrong and one inconclusive on204 matched structures, versus
200/2/2 for released static scoring. This is a simplification of a demonstrated
development benefit, not a new biological accuracy gain.

The reduced method independently reproduces the same canonical25 reference
limits as the five-candidate method. All25 calibration and three consumed crystal
calls remain correct. The two PLM examples retain unknown biological labels.
Every available individual, strict-group and three-fold decision agrees with
the five-candidate method. Eleven individual raw scores change, by at most
1.823652 model kcal/mol; no score is silently treated as identical.

## Full coverage and cost meaning

| Method | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Released static |203|2|2|18|
| Five-candidate adaptive |202|1|1|21|
| Three-candidate adaptive replay |202|1|1|21|

The existing unavailable preparations/searches remain unavailable. Their origin
scores do not substitute for a missing adaptive result. Removing the dependence
on failed older proposal jobs is a separately authorized implementation/recovery
task, not a result of this replay.

Three distinct geometries require12 native GFN2 singlepoints instead of20, and
only the two adaptive searches are needed. That is **40% fewer nominal solvent
calls**, not a measured40% wall-time speedup. This experiment reused every
molecular result and ran zero new MACE, GFN2 or DFT evaluations. It adds no fresh
numerical-convergence qualification; the historical native-solver limitations
remain distinct from the concurrently completed continuation tests.

The gain remains specific: one corrected error, two resolved inconclusives and
one new inconclusive versus released static scoring on common coverage. The
known A0A3F2YLY8 Ca-sample1 error remains. Available structural aggregates were
already correct; this simplification adds no new aggregate accuracy.

## Method and reproducibility

Protocol: `Nikasha_origin_and_two_adaptive_candidates_replay_v1`.
Both metal rows use exactly `{origin, adaptive_Ca, adaptive_La}`. The archived
native OMOL plus native GFN2(ALPB−vacuum) components, physical mappings, charges
and chemical states remain intact. Mathematical row minima and the operational
0.1 kcal/mol origin-retention rule are recorded separately. Synthetic candidate
counts are not populations. No absolute aquo reference or affinity is inferred.

Calibration uses only the original25 canonical structures and the established
class-extrema/minimum-gap rule, frozen before the reduced225-fold replay.
Its limits happen to equal the earlier adaptive limits:
Ca maximum−405459.0155077257; La minimum−405456.46881754015 model kcal/mol.
No fold outcome changed them. Released-band transfer remains199correct,
0wrong,5inconclusive,21unavailable; it is distinct from the new-reference calls.
All these structures are consumed development evidence, not fresh blind tests.

`scripts/adaptive_minimal_pool.py` supplies separate calibration and comparison
operations using existing matrix, reference and strict-aggregation machinery.
The compact [RESULT.json](RESULT.json) pins the full records and lists all eleven
raw changes. Six real-artifact tests pass, including missing required cells,
corrupted labels/reference rejection, label-independent geometry selection and
strict missing-member propagation. No successful molecular output is fabricated.

**Next:** use the smaller adaptive pool as the candidate implementation and
recover the two valid-origin sources excluded by obsolete proposal dependencies.
Keep production unchanged while coverage and practical execution are completed.
See [commands](COMMANDS.md). No promotion, PLM rescore or remote push occurred.
