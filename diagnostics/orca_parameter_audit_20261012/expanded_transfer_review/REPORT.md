# Independent corrected expanded transfer review

12 October 2026. **The corrected expanded candidate retains every supported
primary triple decision and improves three individual development cases over
its own static origins. One Ca-conditioned La-class source remains wrong.**
This is useful, limited accommodation benefit; it does not establish broad
Ca-conditioned robustness or uniformly reduce structural spread.

## Exact population and common-coverage comparison

The frozen canonical calibration is identical to `CANONICAL_COLLECTION.json`,
including all three policies, extrema and bands. No transfer values entered its
construction. The operational and mathematical canonical gaps remain
3.472171047 kcal/mol; the independently calibrated origin gap remains
10.297678912 kcal/mol. Both protocols below use their own frozen canonical bands.

| Population | Origin: correct / wrong / inconclusive / unavailable | Accommodated: correct / wrong / inconclusive / unavailable |
|---|---:|---:|
| Canonical calibration, 25 | 25 / 0 / 0 / 0 | 25 / 0 / 0 / 0 |
| Crystal transfer, 3 | 3 / 0 / 0 / 0 | 3 / 0 / 0 / 0 |
| Prescribed probes, 6 | 2 / 2 / 2 / 0 | 5 / 1 / 0 / 0 |
| Additional fold/context pools, 101 | 101 / 0 / 0 / 0 | 101 / 0 / 0 / 0 |
| All noncalibration pools, 110 | 106 / 2 / 2 / 0 | 109 / 1 / 0 / 0 |
| All prepared pools including calibration, 135 | 131 / 2 / 2 / 0 | 134 / 1 / 0 / 0 |
| Primary unique source/context pairs, 104 | 103 / 1 / 0 / 0 | 104 / 0 / 0 / 0 |
| Declared three-source summaries, 100 | 94 / 0 / 0 / 6 | 94 / 0 / 0 / 6 |

The first four rows partition the 135 pools; subsequent rows are overlapping
summaries, not additive independent tests. The primary 104 pairs comprise the
101 additional pools plus three reused probes, from 98 unique structures in
25 biological groups. Their exact pair membership matches the original transfer
inventory. All 100 original triples and six prior preparation exclusions remain;
no favorable subset was selected. Every prepared pool is available for both
origin and accommodation, so the comparison has common coverage.

Accommodation changes two A0A3F2YLY8 Ca-conditioned probes from inconclusive to
correct and the A0ACD6B9F2 La-conditioned sample-4 probe from wrong to correct.
Its Ca-conditioned sample-4 counterpart remains wrong despite a
+12.901442 kcal/mol shift in the La-like direction. There are no new wrong calls
among these pools. This is not a comprehensive Ca-conditioned test: all 135 pools
contain only three Ca-conditioned sources, alongside 129 La-conditioned sources
and three crystal entries without a conditioning label. All examples have been
used during development; canonical fit and correlated replicas are not new
independent predictive validation.

The earlier failed native scalar cell remains an archived failure. Its newly
available standalone-backend result belongs to this separate corrected protocol;
it does not overwrite or retroactively pass the historical native run.

## Structural spread

Among the same 94 supported triples, accommodation decreases the source range in
51 and increases it in 43. Median range changes from 5.104313 to 5.522710 kcal/mol.
The largest decrease is 16.366158 kcal/mol in an A0ACD6B9F2 triple; the largest
increase is 6.582809 kcal/mol in a Q92WY9 triple. All available median class calls
are unchanged. These overlapping triples are correlated. At the protein grouping
level, ranges over all primary source/context pairs decrease for 12 groups and
increase for 13; that supplementary grouping can also mix context selections.
There is no general reduction in structural sensitivity established here.

## Independent raw checks

The audit imports no scoring-runner code. It reparses actual xTB stdout/charges
and tblite scientific output, verifies receipt/output hashes, source state and
archived MACE energies, and independently rebuilds row minima, origin retention,
calibration extrema, decisions and triple medians. It checks all 1620 input
records (3240 backend evaluations including reuse) and all 810 solvent components.
All reconstructed selections/decisions agree; score residuals are below the
1e-7 kcal/mol comparison tolerance.

All 810 component-v2 gates pass. Maximum independent solvent-component error is
2.124440925e-7 Eh (limit 1e-6); maximum all-atom population discrepancy is
5.411258489e-7 e (limit 1e-5). Finite/convergence, electron/charge closure and
one-thread receipt checks pass. **All 1620 original absolute-energy comparisons
still fail v1**, with maximum absolute residual 1.963768749e-4 Eh. The previously
documented non-SCF dispersion-convention difference is not fitted away or called
absolute-energy agreement.

## Actual cost

| Stage/backend | New calls | Measured child CPU seconds | Sum of call wall seconds |
|---|---:|---:|---:|
| Canonical xTB | 276 | 1979.057753 | 1988.662059 |
| Canonical tblite | 276 | 8106.917830 | 8120.915695 |
| Transfer xTB | 1320 | 9509.873372 | 9554.426605 |
| Transfer tblite | 1320 | 39743.427008 | 39809.425436 |
| Total fresh | 3192 | 59339.275963 | 59473.429796 |

Canonical numerical execution lasted 744.688672 seconds; transfer lasted
3563.979321 seconds. Both used 14 one-thread workers. These are stage timings,
excluding preparation, collection and the interval between stages. An additional
24 input pairs, or 48 backend calls, were reused with no new molecular work.
The maximum recorded individual-child RSS was 129020 KiB, not an aggregate peak.
No GPU or Slurm allocation was used. All MACE evaluations/proposals were reused;
these timings do not measure a complete new-protein scanner run. The audit made
no molecular calls.

## Reproduce and next decision

```bash
python diagnostics/orca_parameter_audit_20261012/expanded_transfer_review/REVIEW.py
```

Final collection SHA-256:
`49bcb0040f31d7f19864591c2392364cb847123ed113095476cdd56c1aa39fab`.
`RESULT.json` retains every source, triple, selected work and unavailable result.
The justified next question is transfer to the missing Ca-conditioned population
and practical source-to-score cost. This result alone does not authorize promotion
or erase the separate shared-pool and representation sensitivities.
