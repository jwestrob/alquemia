# Precision34 retains all calls and recovers both optimizer failures

**34/34 sources are available and correct under both the old reference and the
new canonical25-only reference.** This includes25/25 canonical controls,3/3
consumed crystals,4/4 noncanonical development sources and2/2 previously failed
searches. The stopping change is useful for runtime and coverage. It does **not**
meet strict numerical equivalence on every source:27/32 matched pools pass all
frozen gates. Keep it a separately calibrated candidate; production and the
historical225 ledger remain unchanged.

The only search change is ftol1e−8 Hartree-equivalent. All source coordinates,
q0 forces, charges, inventories, four-mode selections, one-start policy,200
iterations, physical0.8Å displacement and0.8rad bounds are unchanged. Four
already-qualified precision pools were reused exactly;30 new sources ran60
searches and the finite common-pool comparison. No source was chosen after
seeing its new score, and no failure was retried.

## Classification and numerical consistency are separate results

| Prespecified check | Actual result |
|---|---|
| Valid new endpoints | 60/60, plus8 reused endpoints |
| Matched successful-old native-energy gate,0.001kcal/mol | 66/66; maximum0.0000905442 |
| Matched GFN2 medium-cell gate,0.1kcal/mol | 250/256 |
| Matched pooled-R gate,0.2kcal/mol | 29/32 |
| All gates jointly | 27/32 |
| Existing-reference calls | 34/34 correct |
| Own canonical25 calibration | 25/25; gap7.325388850kcal/mol |
| Own-reference consumed transfers | 3crystals+4noncanonical+2recovered,9/9 correct |

The old canonical gap was7.240239674. The small gap change is not evidence of
added physical information. New bands (both mathematical and operational):
Ca≤−405463.71230032056; La≥−405456.3869114709 modelkcal/mol. Only the25 designated
canonical sources determine them; all other sources remain outside calibration.
All data are consumed development; structural samples are not independent
biological observations. No within-lanthanide or measured-affinity claim follows.

Five sources fail at least one frozen numerical gate:

| Source | deltaR, kcal/mol | Largest absolute medium-cell difference | Failed gate |
|---|---:|---:|---|
| C5AXV8 | +0.081906 | 0.104415 | cell |
| I0JWN7 | −0.906285 | 1.054521 | cell andR |
| A8R3S4 | +0.202652 | 0.235520 | cell andR |
| P16027 | −0.000057 | 0.117166 | cell |
| Q60AR6 canonical | +0.417686 | 1.038771 | cell andR |

All their classifications remain correct. Their native MACE endpoint differences
are below0.001kcal/mol, while native-GFN2 medium differences are larger; this
locates the limiting numerical sensitivity in the solvent-scoring calculation,
without establishing a unique SCF-branch mechanism. No values were smoothed,
restarted or selected to make these gates pass. The largest actual mapped-atom
change among fresh matched successful endpoints is0.0049545Å. Candidate choices
change for I0JWN7La, P38539La and canonicalQ60Ca.28/68 endpoints retain a boundary
flag; bounded candidates are not asserted unconstrained minima or populations.

## The actual failures are recovered

| Source endpoint | Old iterations/evaluations | New iterations/evaluations | New search wall | New R/call |
|---|---|---|---:|---|
| Q92WY9 Ca-conditioned sample2,La | 200/2061,failed | 13/13 | 3.589s | −405447.264232,La-supported |
| Q60AR6 La-conditioned sample0,La | 200/2052,failed | 15/17 | 3.641s | −405485.637602,Ca-supported |

Both satisfy the same final geometry gate and retain boundary flags. Their old
successful-pool comparisons stay null. Relative to the actual **unsuccessful**
old final iterates, mapped-atom differences are0.00007096Å and0.000000740Å,
with native-energy differences0.000000120kcal/mol and0 respectively. These are
supporting diagnostics, not fabricated successful old reference results.

## Cost and execution

All240 new nativeGFN2 endpoints and all new MACE calls succeed. The60 new searches
use907 optimizer function evaluations versus5355 old, and summed search wall
230.700s versus1430.420s (the latter includes both failed old searches). The907
optimizer evaluations include60 reused q0 evaluations: **847 new search MACE
calls+60 fixed cross-MACE calls=907 new molecular MACE calls**. No new q0, DFT,
training, additional starts or alternate stopping tolerances were run.

GPU1210813:350s×32CPUs,1requestedGPU,200000MiB. CPU1211010:666s×64CPUs,128GiB,
eight8-rank native workers. Total **53824 allocatedCPU-seconds/350 requestedGPU-
seconds**. This includes startup, validation, preparation and collection inside
jobs. Local setup/tests/reporting are unmetered. The four reused precision pools
previously cost8480CPU-s/75GPU-s and are excluded from this new-work total; all
older q0/method-development costs remain historical reuse. Existing repeated
manifest validation is appreciable overhead, distinct from molecular search time.

## Records and recommendation

Full comparison, all component differences, actual selections and failed gates:
`workspaces/slsqp_precision_expansion_20260923/COMPARISON_v1.json`, SHA
c276b5f9ad35d7411f6b50e9a19670588607cd738b331696d41363a601c42b95.
New `REFERENCE_v1.json` SHA
1f8470bbdf7a056098ca261aa10b72bdc940cf71ab860fbbc54090bda1e73250.
`SEARCH_AUDIT_v1.json`, `COSTS.json` and all raw execution receipts remain pinned.
Nine actual-fixture tests pass in36.889s (zero skips) and cover source membership, original failure retention,
unchanged method/selector, exact four-case reuse, mapped displacement algebra,
canonical-only reference and the failed equivalence gates. No fake science or
executable substitutions were used.

Proceed to a complete structural-transfer evaluation of this separately
calibrated numerical candidate. Do not silently replace the old protocol or
call the solvent scoring numerically identical. No additional molecular work
was launched after these results; the possible one-rank runtime requires its
separate qualification. The repaired-H MMOL1770 supplement remains separate.
