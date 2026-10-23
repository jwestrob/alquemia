# Actual dry-acetate component review — job18733

**Pass:** 16/16 raw analyze records and all18 declared numerical checks replay independently. The original admission helper ran once; independent Decimal arithmetic and source/input/receipt checks passed534 checks. No native calls, repeat controls, source changes or quantum/current-fit reads were made.

Largest absolute total FD-minus-analytic derivative: `0.00671648272379` kcal/mol/Å. Largest two-step total change: `0.005038` kcal/mol/Å. Both are below the unchanged0.1 resolution. Maximum raw component-sum remainder: `2E-8` kcal/mol; derivative component-sum remainder: `0.000001` kcal/mol/Å; all remain inside their propagated printing bounds.

## Retained derivative values

All entries are dE/dq in kcal/mol/Å; generalized forces have the opposite sign. The JSON retains both steps, every component, signed differences and print bounds.

| Motion / channel | h=.010 | h=.005 | Saved analytic |
|---|---:|---:|---:|
| radial / La | 98.178206500 | 98.183174000 | 98.184829964 |
| radial / Dy | 119.508327000 | 119.513365000 | 119.515043483 |
| radial / paired | 21.330120500 | 21.330191000 | 21.330213519 |
| rotation / La | -20.387989000 | -20.390410000 | -20.391218258 |
| rotation / Dy | -7.316291500 | -7.316093000 | -7.316026701 |
| rotation / paired | 13.071697500 | 13.074317000 | 13.075191557 |

At h=.005, the model Dy-minus-La radial derivative is `21.330191 = 14.709231 vdW + 6.620960 polarization`; the rotational derivative is `13.074317`, with `7.147166 vdW + 5.927152 polarization` (1e−6 printed closure). Permanent-multipole differences and rigid intraligand component derivatives are printed zero in this panel. Individual rotation totals are La−20.390410 and Dy−7.316093.

These are properties of the selected17621 Hamiltonian and this motion. Centroid rotation changes donor distances as well as orientation; a large model component is not a causal decomposition of the difference from quantum forces. No physical accuracy or unique mechanism is qualified. Native induction residuals were not printed under the unchanged key/command.

## Costs and provenance

Scheduler: COMPLETED, exit0:0, 11seconds ×1CPU =11 allocated CPU-seconds. Application: 8.89150seconds wall, 0.97824seconds parent CPU, 2.794175seconds child user CPU and3.885449seconds child system CPU; native elapsed sum7.42574seconds. These measures have different scopes and are retained separately.

REVIEW.json binds exact actual result, releasea002c79d, sourcebfafefed, manifest45311d30, original analyzef67420e7, source review71dab138, all raw output pins through admitted rows, and an immutable completion-event snapshot.
