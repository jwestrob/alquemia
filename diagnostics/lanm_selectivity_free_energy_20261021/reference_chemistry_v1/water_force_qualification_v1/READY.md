# Ready for root review: matched water force panel

**48 cells are prepared; no native calculation or submission was made.** Use `prepared_v2/PANEL.json` and `READY.json`. The pre-review source and its prepared_v1/validation_v1 artifacts are preserved; all 48 final decks and coordinates remain byte-identical.

- Source: `2871d929c576281fba152de0a7f58069e8c8e88068cd173f970c11fe0952db7c`.
- Panel: `c11f1b12e5f70b92a32c8912b31d3c3e10a275d2d2cc7a226f432c5590cf0cd2`.
- Validation: `310a7022e3fd0d25aca09af64bdd12044cee9bb876cc75ab2fe8d4d0bcae27a9`.
- Final static preparation: 1.013 s; saved validation: 3.726 s, one local thread. All commands/results are recorded in PREFLIGHT.json and VALIDATION.json. Validation passes 12 actual native output replays and 65 synthetic controls; zero molecular calls.

## First native cell

`cc__La__plus_010__complex`: canonical CCSD(T), 18 correlated electrons, four centers/90 orbital functions, four MPI ranks. Exact launch arguments are in READY.json. Default scheduler memory/time and the established 75%-of-grant MaxCore rule are retained. This largest water inventory precedes the remaining 47 cells; the completed scout must be reused, not repeated.

## What the collector will resolve

Both fixed step sizes report complex and counterpoise derivatives, separate PBE0 electronic/D4 contributions, model residuals and Dy−La contrasts. Individual and paired derivatives have separate 0.1 kcal/mol/Å two-step checks. The smaller-step summary is withheld when that check fails. Saved PBE0 analytic forces are compared directly with uncorrected PBE0-D4 finite differences; neither is equated to a counterpoise derivative. Central/displaced final-integration and COSX signatures are compared explicitly, with raw native stages retained.

This is a same-basis frozen-f method diagnostic. It establishes no physical truth, basis convergence, free energy, parameter fit or protein-sampling approval.
