# Saved pair/collective energy cancellation

Declared before arithmetic. Four consumed parents only: central W8/W9, La/Dy.
Use frozen selected46 aa06677d, evaluator baseline/static manifest8f29dde6,
literal34 pair coefficient dictionaries from targets61d73a16, admitted common-
basis collectiond60a2f8d and its exact original D4 authorities9ec1b009.
The eight-trimer result0c017d17 motivates this diagnostic; its selected-candidate
join is retained. No fresh outcomes, native calls, fitting or selection.

For model energy E in kcal/mol, use the exact existing fixed-coordinate cycles:

```
C = E(M+shell) - E(M) - E(shell)
P = sum_i [E(M+Wi) - E(M) - E(Wi)]
B = C-P
```

Retain all four baseline and selected C/P/B values and selected-minus-baseline
changes. The literal coefficient dictionary subtracts the ion N times in P.
All unchanged ion/water/shell model energies are saved static terms.

The common-parent-basis electronic C/P/B is available for every parent. Absolute
per-parent **electronic-plus-D4 targets and their total-energy errors are not
available**: the real-shell D4 term was not calculated, and ghost shell/water
decks deliberately omit D4. Do not set missing terms to zero or compare model
total energies against electronic-only references as total-energy errors.

The exact Dy-minus-La differences at each N are available. Their pure-shell
and pure-water D4 terms cancel because the two parents have identical literal
water coordinates and chemistry. Keep all ghost-basis electronic terms; they
do not cancel by assumption. For these paired quantities only, add once:

```
Delta C_D4 = D4(MDy+shell)-D4(MLa+shell)
Delta P_D4 = sum_i[D4(MDy+Wi)-D4(MLa+Wi)] - N[D4(Dy)-D4(La)]
Delta B_D4 = Delta C_D4-Delta P_D4.
```

Use the existing reference conversion627.509474 kcal/mol per hartree. Display
electronic, D4 and their sum separately. For paired W8, paired W9 and their
original coordination difference W9-minus-W8, report model-minus-reference
errors and the exact relation error(C)=error(P)+error(B). Opposite signed P/B
errors describe cancellation; no new physical threshold is introduced.

Recheck only input/source joins, full denominators, literal geometry matches,
cycle closure and replay of the old paired aggregate. No duplicate raw native
admission or broad test suite. Record one-thread wall/CPU cost. Retain the
provisional ECP-PBE0 basis/state/SCF limitations and the unavailable per-parent
D4 targets. B contains all orders beyond independent pairs; it is not uniquely
a three-body term, charge transfer, polarization, or an isolated ninth-water
insertion. Passing the fitted25 E/F gates need not mean accurate absolute C/P/B.
