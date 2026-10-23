# Joint energy/force and response objective

Root prospectively declared `J=(J_old_EF+J_pair+J_C+J_chi)/4` in the parent
ROOT_DECLARATION. This adapter performs pure arithmetic. It neither evaluates a
parameter point nor chooses coefficients, domains, optimizer or stopping rules.

The exact old413-vector evaluator and body helper are passed as pinned modules.
Its returned result is retained intact, including all10 energy and15 force gates,
four0.25 weights, scales and correlated engineering residuals. No new physical
acceptance flag is inferred from the augmented loss.

Pair errors are Cartesian Debye errors divided by1D. Average their squares over
three components and all literal pairs within each parent, then average the four
parent contexts equally. W8 has8 pairs/context and W9 has9; W9 receives no extra
parent weight. C averages the four full-shell vectors and their three components,
also at1D. Chi averages four state errors divided by1a0³. Its reference is always
the h=0.001au **energy curvature**, for every state. Weights/scales are engineering
choices, not physical tolerances, variances or uncertainty likelihoods.

The weighted vector concatenates half the old413-vector, half the102 normalized
pair residuals, half the12 normalized C residuals and half the4 chi residuals:
531 entries. The returned block objectives precede the outer factor0.25. Pair
context means and every individual error are also returned. P/B, derived field
contrasts and all8 E3/mu3 cycles do not enter the loss. All28 original field rows
are retained descriptively, with exactly4 marked `used_in_objective`.

## Callable contract

`targets_from_inventory(inventory)` extracts the pinned TARGET_INVENTORY data.
The saved TARGETS.json is the ready numerical target object. Caller verifies its
pin, the exact source/helpers and data authorities; the pure function does no I/O.

`evaluate(data, model, response, targets, old_objective, body)` returns the scalar
objective, all four block objectives,531 weighted residuals, complete unchanged
old EF result, response details and descriptive field rows.

- `data`: original compact ten-cell development DATA.
- `model`: exact ten IDs, each `{energy: kcal/mol, forces: N×3 kcal/mol/Å}`; forces
  must come from the same conservative candidate energy. Missing forces fail.
- `response.pair_dipoles_Debye`: exact34 literal pair IDs, each a three-vector.
- `response.C_dipoles_Debye`: exact four parent IDs, each a three-vector.
- `response.chi_states_a0_cubed`: La_dimer, La_trimer, Dy_dimer, Dy_trimer scalars.

Response inputs are energy-conjugate physical predictions in original laboratory
frames, not preweighted errors or raw induced arrays. Pair and C quantities must
use their original literal constituent cycles. Finite-field/analytic candidate
susceptibilities must first satisfy their separately declared numerical contract;
the objective does not qualify an implementation or select a field estimator.

Existing1kcal/mol energy and2/5kcal/mol/Å force gates remain required separately.
No physical dipole/chi error limit is created. Eight consumed trimer cycles remain
regression outside this fit; reserved radial/donor/protein outcomes stay excluded.

## Bounded saved controls

One thread, no native programs or optimization. Replay the saved17676 EF result
and response baseline; check independent block sums,531-vector normalization,
equal parent weighting despite8/9 pairs, fixed estimator choice, exact-reference
zero, excluded-block independence, and missing/nonfinite/overflow rejection.
Synthetic response errors are algebra controls, not candidate parameter points.
