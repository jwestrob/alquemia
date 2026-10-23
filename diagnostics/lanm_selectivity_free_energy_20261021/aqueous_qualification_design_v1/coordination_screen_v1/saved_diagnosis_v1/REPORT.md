# W9 failure predates fitting; the force tradeoff worsened

The original model already fails the balanced W8→W9 energy test. Fitting improves
that error by0.941981 kcal/mol but worsens both individual W9 force errors while
substantially improving W8. This supports testing the two coordination geometries
together in development; it does not yet demonstrate that the functional form
cannot represent them or identify a unique missing physical interaction.

## Balanced energy and native components

All entries are kcal/mol for `(Dy W9−La W9)−(Dy W8−La W8)`.

| Quantity | Original pre-repair | Frozen17605 | Fit-induced change |
| --- | ---: | ---: | ---: |
| Total model contrast | −6.427418 | −5.485437 | +0.941981 |
| Reference contrast | −0.534207 | −0.534207 | 0 |
| Model−reference error | −5.893211 | −4.951230 | +0.941981 |
| Native van der Waals contribution | −6.977397 | −6.067869 | +0.909528 |
| Native polarization contribution | +0.549979 | +0.582432 | +0.032453 |
| Permanent multipoles and bonded terms | 0 | 0 | 0 |

The printed native components close to the total within their eight-decimal
precision. Most of this *model contrast* is the van der Waals term. There is no
matching unique decomposition of the quantum reference, so this accounting does
not prove that the physical discrepancy is solely repulsion, dispersion or a
missing polarization effect. Radius, well depth and damping changed together.

## Rigid-force transfer

Entries are category-balanced RMS / largest body-component error, kcal/mol/Å.
Existing gates remain RMS≤2 and maximum≤5 separately for each row.

| Geometry and comparison | Original | Frozen17605 |
| --- | ---: | ---: |
| W8 La | 3.083 / 9.441, fail | 0.453 / 1.066, pass |
| W8 Dy | 2.066 / 6.156, fail | 1.336 / 3.813, pass |
| W8 Dy−La | 1.028 / 3.286, pass | 1.018 / 3.226, pass |
| W9 La | 2.725 / 7.754, fail | 3.811 / 7.961, fail |
| W9 Dy | 2.047 / 6.551, fail | 4.686 / 9.949, fail |
| W9 Dy−La | 0.968 / 2.006, pass | 1.054 / 2.268, pass |

The W9 error involves the whole shell. After fitting, each water's outward
rigid-net-force error is positive: La4.398–5.970 and Dy6.513–8.027 kcal/mol/Å.
The first eight waters contribute summed outward errors+39.722/+56.578 for
La/Dy; the ninth contributes+5.090/+6.513. The largest candidate body component
is the metal z-force error,−7.961La/−9.949Dy. Before fitting, the per-water radial
errors were negative. These are force projections at fixed coordinates, not
energy-component forces or equilibrium radial pressure.

Thus the fit shifted radial forces outward at both sizes: this repaired W8 much
more successfully than W9. Shared errors still cancel in the Dy−La force
comparison, explaining why that passing check cannot rescue the individual tests.
The original and fitted models each pass only one of the four W9 screen gates.

## Smallest next falsifiable development test

Use the same six water parameters in one separately declared finite joint test
on the **ten already-consumed pure-water cells**: water-dimer minus/central/plus,
W8 central and W9 central, each for La and Dy. Preserve every existing individual,
paired and energy gate; adding W9 must not discard the previous radial or W8
requirements. All ten have compatible uncorrected ECP-PBE0-D4 references under
the same per-metal input recipe, charge+3 and effective multiplicity1; physical
La1/Dy6 and ECP cores46/55 remain explicit. No new quantum campaign is necessary
to ask whether this candidate family can satisfy these development data jointly.

The first matched six-column parameter Jacobian can itself diagnose whether W9
residuals point into directions that would spoil the original eight cells. Twelve
central probes affect five cells of the perturbed metal each: at most120 native
programs before exact cache reuse. If a finite joint fit is chosen, reuse that
initial Jacobian in it rather than launch a separate duplicate diagnostic. No
objective, fit configuration or candidate is bound here. A finite failed search
would not prove global non-representability.

W9 also rearranges the original eight waters: their oxygen z coordinates change
from±1.25Å to±0.875Å. Neither this saved comparison nor a joint fit can separate
water count from that geometry change. A later matched-removal control could
separate those effects if the scientific question requires it; no such geometry
or calculation is prepared here. The existing W9 references remain provisional
frozen-f DFT without a W9-specific finite-difference force qualification.

## Evidence and cost

Eight model records were replayed from actual native outputs: four original and
four candidate, all with exact typed coordinates, bonds, parameter bytes and
native executable/parser lineage. The original key is retained as the candidate
key prefix. The four admitted quantum reference records are reused.57 saved
checks passed in0.946s wall/0.905s CPU on one thread; zero native calls. The
metadata-only ten-cell inventory took0.099s wall/0.078s CPU. An earlier saved
report attempt failed only on an absent inventory field after its arithmetic
checks; source and failure receipt are preserved in `attempt_v1` (1.117s tool
wall time, CPU not measured).

Workspace artifacts: `saved_diagnosis_v1/RESULT.json` (`f5cb2cb5`) and
`INVENTORY.json` (`54b64bb5`). Analysis source `ANALYZE.py` (`e340244b`), declared
specification (`dc22887b`). The actual17670 failure, all completed fit results,
current donor work and reserved donor outcomes remain unchanged.
