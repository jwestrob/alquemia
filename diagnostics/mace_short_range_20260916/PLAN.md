# Whole-protein short-range learned descriptor

Under the active MACE discriminator goal, declare this component experiment
before inspecting these component contrasts. The whole-protein MACE/GB energy
failed all three predeclared comparisons; the matched decomposition identifies
both a local alpha failure and a PQQ reversal caused by the tested global term.
No direct-score threshold or parameter is changed to rescue that model.

## Concept and fixed score

Evaluate one new descriptor:

    R_short = interaction_energy(Ca) - interaction_energy(La).

Use the `interaction_energy` output already saved from both medium and large
whole-protein MACE calculations on all five development cases. Verify its actual
meaning against the installed model source before interpretation. Also report
its corresponding local-core value for the medium model on both declared
hydrogen preparations, and full-minus-core on matching global_H coordinates.
Retain other component contrasts only as an accounting check, not as a search
for fitted weights or a favorable combination. Do not reverse the sign rule.

Motivation: the learned local term describes whole-protein neighborhood chemistry
without artificial fragment caps. The subsequent global field/charge updates
showed strong checkpoint dependence and their tested total contribution harmed
ordering. A finite-neighborhood learned score might still be a useful empirical
discriminator. It is not a replacement solvation functional, a complete electronic
energy, or a binding free energy. Whole-protein input does not imply sensitivity
to arbitrarily distant atoms. No label-fitted coefficients, environmental scalar
or post hoc mixture is introduced.

Medium remains the primary candidate; large is a declared model-sensitivity
check. Larger R_short is predeclared as more La-like on this descriptor's own
scale. It inherits neither the baseline aquo offset/bands nor any absolute zero.
Same cases, groupings and three ordered differences as the completed global
panel:4MAE−1H4I,1F6S−GGR,6IP9−GGR. All are consumed development observations.
Preserve both alpha structures; do not select the favorable one.

This first stage reads saved outputs only: zero new MACE, GB or DFT evaluations.
Check identity/protocol/coordinate compatibility and exact component algebra.
If the descriptor looks useful, a separate declared implementation/validation
stage may extract it without evaluating the expensive long-range blocks, then
calibrate/evaluate it on a larger real panel. Passing these three comparisons
alone cannot complete the goal or establish broad predictive accuracy.

Whether promising or failed, retain all five values from both checkpoints and
all medium local values. If it fails, do not choose a component combination or
threshold after seeing these cases and present that as validation.
