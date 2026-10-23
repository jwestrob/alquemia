# Joint W8/W9 development: independent design and source review

No design or incremental driver/objective blocker was found. This is conditional
on the final evaluator, input binding, and actual zero-call preparation review.
No native calls or new numerical controls were run.

- The residual count is **413 =8 +3×(3×9 +51 +57)**: eight energy terms and
  three force kinds across the three dimer geometries, W8, and W9. All **10 energy
  controls and 15 fixture force controls** remain individually required.
- Four blocks each have weight 0.25. Energy divides equally between four paired
  W8-referenced contrasts and four individual radial works. Two paired radial
  works remain acceptance checks only. Each force block averages five fixtures
  while retaining the original category/body projection. The new aggregate loss
  cannot be compared directly with the old weighted loss; the old eight-cell
  score is separately replayed and reported.
- The exact17605 selected six-vector is inside all original intervals. The new
  source requires zero-call replay of its eight cells plus the two17670 cells.
  Preparation must preserve the same uncorrected ECP-PBE0-D4 target, states and
  geometry; CC/CP diagnostic quantities must not enter these force targets.
- Twelve probes plus four trials per step give at most384 new vectors in24steps.
  Active-metal reuse gives at most120 derivative and80 trial programs per step,
  or4,800 total. The5,120-program hard ceiling is conservative. Eight single-thread
  native workers fit the eight-CPU allocation.
- The new driver retains exact complete-batch failure handling, probe construction,
  actual clipped spans, four-trial ordering, strict decrease/tie behavior and
  trial-only selection. Its entire outer-loop AST is unchanged except the explicit
  joint-pass stop label. Existing45saved controls therefore remain applicable;
  no duplicate solver experiment is needed. NumPy `rcond=1e-6` is a **relative**
  singular-value cutoff. A recorded Jacobian belongs to its origin and must not
  be presented as final-point rank after selecting a different trial.

The residuals share energy anchors and La/Dy force data. Their weights define an
engineering objective, not independent observations, a statistical likelihood,
or parameter uncertainty. A pass shows representability on ten consumed cells;
a local stop or bound contact does not prove global impossibility. W9 changes
the eight original water placements as well as water count. Previously exposed
orientation results cannot be reused as fresh blind validation.

The17670 pre-submit root-review artifact remains unavailable after it was
overwritten. Its explicit provenance notice must accompany the intact native
evidence; this review does not reconstruct the lost file. Any changed water
candidate must also remain separate from donor17621, whose water parameters are
frozen17605, until a separately declared matched donor evaluation.

## Evidence

The matching workspace `independent_review_v1/DESIGN_REVIEW.json`
(`a2e1b1f0e81a9b2e7ffdb1a2ba88f25e745aee85003e6825494892c9ab7403da`)
binds declaration897e239f and the prior source/evidence. `SOURCE_REVIEW_v1.json`
binds joint REFINE3c24c3c4, OBJECTIVE60274b54, unchanged direction source741167a8,
and original saved-control receiptd8e9aa02. Final evaluator/config/preparation
review remains pending; these source checks do not release execution.
