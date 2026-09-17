# Exact derivatives of the masked descriptor: engineering and sensitivity

Declared after the GGR robustness and multisite results, before any derivative
implementation or new derivative/displacement evaluation. Active goal and
blanket authorization apply. This is a new development experiment, not a rescue
of the failed affinity comparisons. Baseline and masked scalar definitions stay
unchanged; no relaxation correction, fitted stiffness or entropy is authorized
by a numerical gradient pass alone.

## Question and method

Can we obtain affordable, exact coordinate derivatives of the existing masked
whole-chain descriptor, then use them to identify which physical displacements
change its score? Energy-only inference now works. The qualified edge/product
adapters explicitly reject gradients; do not merely remove that guard and claim
support. Add a separately versioned backward-capable adapter with activation
checkpointing/recomputation where needed, preserving the same learned weights,
charge-feature mask, spin feature, graph, native algebra and float64 arithmetic.
Use analytic automatic differentiation, not numerical DFT derivatives.

Report grad(R_mask) = [grad(T_Ca)-grad(T_La)] times the established conversion;
the fixed disconnected constants have zero coordinate derivative. Keep the
negative-gradient convention distinct from the gradient. These derivatives are
of a learned masked descriptor, not validated nuclear forces or gradients of
the ORCA/CPCM baseline. No force-field or environmental gradient is added.

## Finite inventory, conditional on numerical checks

Stage A uses the exact real 73-atom 1H4I core coordinates from
`charge_ablation_development_v2/manifest.json`, tasks
`bridge_1H4I_{Ca,La}_mask_native`. Ten forwards total:

- Two native, unbatched analytic-gradient centers.
- Two corresponding batched/checkpointed analytic-gradient centers.
- Two rigidly rotated batched gradient centers using the existing fixed
  `mace_hybrid.rotation()` transformation.
- Four batched energy-only signed displacements: each metal at plus/minus
  0.01 A along the center metal-to-nearest-oxygen direction. Choose the oxygen
  from the fixed real center coordinates before any displaced output; ties
  resolve by source index. Only the metal moves.

Stage B runs only after Stage A passes. Six forwards on the existing intact
GGR1GLG preparation: two batched analytic-gradient centers and four signed
energy-only displacements of the metal by plus/minus0.01 A toward source
GLN140/O. Protein, protonation, water inventory and all other atoms stay fixed.
Compare center energies to actual archived masked bound endpoints. Export all
source-mapped endpoint gradients and grad(R_mask), plus its projection onto
this physical metal direction. Do not label an unexecuted stage successful.

Sixteen new model forwards if both stages run; no DFT, training, Hessian,
optimization, trajectory or new biological comparison. Existing saved scalar
results cannot supply these derivatives. Record all attempts and costs.

## Acceptance and interpretation

Native/batched and rotated energies must agree within0.01modelkcal; Cartesian
gradient comparisons within0.001eV-equivalent/A after the proper rotation.
The actual odd signed energy change must match h times the analytic projected
gradient within max(0.01modelkcal,1% of the absolute predicted odd change), for
each endpoint and the paired score. Those tolerances are fixed before outputs.
Keep negative curvature/large gradients if encountered; no clipping or spring
fitting. Finite displacements check derivatives, not a relaxation energy model.

Use the current A5000 allocation convention16CPU/64474MiB, pinned MACE software,
existing runner and isolated versioned implementation snapshots. First measure
backward memory/time on the real core. Whole-chain memory feasibility is open;
retain any failure before choosing a documented technical recovery. These are
development checks; routine scoring retains its existing two energy-only calls.

Next scientific decision after this experiment: whether a bounded, physically
defined structural-response score is defensible and useful. Gradient availability
alone does not validate curvature, metal-adapted geometries or an affinity score.
