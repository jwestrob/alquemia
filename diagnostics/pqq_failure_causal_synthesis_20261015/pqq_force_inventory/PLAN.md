# Saved PQQ residual forces — frozen before numerical loading

Root selected this finite analysis on15 October2026. Use precisely the eight pinned MACE vacuum-force arrays in INVENTORY.json: ExaA samples0/3 × origin/adaptive_La × Ca/La. DFT nuclear-force fields remain unavailable for all eight; their completed single-point energy comparisons are separate evidence. No scientific force/coordinate array has been numerically loaded for this analysis at declaration.

## Directions

The intact PQQ is27 physical atoms:24 heavy atoms and its three modeled H; no synthetic cap belongs to a direction. The named quinone pairs are generic C5–O1 and C6–O2, corresponding to CCD C5–O5 and C4–O4. For each stretch keep C and every other atom fixed, displacing only terminal O along the current C→O unit vector. Positive generalized force favors elongation. Do not use the earlier proposed symmetric C/O direction.

Use all27 atoms for rigid translation and rotation. The center is the centroid of24 heavy atoms. The right-handed chemical frame is e1=C5→O1, e2=C5→C6 with its e1 component removed and normalized, e3=e1×e2. Report the three raw total-force components and three raw torque components in this frame; torque is about the same centroid. Retain global vectors and raw direction norms.

Normalize each direction by the square root of the sum of squared heavy-atom displacements, giving norm1. Translation therefore has normalization sqrt(24); the rotations retain their actual heavy-atom lever-arm norm. Include H displacements and forces in each rigid projection despite excluding H from its normalization. Report all normalized projections in eV/Å of collective heavy-atom displacement. Keep raw torque in eV/rad, not a force in eV/Å.

## Checks and comparisons

Before projection verify all pins, exact request/XYZ/state/model identities, finite192×3 forces, the negative total vacuum-energy-gradient convention, the complete name/element/source map and current donor-mode exclusion of PQQ. Require paired Ca/La coordinates identical with only the intended first element differing, and fixed PQQ coordinates across origin/adaptive_La within each source.

On the actual saved arrays, chemical-frame orthonormality and heavy-direction normalization residuals must be below1e−10. Jointly rotate every coordinate and force by+90° about z and translate coordinates by(11,−7,3)Å; recompute the chemically defined frame/directions and require every normalized projection to remain invariant within1e−9eV/Å. The transform is an algebraic self-check, not a new model evaluation. These tolerances are fixed by root before projection; do not alter them after inspection.

Preserve per-cell loads, Ca−La differences and arithmetic Ca/La mean loads at each source/geometry. Report proposal−origin changes for each metal and the corresponding differential/common responses, plus source3−source0 comparisons. Arithmetic common load is not a strain energy. Do not select axes or subsets using their resulting force magnitudes.

## Execution and limits

One shared Slurm CPU, with numerical libraries restricted to one thread. Use the existing pinned NumPy environment from the completed saved-force analysis4134. Zero MACE/SCC/DFT/optimization calls, no new coordinates for scoring, no class bands, no predicted energy correction and no production change. Root alone submits.

These are residual forces of the proposal Hamiltonian. They can reveal excluded PQQ degrees of freedom carrying load, but cannot uniquely establish bad structure, a redox mismatch, beneficial relaxation, a composite-solvent gradient or classifier improvement. Force-based causality and DFT validation remain unestablished.
