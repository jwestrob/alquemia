# Next modeling step after the arrangement table

The approved direction is local bound-water motion with MACE proposals checked
against native DFT. The arrangement experiment itself keeps all water oxygens
fixed. It does not finish the free-energy occupancy model.

Existing saved analytic gradients already permit physical translation and rotation
projections for every retained water. The full-water configurations show directional
agreement for whole-water translations (cosines0.935–0.994), but nonzero DFT forces;
this supports testing local motion, not assigning a harmonic entropy immediately.
Use the final MOTION.json for the complete state inventory and actual errors.

The complete occupancy gradient table identifies a particularly weak direction:
6IP9La pattern110 waterA322 has a translation-gradient cosine of0.0617. Therefore
a concrete next check should use the paired full1F6S endpoints and the paired
6IP9 pattern110 endpoints (four centers, chosen by physical model disagreement,
not biological score). Translate each variable H2O rigidly along its own metal–O radial direction by
+0.05 and−0.05A, with protein, outer waters, water orientations and internal
geometry fixed. That is8 native EnGrad checks plus matched cheap MACE calls.
Check smaller MACE displacements at±0.025A for inexpensive curvature stability.
All moves retain water identity/count and the native r2SCAN-3c/CPCM policy.
Do not interpret collective radial motion as the complete6n-dimensional basin.

Use DFT energy and gradient at the center and test a MACE curvature approximation:

E_pred(u) = E_DFT(0) + g_DFT(0).u
            + E_MACE(u) - E_MACE(0) - g_MACE(0).u.

This includes the target CPCM force at the center and approximates the change in
force with vacuum MACE; it is not a self-consistent environment-corrected potential.
An initially defensible check reuses the earlier local-curvature tolerances:
anchored energy error at most max(0.02kcal/mol,25% of the target even energy).
Freeze exact displacement inventory and acceptance rules in its own manifest
before evaluating these new points. This file is a plan, not an execution receipt.

Two current6IP9Ca occupancy patterns also have distinct converged MACE orientation
minima (~2kcal/mol apart). Their alternative orientations need DFT checks before
calling the selected configurations DFT minima or using them as complete basin
ensembles. These are recorded alternatives, not reasons to choose a favorable
classification or to overwrite the24-endpoint arrangement table.

If inexpensive curvature passes physical checks, extend to coupled water
translations/rotations (up to18 physical coordinates here), with correct masses,
moments of inertia, basin multiplicity and coordinate measure. Retain unstable
modes or trust-region violations as failures; do not clamp them into a favorable
entropy. Internal bound-water vibration and non-electrostatic solvent terms must
also be reconciled with the already computed liquid-water reference. No numerical
DFT gradients or whole-system Hessian are needed for the initial local checks.

Status: proposed follow-on; none of these new displacements has been executed
as part of the20-endpoint arrangement experiment. The broader direction is
approved; exact follow-on execution must preserve these physical distinctions.
