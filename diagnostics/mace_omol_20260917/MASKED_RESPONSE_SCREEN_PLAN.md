# Next experiment: does the masked descriptor follow real DFT response?

Declared after the exact derivative qualification passed, before any new
masked-OMOL evaluation of these response geometries. Jacob's active-goal
blanket authorization applies. This is method development on consumed cases.

## Decision this experiment will support

Can the new masked OMOL model supply useful local response for the discriminator,
or do we need to retain DFT first-derivative anchors (or reject this response
route)? The whole-GGR derivative is now affordable and numerically correct, but
that does not make the learned masked potential a validated mechanical model.
Test its actual directional forces and deformation energies against the saved
DFT calculations before performing optimizations or adding score corrections.

Earlier POLAR+GB curvature checks passed but paired relaxation failed because
La was unstable or outside the declared trust region. Do not repeat those
calculations or present this as the first curvature test. This experiment tests
the distinct OMOL masked model; it has no charge output or GB term. Agreement
with solvated DFT, if observed, is an empirical approximation, not a claim that
vacuum OMOL and CPCM DFT have identical Hamiltonians.

## Fixed inputs and runs

Use the exact archived `workspaces/mace_mechanics_20260916/prepared_v2/manifest.json`
(SHA256 361d2ac85b93360c8560a761e0d76d622b8f314a7c238013f8fb7faa184ae99c),
including source graphs, physical displacement definitions and core Jacobians.
All four representations: GGR_extended, GGR_connected, ALPHA_1F6S, ALPHA_6IP9.
GGR representations are one structure; both alpha structures are one qualified
biological comparison. No fresh blind cases or new labels.

For each representation and each metal, evaluate exactly these five existing
points: center, s_minus, s_plus, theta_minus, theta_plus. The s coordinate is
metal translation of +/-0.02 A along its source-defined direction; theta is
rotation of the complete selected peptide unit by +/-1 degree. Read actual
coordinates, cap handling and mappings from the archived preparation; verify
them rather than reconstructing from this prose. Forty new OMOL calls total:
eight analytic-gradient centers and thirty-two energy-only displacements.
Keep original H positions, charges, water inventory, caps and all coordinates
byte-identical to the DFT inputs. Known original-H preparation limitations remain
explicit. No whole-chain subtraction with the differently H-normalized current
GGR preparation is permitted in this screen.

Reuse DFT energies and analytic gradients from the existing sensitivity and
mechanics receipts. No new DFT, GB/PB, training, Hessian, trajectory, optimization,
relaxation energy, entropy or biological classification is calculated.

Use the pinned OMOL checkpoint, mask and qualified autograd adapter v3 in
float64. Source spin multiplicity is preserved. Use the existing manifest,
snapshot and allocation runner: one A5000, 16 CPUs, 64474 MiB host RAM. Measured
core inference is about0.5–0.75 seconds per task; process startup dominates the
roughly several-minute projected allocation. Record actual failed attempts and
costs. Do not change the production scorer or its reference.

## Comparisons and criteria, frozen before outputs

Project the eight analytic center gradients through the actual physical
Jacobians, including link-atom chain-rule terms. Keep gradient versus force and
Hartree/Bohr versus eV/A conventions explicit. Compare each model's analytic
projection with its own odd signed energy difference; tolerance is
max(0.01 model kcal, 1% of the absolute predicted odd change), as in the passed
numerical derivative qualification. Report any curved-path even contribution;
these are directional secant curvatures, not Cartesian Hessian eigenvalues.

For each endpoint and paired R=Ca-La, report:

- Each signed deformation energy relative to its own center and its DFT error.
- Odd (first-response) and even (curvature) terms against DFT.
- The DFT-anchored prediction: DFT center/gradient plus the learned even term.
- Both GGR representations, retaining their difference rather than selecting one.

Reuse the earlier response-screen curvature and anchored-error criteria:
even error <=max(0.005 kcal/mol,25% of absolute DFT even term); anchored signed
error <=max(0.02 kcal/mol,25% of absolute DFT even term). Add a separate direct
response criterion: each unanchored signed deformation-energy error
<=max(0.02 kcal/mol,25% of the absolute DFT deformation energy). These are
screening tolerances at the declared small amplitudes, not universal error bars.
Every endpoint, direction and representation remains in the denominator.
No fitted scale, springs, eigenvalue clipping or threshold calibration.

If direct responses fail but DFT-anchored responses pass, pursue a matched
hybrid response with DFT anchors; do not claim unanchored minimization validated.
If curvature fails, this cheap curvature is unsupported. Even an all-pass screen
still requires coupled-coordinate and independent displacement validation before
any relaxation correction. Negative curvature and out-of-trust optima remain
failures; the old trust region is not enlarged. Prediction accuracy and complete
binding free energies are not conclusions of this experiment.
