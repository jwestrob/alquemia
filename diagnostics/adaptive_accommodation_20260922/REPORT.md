# Restricted donor searches leave a larger response untested

**The existing force archives justify testing a small adaptive donor-coordinate
set. They do not show that broader relaxation improves discrimination.** No new
energies, optimization, labels, thresholds or production changes were made.

## What the real forces say

All 30 original contexts and all 225 primary folds are retained. The original30
include 28 consumed references and two unlabeled PLM examples. Primary folds are
correlated structural repeats. Their 208 prepared pairs supply 415 paired-endpoint
origin/final force records: the original missing q0 endpoint remains missing,
and 17 unprepared folds remain in the denominator. Joint origin diagnostics are
available for 30/30 and 207/225 cases.

For each existing physical mode, the mapped Cartesian derivative includes actual
cap chain rules. Divide its generalized gradient by the Frobenius norm of that
mode's **physical heavy-atom** displacement Jacobian. The result is force per
angstrom of total heavy-coordinate displacement, not energy, stiffness or a
population. This makes rotations with different lever arms comparable without
mixing radians and angstroms or diluting a local motion by protein size.

| Final vacuum-OMOL diagnostic | Original30, 60 endpoints | Primary folds, 415 available endpoints |
|---|---:|---:|
| Median largest load in previously active modes | 0.0382 | 0.0337 |
| Median largest load in omitted angular modes | 12.7398 | 12.8280 |
| Omitted angular load exceeds active load | 60/60 | 415/415 |
| Median metal-translation vector norm | 28.2075 | 23.0687 |
| Archived final boundary flags | 0/60 | 8/415 |

Loads are kcal/mol/angstrom under the declared geometric norm. Final Ca and La
geometries differ: these are separate residuals, not a differential gradient at
one structure. At the **common origins**, the largest angular Ca-minus-La load
lies outside the former active set in **23/30** original cases and **180/207**
complete fold pairs. Glu chi1 supplies that largest load in 18/30 and 137/207;
Asn chi1/chi2, Asp chi1 and Glu chi2 also appear. No biological label entered
these calculations. Force size alone cannot distinguish useful compliance from
a stiff, highly constrained direction or Hamiltonian error.

## Reusable coordinates and proposed selector

PQQ source maps contain 10–12 modes: three metal translations and 7–9 original
donor torsions. Implemented residue templates support Glu/Gln chi1–chi3 and
Asp/Asn chi1–chi2. Actual PQQ records here contain Glu, Asp and Asn motions.
Existing non-PQQ mappings additionally support complete peptide crankshafts;
proline or missing amide neighbors are explicitly unsupported.

The standalone preview selects up to four common angular primitives. It alternates
the largest differential load and largest individual Ca/La load after projection
out of already selected directions. Modified Gram-Schmidt uses the same physical
heavy-coordinate metric; residual norm <=1e-5 is numerically redundant. Both
endpoint covectors undergo identical projections, with deterministic ID ties.
This avoids purely common-mode strain monopolizing selection while retaining
individual strain that cancels in the first derivative of the contrast.
No curvature, fitted weights, labels or score corrections are introduced.

Metal translations must be treated as one isotropic three-dimensional block,
not ranked as arbitrary lab-frame x/y/z axes. They are reported but excluded
from this initial four-angular-mode preview. Synthetic caps have no independent
degrees of freedom. Added scaffold/second-shell residues and PQQ remain frozen.
The existing water `WaterCoordinates` supplies six rigid coordinates per real
water, but is not integrated here: a heavy-only metric cannot measure water-H
orientation, and the PQQ inputs are dry. Do not invent waters or treat the failed
wider-water-basin entropy model as validated response.

## Solvent and prior negative evidence matter

The completed seven-context analytic response supplies a separate exact-source
comparison. Solvent changes the four-mode preview in **4/7** cases and the top
differential motion in **2/7** (4MAE and 1GLG). In 4MAE the top motion changes from
extra-Asp301 chi2 to Glu172 chi2. Fullfold archives have vacuum-OMOL forces;
their GFN2 singlepoints do not supply composite gradients. The selector is
therefore vacuum **proposal guidance**, not solvent-aware stationarity.

The previous complete-coordinate pilot already optimized these broad supported
motions, narrowed every tested native separation and ended at all ten 0.20-Angstrom
boundaries ([report](../coordination_preparation_20260919/REPORT.md)). That is a
negative utility result, not missing engineering. A small adaptive proposal plus
the parent's common geometry pool is a distinct test, whose benefit remains open.

**Recommendation:** after common-pool scoring, test the frozen four-angular-mode
rule as one contained challenger. Reuse `Kinematics.evaluate(q)` and its mapped
Jacobian; do not copy terminal-carboxylate +/-0.8-radian bounds onto proximal chi1/
chi2. Declare physical displacement limits before new energies. Preserve common
Ca/La active coordinates, individual candidate failures and the complete source
denominator. Do not add a force scalar to the score.

## Artifacts and checks

- Plan: [PLAN.md](PLAN.md); runnable replay: [COMMANDS.md](COMMANDS.md).
- `workspaces/adaptive_accommodation_20260922/force_projection_v1.json`:
  all raw/normalized mode loads, source pins, residuals, previews and failures.
- `workspaces/adaptive_accommodation_20260922/composite_comparison_v1.json`:
  all seven existing analytic comparisons, with separate Hamiltonians.
- Five real-fixture geometry/algebra tests pass in **0.898 s**, zero skips:
  actual derivative replay, rigid-axis invariance, redundant-coordinate removal,
  truncated-force rejection and common-measure composite replay.
- **New molecular calls: 0; new cluster allocations: 0.** No new native response
  qualification, predictive improvement, entropy or affinity correction claimed.
