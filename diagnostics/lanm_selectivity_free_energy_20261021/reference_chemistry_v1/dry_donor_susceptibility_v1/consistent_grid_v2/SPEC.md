# Uniform matching-COSX donor response recovery

## Prospective question and finite inventory

Does making the late-SCF COSX grid equal to the unchanged final grid resolve the dry-donor energy–dipole inconsistency under the original numerical criteria? The original ten SPs remain preserved and unqualified (14/16 checks); both Dy odd-energy slopes failed. Saved diagnosis `99997dce` assigns 74.36%/91.80% of those discrepancies to the final integration contribution. This is a diagnostic decomposition, not a permitted subtraction or proof of cause.

Run exactly ten new SPs: the same two consumed central La/Dy geometries, each at zero and the four signed fields −0.001, −0.0005, +0.0005, +0.001 au along the original fixed tangent. **Dy's complete five-cell block runs first.** Its native/property/grid admission, all six original numerical checks, and an independent actual cost review must pass before the matching five La cells. No old/new recipe mixture enters the paired result. No old SP zero, EnGrad center, rotation, model call, extra field or retry enters this inventory.

## Only numerical recipe change

Relative to each original prepared SP input, insert exactly:

```text
%method
 AngularGridX 2,4,4
 IntAccX -1,-1,-1
end
```

The donor had no `%method` block. Do not copy the historical Mex `DoEQ false` setting. Preserve all other tokens, including MaxIter300, PModel, NoAutostart, PBE0-D4, RIJCOSX, DefGrid3, VeryTightSCF, the exact ECP/GTO/Aux-J assets and effective/physical state labels. Resource fields follow the grant as before. No gradient is requested. Keep the original eight real centers, charge +2, restricted effective singlet, 40 explicit electrons, 189 orbital functions and La580/Dy598 Aux-J functions.

Arendt's saved grid evidence (`c5a2067d`, report `7baced77`) pins the executed Mex syntax and all original donor grid fingerprints. The old Mex calculation attained matching grids but failed its derivative test. Its evidence supports syntax and activation only. The new native output must reproduce the original XC, grid1 and grid3 physical fingerprints and make grid2 exactly equal to original grid3, including angular/radial rules, IntAcc, pruning, weights, cutoffs, contraction, rotational invariance, H/He reduction, fitting and point count. Final grid3 reintegration remains present. Runtime timing/batch partitioning is excluded. La point counts are XC107266/grid1 16780/grid2=grid3 52712; Dy107343/16804/52765. These counts alone are insufficient for admission.

## Unchanged estimators, bounds and scientific limitations

The inherited original SPEC and SOURCE_INVENTORY remain pinned. The tangent is `(0.7562812564509392, -0.15594410634600103, -0.6353897204369424)`. All fields and property dipoles use laboratory origin zero. Preserve the +2 charged-origin treatment, effective nuclear-charge reconstruction, full vector/property tokens and electronic/final/D4 joins. Use property SCF electronic energy for curvature, retaining total and D4 separately; D4 remains constant across each fixed-geometry field block. No correction is subtracted from any energy.

For each step s, retain chi_E = −[E(+s)+E(−s)−2E(0)]/s², chi_mu = [mu(+s)−mu(−s)]/(2s), and mu_E = −[E(+s)−E(−s)]/(2s). Preserve all six per-metal checks: two energy-versus-moment curvature comparisons and both curvature step comparisons at 0.5 a0³; both odd-energy/zero-dipole comparisons at 1e−4 ea0. Preserve four paired curvature checks at 1 a0³. Ten slots and sixteen checks remain explicit; missing or failed slots cannot pass. All estimates and coefficient-propagated print bounds remain reported without averaging, extrapolation or threshold changes.

Native Energy Check and all false printed density/DIIS/orbital criteria remain visible. Matching integration does not establish SCF stationarity, a physical electronic state, PBE0 accuracy or a unique error mechanism. Failure still leaves numerical error and finite-step response unresolved. Old La admission does not exempt new La from the matching recipe. No paired response target or corrective coefficient is admitted from an incomplete or failed new block.

## Execution, reuse and costs

Root owns release/submission. Each SP uses four single-thread MPI ranks on one node in standard_expand, default memory/time, MaxCore=floor(0.75*actual MiB grant/4). Exactly one invocation per cell and ten maximum; native failure stops further claims. Preparation/collection use one local numerical thread and zero native calls. Completed historical files remain immutable. Reuse only geometry/assets, source parsers/estimators and the independently admitted model40 comparison; do not reuse old quantum response energies in the new ten-cell estimate.

The original ten field SPs used 143.6–180.4 s per call and 6320.94 allocated CPU-s proxy (1.756 CPU-h). The new denser late-SCF grid can cost more. A rough planning estimate is 2–4 CPU-h for ten; this is not a guaranteed bound. The actual complete Dy5 native/cost receipt controls continuation. Default scheduler grants are unchanged. The source/preparation review and actual Dy gate are required before the root releases their respective stages.
