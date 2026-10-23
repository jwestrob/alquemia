# Frozen-point water well-depth diagnostic

Root scientific decision, 21 October 2026. This continues the active, user-authorized
La/Dy exchange goal. It does not change a completed fit or release a new candidate.

## Evidence and question

Water continuation17593 completed14 additional iterations and1472 native programs,
stopping when a complete batch of four proposed steps gave no improvement. Actual
result SHA256 `fa2ef37ca71e79b27fdd8c27f8e7b694682ec2416e6ba477c7e64b77c9e79e09`;
saved review `c563c3904e8d175f9c97860fa7274c71e89507e0d4ed11750710933a5fb01003`.
All nine energy controls, all Dy and all paired rigid-response controls pass. The
three La dimer RMS errors remain2.374973,2.116768 and2.083147 kcal/mol/Å, above the
unchanged2 threshold. No maximum component exceeds5. The selected objective is
0.4194902374874799. The final measured four-column Jacobian is at this exact point;
rank4, condition62.262, with approximately2.1e-12 predicted local loss reduction.
This supports examining an additional interaction parameter, without claiming a
global optimum or proving that the four-parameter family cannot work elsewhere.

Question: do the two previously fixed metal–water pair well depths provide local
response directions that can address the remaining residuals beyond radius and
polarization damping?

## Frozen inputs and finite work

Use the existing eight consumed water states and the unchanged weighted residual
vector, force definitions and energy/force controls. Freeze the actual selected
radius/damping vector:
`[3.832755583060068, 0.31164088100092946, 3.525530001514476, 0.32135326991694113]`.
Keep all geometry, charges, multipoles, intrinsic polarizabilities, bonded terms,
water-H and non-water pairs unchanged.

Baseline explicit water-O pair epsilons are La0.24429636467571708 and
Dy0.37332330078584813 kcal/mol. Perturb each independently by±5%, giving four
parameter vectors. Each changes four physical cells of one metal; reuse the other
metal exactly. Each new cell uses the existing native analyze and analytic-force
programs: **16 new cells, at most32 native programs**, four single-thread workers.
Use Slurm standard_expand with default memory/time and the established engines.

Reuse the final eight radius/damping derivative evaluations only after verifying
the exact selected origin, original±0.005 normalized parameter steps, actual
clipped widths, inputs, native producer records and weighted residuals. If that
admission fails, stop this preparation for review; do not silently increase the
declared32-program inventory. There is no central rerun or Cartesian finite-
difference panel in this diagnostic. No selected new parameter vector is produced.

## Interpretation and stopping

Follow the existing PROPOSAL.md: compare the old four-column span with the two
epsilon columns at this same point. Preserve coordinate normalization, singular
values and the same absolute singular threshold for the orthogonal increment.
Report projected residuals by all four loss blocks and individual fixtures. Verify
epsilon linearity against the saved center at the existing energy/force replay
scales. A projection is a local prediction, not achieved improvement or physical
validation. Preserve every probe's actual score and control failures without
selecting a winner.

Stop new native starts on a native/nonfinite/polarization failure; retain completed
and running sibling evidence and the complete eight-cell denominator. Finish after
the four declared vectors. Any subsequent bounded six-parameter fit requires a
separate declared range, trust convention and finite budget within the active goal.
Donor transfer and fitting remain unbound while the water candidate fails. The
reserved24 outcomes remain unopened and production is unchanged.
