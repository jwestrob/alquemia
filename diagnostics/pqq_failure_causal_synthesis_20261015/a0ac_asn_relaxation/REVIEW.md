# Independent review of rebased-start preparation

No material static blocker found in the first complete preparation:

- `PREPARE.py`: `a60ca8a0aeba2ca51d2347f9298512e3a1f31a0bbef65493589bbc8e31a021c3`
- `PREPARE.sbatch`: `616add77e0858b96798b99f2931057a3f6c01b0c248a393bf007d95f20ce8d4f`
- `PLAN.md`: `ce8b3665d02ea2dc901f2201b6d5dc6def1dfed4085dbf3be203f2042d634963`

The code grafts only the four mapped amide sites into the physical origin and
sets the204 scored coordinates to the actual admitted flip. It preserves source
identities, physical exterior, topology, modes and cap-link definitions. q0 must
reproduce the saved flipped coordinates and unchanged caps. The existing analytic
Jacobian/covalent checks run on Slurm; cached flip forces are projected through
the new map and independently contracted with the scored-coordinate Jacobian.
Old-origin gradients are not reused. Unsupported selected motion outside mapped
physical sites or across cap anchors fails closed.

Actual source/task metadata supports the accessed schema. The original paired
four mode IDs and numerical settings remain fixed; La-conditioned0/2 do not gain
new modes. The separate0.8Å/±0.8rad domain is explicitly relative to the admitted
flipped start. No old-domain admission is retroactively changed.

The contact union keeps the local204 sites first and the exact original omitted
source identities afterward, including hydrogens. It carries the saved moving
sites and cap identities forward for subsequent physical-contact checks. All ten
source outcomes remain visible, with failed preparation explicit. No new model
call occurs in preparation; all twenty actual flip q0 force/energy records are
reused with the new projections after the original sixty baseline checks.

This static clearance supports the finite one-CPU Slurm preparation only. Actual
q0/Jacobian/geometry gates remain unexecuted at review time, and the separate
optimizer/scoring runner still requires review. The first-source paired scout
and release of the remaining18 searches remain root-owned. No scientific imports,
tests, numerical calculations or submissions were performed by this reviewer.

## Actual preparation and narrow recovery

Original4814 passed5/10 sources. Saved-map diagnostic4819 established identical
old/new computed caps and exact raw inputs/anchors in all ten; only the raw-versus-
computed representation comparison failed at existing1e−16Å roundoff. The v2
delta was independently reviewed by Gauss: exact raw/raw and computed/computed
comparisons, unchanged1e−12Å full q0 checks, pinned original failure/diagnostic.
Actual4821 then passed10/10 sources and20 prepared tasks, with zero energy calls.
See [the complete preparation report](PREPARATION_REPORT.md). This closes the
preparation issue; it does not replace the separate optimizer-runner review.
