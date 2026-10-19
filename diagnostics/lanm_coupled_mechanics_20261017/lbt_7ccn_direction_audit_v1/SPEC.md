# 7CCN-La conditional next-direction audit

## Frozen question

Does the direction the unchanged optimizer would construct from the exact 80-step 7CCN-La terminal checkpoint, **if continuation were permitted**, have a downhill local derivative under the admitted whole-PBE0 plus the same fixed-FF ALPB gradient? The original optimizer actually stops at its80-step cap before constructing this direction. This audit does not change that guard or authorize continuation.

Bind collector15173 (`fff95c2f`), its after-La force report (`e6b1cb93`), search12968, checkpoint (`4e3a6123`) and terminal POINT (`cf59ef7a`). Full256 physical IDs, coordinates, q−2/La singlet state and common ALPB receipt must match exactly. Reuse the admitted reference arrays, native receipt and pinned original component arrays. Preserve the two original solvent-conversion orders separately: model uses raw*(HA/BOHR); reference uses (raw*HA)/BOHR. No numerical tolerance is added.

## Direction and outputs

Call the exact pinned optimizer rigid projection and original L-BFGS helper with all eight archived history pairs. Reproduce optimizer lines216–224: project model gradient, construct and project L-BFGS direction, test sum(pg*d)>=−1e−14, apply the original −pg fallback and local history clear if triggered, test maximum atom norm, and scale by min(1,trust/maximum). Trust remains the saved0.1Å. Preserve raw/projected/pre/post-fallback vectors and branch values. No alternative direction is selected and no trial coordinates are formed or evaluated.

For the single resulting direction d, report dE/dalpha for X+alpha*d (kcal/mol, dimensionless alpha), and dE/ds along d/||d|| (kcal/mol/Å). Use both full saved gradients; reference is not rigid-projected for its derivative. Signs are descriptive, without new accuracy/admission thresholds. Preserve full256 model/reference forces, signed model-minus-reference errors, atom and residue max/RMS/L2 summaries and projected contributions. Display all atoms/residues; rank by error norm/L2 with deterministic index/residue-ID ties. These are force summaries, not residue energy decomposition.

## Admission, limits and interpretation

All numerical imports and execution require Slurm; one CPU, original MACE Python environment, BLAS threads1, no GPU. Fresh output only. Exact SHA/source/state/geometry/component/force-table joins precede direction construction. Any failure produces unavailable with partial evidence; no retry. Zero DFT/MACE/GFN/ALPB/MM calls, model loads, geometry trials and optimizer steps. No production edits or optimization extension.

Native reference admission is inherited, not a new derivative qualification. Preserve all printed SCF residual flags, vacuum-dianion, terminal-chemistry and frozen-f limitations. A negative local slope does not establish finite-step energy lowering, force accuracy or a stable minimum. Other reference cells are outside this one-direction audit; missing Dy model forces remain unavailable. Root reviews the actual result before any separate mechanical decision.
