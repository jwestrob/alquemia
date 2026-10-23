# Zero-field model redistribution counterpart

## Frozen question and scope

Does the unchanged failed-science 17676 AMOEBA water model reproduce the charge-neutral C/P/B molecular dipole redistribution of the consumed, fixed W8/W9 La/Dy geometries? Use the admitted selected six parameters and original SI water, all molecular terms unchanged. No EXCHPOL, external field, fitting, dynamics, Cartesian derivative check, reserved geometry or target change.

The finite inventory is 78 model calls: four complete clusters, 34 literal metal–water dimers, 34 corresponding free waters, four intact water shells, and two bare ions. Each ion is computed once and reused in both coordination cycles. Every fragment is an exact atom-mask extraction in its parent coordinate frame, with retained coordinate strings, types, masses, water connectivity and multipole-axis atoms. Previous four cluster and 34 dimer raw energies are admitted and replayed within 1e−4 kcal/mol. The two bare-ion energies and induced dipoles must be zero within the source/printing tolerances; they are not missing-data substitutions.

One native `analyze system.xyz EM` per cell, with reporting-only `polar-print` and `save-uinduce` appended to the exact selected key. One thread per call, at most four workers in one four-CPU Slurm task, default memory/time. Hard ceiling 78 program starts, zero testgrad/internal-FD evaluations. Preparation and tests perform zero molecular calls. Root binds a reviewed preparation and authorizes the run through a separate release; no automatic submission.

## Admission and numerical precision

Require successful native receipt, native version, finite total/components and printed component closure ≤1e−6 kcal/mol. Admit every printed induction summary only if finite and nonnegative, ≤1e−8 Debye; exact zero is valid for isolated fragments. Require precisely two induction summaries for EM (energy analysis and moments), positive finite intrinsic polarizability at every expected site, one complete n-site induced file, exactly ordered indices/names and three finite components per site. Any native/parser/admission failure sets global stop before further starts; retain all 78 slots and drain already running siblings.

Native induced files contain Debye to six decimals. Per component, a sum over n sites has print-rounding bound n×0.5e−6 Debye. For a cycle with coefficients c_s, propagate `Σ|c_s| n_s×0.5e−6`; this is a printing bound, not a rigorous induction-error bound. Use compensated `math.fsum` and keep the unchanged 1e−8 native residual separately. Native total M components have three decimals: propagate `Σ|c_s|×0.0005 Debye`. Convert each charged fragment's printed center-of-mass moment to the common laboratory origin with native masses and `D Q COM`, D=4.803206802 Debye/(e Å). Allow a further 2e−9 Debye per weighted row for this deterministic conversion. Require the corrected M-cycle and induced-cycle components to agree within the sum of these predeclared bounds. Retain both values and exact differences; do not use rounded M as the primary observable.

## Cycle and cancellation proof

For each metal and N=8 or9, use shared parent-origin dipoles:

```
C = full − ion − shell
P = Σ_i(dimer_i − ion − water_i)
B = full − Σ_i dimer_i + (N−1)ion − shell + Σ_i water_i = C−P
```

All are vector quantities in Debye. Retain each parent result and the `(Dy9−La9)−(Dy8−La8)` contrast separately. The permanent contribution cancels in C and every term of P and thus B: each parent atom has identical type, permanent multipole record, bond neighbors, and local-frame atoms/coordinates in every occurrence; the coefficient sum for each atom is zero. Water O uses the retained two-H bisector, H uses its O/other-H Z-then-X frame, and each metal has charge+3 with no permanent dipole/axis atoms. Preserve explicit signatures and atomwise coefficient counts in the manifest; do not infer cancellation merely from net charge. All rows keep the exact inherited key/prm. Source pins for kmpole/rotpole/moments and native atom masses support this proof.

Accordingly use induced sums for primary C/P/B; full low-precision M after origin correction is an independent consistency check. A total vector can be symmetry-suppressed, especially W8. This measures zero-field redistribution, not polarizability, atomic density, a unique charge-transfer/polarization partition, force accuracy, or a model rescue. QM comparison remains a separate saved-data integration; original energy failures and targets remain unchanged.
