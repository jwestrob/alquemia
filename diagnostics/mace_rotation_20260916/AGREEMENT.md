# Approved rotation-attribution investigation

2026-09-16. Following the completed memory-work report and recommendation to
investigate rotation sensitivity, Jacob: **“Proceed and investigate!”**

Question: does the observed orientation dependence come from the memory rewrite
or the unchanged realspace dipole representation? Preserve the same four consumed
1H4I qm33/qm36 La/Ca cores, medium checkpoint, float64, charge/multiplicity,
protonation, no-water/vacuum state and numerical parameters. Reuse archived
unrotated results; rotate by the previously declared 37 degrees about [1,2,3]
and the source metal. Eight new core energy/analytic-force calls: four with the
original interface-repaired backend, four with blocked execution. No new DFT,
full-protein call, large checkpoint, training or biological case.

Also inspect each backend's field and energy modules with the four archived
computed densities. Compare original coordinates/dipoles with the same physical
rotation. Separately co-rotate the auxiliary displaced-charge axes, expressing
dipole components in that co-rotating frame, to test whether fixed laboratory
axes are the source of the finite-stencil error. This is a frozen-distribution
representation diagnostic, not an alternate model score. Retain exact offsets,
widths, self terms and exclusions. No offset sweep or replacement kernel here.

Use existing runner, immutable manifests, resource receipts and original frozen
physical tolerances (energy/contrast 0.01 kcal/mol; force 0.001 eV/Angstrom;
charge 1e-5 e). Original-vs-blocked tolerances remain energy/force 1e-6 in eV
and eV/Angstrom, density 1e-8. For co-rotated stencil identity use absolute
1e-8 plus relative 1e-10, matching prior kernel tests; report exact residuals.
No tuning from classifications, which are unavailable for this vacuum experiment.

Expected allocation: one A5000, 16 CPU cores and 64,474 MiB requested host RAM;
minutes based on observed original core costs. No project time/compute stopping
budget; cluster seven-day QOS limit remains. Scientific variants beyond this
attribution investigation require discussion; routine implementation fixes and
retries preserve the agreed inputs and methods. Baseline unchanged.
