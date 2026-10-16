# Mex local field: four converged endpoints and analytic gradients

The 209-atom Mex EF3 region converged for both La and Dy, with and without the
actual 1,840-site protein/solvent field. All four numerical and state gates pass.
This supplies an affordable local environmental-response component. It does
not yet qualify a combined protein energy or establish metal preference.

| Metal | Isolated SCC iterations | Embedded SCC iterations | Embedded minus isolated (kcal/mol) |
|---|---:|---:|---:|
| La | 48 | 74 | -118.055578090860 |
| Dy | 50 | 77 | -118.096897409121 |

The field stabilizations are similar at this source geometry. They are neither
binding energies nor evidence for a reliable La/Dy classification. The existing
original spin-independent, f-in-core GFN2 approximation remains explicit; it
has not acquired open-shell accuracy by converging. The physical local La
singlet and Dy sextet are separate metadata.

Both native QM and point-charge gradients were collected and mapped through
both cap anchors to the 2,049 physical sites. Translation-gradient sums are
recorded in RESULT.json; external gradients have eight-decimal serialization.
Independent directional finite differences are still required. The field uses
the actual source charge distribution and explicit element-specific GFN2
hardness, with no point-charge self-energy or second QM–MM Coulomb term.

Preparation job 12569 used one CPU for one second and no molecular calls.
Job 12570 used four CPUs for 18 seconds: 72 allocated CPU-seconds. Its application
ran 17.219 seconds; observed child CPU time was 56.626 seconds and peak child
RSS was 93,840 KiB. Scheduler batch MaxRSS was 371,688K. No GPU, MACE, solvent,
DFT or optimization was involved in these four calls.

Next: assemble two static whole-system states from local MACE, this field
correction, the already available classical ledger, and one explicitly defined
fixed-charge global solvent term. Reuse accepted components. Check the combined
energy derivative, including boundary and environment motion, before relaxation.
Exact commands, source/state/runtime pins and unabridged gradients are in the
workspace files linked by RESULT.json.
