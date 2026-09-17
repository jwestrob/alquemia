# Large-checkpoint comparison under standing pilot approval

2026-09-16. Jacob approved all contained discriminator pilots and superseded
separate launch approval; source ../mace_analytic_20260916/AGREEMENT.md.
The analytic medium pilot is complete and passes its core/full rotation and
charge checks. Its hybrid partition residual remains -5.023220584 kcal/mol.

Question: does MACE-POLAR-1 large improve the same consumed partition comparison,
and can it evaluate the same 9,141-atom protein affordably? Change only the
pretrained checkpoint, retaining the analytic kernel, float64, state/geometry,
assembly, vacuum/no-water policy, charge/spin and old DFT endpoints. Larger
network receptive range is an intrinsic checkpoint difference, not a tuned
physical cutoff. This is method development, not a blind biological benchmark.

Official checkpoint URL:
https://github.com/ACEsuit/mace-foundations/releases/download/mace_polar_1/MACE-POLAR-1-L.model
SHA256: 9f65f8dc6ddaff1d631e299cb531376a7da5e68d1bef04f34a2d5073d5ef114b.
Reuse the isolated pinned software installation; no package/source update.

First test analytic fields/derivatives, covariance and unchanged large weights
using the same four real archived densities and the large model's own widths
and normalization. Then eight complete core calls: same four 1H4I cores, primary
and the same 37-degree rotation. If core numerical gates pass, four full calls:
La/Ca primary and rotation. Twelve MACE calls, zero new DFT. No fitting or new
labels, preparation variants, solvent, geometry search or default promotion.

Keep thresholds: 0.01 kcal/mol rotation energy/contrast, 0.001 eV/Angstrom force,
1e-5 e charge. Kernel tests: absolute 1e-8 plus relative 1e-10. Report the
2 kcal/mol partition check separately; do not adjust it after observing results.
First use one A5000, 16 CPUs, 64,474 MiB requested host RAM, native mode.
Use 256-atom pair tiles, 2,048-edge and 128-node neural blocks (cores 128/17)
to conservatively accommodate the larger network without changing neighborhoods.
Same-model memory-only recovery may use the previously approved H200 allocation
or proportional shares if needed; no change of physical model to avoid an OOM.

Measured medium full cost is ~58 seconds per endpoint. Large runtime/memory
are unmeasured and reported from actual receipts; expected order is minutes.
No project CPU/time stopping budget. Preserve failed attempts and explicit
unavailable outcomes. Baseline stays available/default. New method version
and side-by-side medium/large comparisons; no inherited aquo reference or class.
