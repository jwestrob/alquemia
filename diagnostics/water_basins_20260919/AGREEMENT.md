# Approved finite water basin study — 2026-09-19

Jacob: “let's get this and the water occupancy/entropy process each going in parallel.”
Root reviewed and approved this concrete scope before any new scientific execution.

Use the four consumed 1F6S11/6IP9110 Ca/La latest actual native DFT centers.
Hold scaffold, metal, outer waters, internal water shapes and inventory fixed.
At298.15K, compute the finite12-dimensional rigid-water integral of the existing
Cartesian-DFT-gradient-anchored native OMOL potential. Two independent scrambled
Sobol sets of1024Gaussian importance draws per endpoint; full coupled curvature
provides the proposal covariance, inflated by1.5 in standard deviation. Gaussian
normalizer and rejected draws remain in the denominator. Only samples in the
largest domain require new energies; maximum8192sample energies plus4anchor checks.

Domains per water COM radius0.30/0.45/0.60Angstrom and SO(3) rotation radius
0.50/0.80/1.10radian are fixed before outputs. Both metals share the same physical
midpoint COM/orientation domain, after mapping equivalent water hydrogens to the
nearest physical orientation. Actual endpoint atom ordering and shape are retained;
preexisting sub1e-9Angstrom shape serialization differences are recorded. Distinct
water COM spheres must not overlap. Rotation radius belowpi/2 avoids duplicate
H-exchange orientations. No factorial or arbitrary symmetry factor is inserted.

The measure is product d3COM times SO(3) Haar volume
sinc(theta/2)^2 d3rotation, explicitly divided by1Angstrom^3 per water for a
reported dimensionless configurational integral. Its absolute additive gauge is
not entropy; common kinetic/measure factors cancel only for matched water counts.
No occupancy probability or absolute entropy is released. Bound internal-water,
quantum and non-electrostatic solvent terms remain missing, not zero.

Necessary numerical gates, frozen before sampling: independent scramble F
agreement<=0.10kcal/mol and effective sample size>=128per endpoint/domain.
Report half/full refinement and outer10percentdomainweight. Passing these does
not establish global basin coverage. Failed domains remain explicit.

Select exactly4distinct native validation geometries per endpoint, using only
sampled geometry and cheap weights: weighted-coordinate medoid, weighted50th and
90thnormalized radial extent, and highest-importance-weight point in outer10percent
domain. Native r2SCAN-3c/NoAutostart/CPCM(Water)/DefGrid3/TightSCF EnGrad,
16new endpoints total. Necessary sampled-target checks: representative-weighted
absolute DeltaE error<=0.25kcal/mol and maximum error<=0.50kcal/mol. Selection is
not an unbiased DFT-reweighting sample and cannot qualify fullDFT thermodynamics.
No threshold tuning, new recentering, training or fullDFT Hessian.

Existing runners, isolated frozen dispatch and existing environments. OneH200,
16CPUs,200000MiB for sampling; native16MPI ranks per endpoint/four concurrent.
Estimated tens of GPU-minutes and20–35minutes native node time from prior receipts;
actual costs measured. Finite manifests specify scientific work; no arbitrary
compute/time-budget stopping rule. Baseline/default/PQQ remain unchanged.
