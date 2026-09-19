# Wider water basins: sampling works, single-anchor energies do not

**The existing DFT-gradient-anchored MACE potential does not reproduce native DFT
accurately enough across thermal-size water motions to supply an occupancy or
entropy correction.** All 16 native calculations completed, but none of the four
endpoints passes the predeclared broader-domain energy checks. This limitation
appears in both tested alpha-lactalbumin structures (1F6S and 6IP9). It does not erase the previously measured
success of MACE-assisted local water relaxation/preparation.

The useful distinction is now experimental: we can integrate a finite physical
water domain cheaply, but extending the local energy model across that domain
introduces consequential errors. More samples of the same potential would not
resolve its native DFT disagreement. The scanner/default/PQQ pipeline is unchanged.

## What ran

Four consumed development states:1F6S11Ca/La and6IP9110Ca/La, each with two mobile
rigid waters.6IP9110 is a selected two-water arrangement, not the full-water
6IP9 alpha-lactalbumin preparation. Both structures belong to one protein group;
this study contains no GGR endpoint and no new biological discrimination comparison.

We reused their latest actual DFT energies/analytic gradients and matching native
OMOL0-100M energies/forces. The approximation remains

`V(x)-E_DFT(anchor) = E_MACE(x)-E_MACE(anchor) + (g_DFT-g_MACE)·(x-x_anchor)`.

Protein, metal, outer waters, internal water shapes, inventory, charge and native
r2SCAN-3c/CPCM(Water)/DefGrid3/TightSCF method were fixed. The four independent
12-dimensional integrals used two scrambled Sobol sets of1024 Gaussian importance
draws each at298.15 K, with full coupled cheap curvature as the proposal only.
We evaluated6250 in-domain configurations plus 4 anchor reproducibility calls;
1942 rejected candidates stayed in the 8192-draw denominator. No DFT Hessian,
trajectory, reoptimization, spring fit or numericalDFT gradient was used.

Ca/La share midpoint-centered physical domains, with nested per-water COM radii
0.30/0.45/0.60 Angstrom and rotation radii0.50/0.80/1.10 rad. The integration retains
the exact SO(3) Haar Jacobian. Both1F6Swaters needed equivalent-H label mapping:
nominal3.118/3.134rad endpoint differences became physical0.0748/0.0607rad differences.
Actual endpoint ordering and coordinates were retained. Distinct-water domains
are disjoint and each orientation chart avoids H-exchange duplication; no arbitrary
factorial/symmetry correction was inserted. Largest inherited paired-shape
serialization difference is4.26e-11Angstrom.

For each endpoint we selected four distinct geometries before seeing native
results: weighted-coordinate medoid, weighted50th/90thextent, and highest-weight
outer-domain point. These16 native analytic-gradient calculations ran through the
existing ORCA manifest runner, overlapping remaining GPU sampling.

## Actual results

| Endpoint | Representative-weighted absolute error | Maximum error | Native energy gate |
|---|---:|---:|---|
|1F6S11Ca|0.56701|0.71975|fail|
|1F6S11La|0.33096|0.73101|fail|
|6IP9110Ca|0.42882|0.75950|fail|
|6IP9110La|0.68470|2.15629|fail|

Errors arekcal/mol. The frozen necessary limits were0.25 weighted and0.50 maximum.
For a concrete consequence, the 6IP9 La outer representative costs 0.90324 kcal/mol
in actual DFT, whereas the anchored MACE model predicts 3.05953. At298.15 K this
underestimates that configuration's relative Boltzmann weight by about 38-fold.
This is a pointwise weight comparison, not a measured occupancy or whole-integral
error. Errors have both signs, so a universal additive offset cannot fix them.

The checked representative geometries have metal–waterO distances2.290–2.787 Angstrom
and nearest waterO–other-heavy-atom contact at least 2.474 Angstrom. They do not show
an obvious gross-clash explanation for the energy-model disagreement. The current
comparison cannot uniquely separate MACE curvature error from nonlinear differences
between its vacuum potential and the native continuum-solvatedDFT target.

Numerical integration fares better: both wider 1F6S domains pass their prescribed
ESS>=128 and independent-scramble F agreement<=0.10 kcal/mol gates. 6IP9 remains
more diffuse: outer-domain weight is15–19%at the largest radius; one largest-domain
ESS is115and fails the128gate. Small scramble differences do not repair low ESS,
and none of these checks establishes complete global basin coverage.

| Pair | Small-domain cheapCa−La Fconfig | Middle | Largest |
|---|---:|---:|---:|
|1F6S11|−0.18059|−0.17298|−0.20510|
|6IP9110|−0.73453|−0.99307|−1.07464|

These are **conditional cheap-potential configurational terms**, with the explicit
common coordinate measure; they are not qualified DFT corrections or calibrated
scores. Their common measure/kinetic factors cancel only in matched water-count
pairs. The apparent pair difference must not be reported as improved biological
accuracy. Full unrounded values and every failed check are in
[the tables](export_v1/TABLES.md) and [the record](export_v1/result.json).

## What remains unavailable

No occupancy probabilities, absolute entropy, DFT-reweighted integral or promoted
score were produced. The representative checks are not an unbiased reweighting
sample. Bound internal-water relaxation/vibrations, quantum consistency,
non-electrostatic solvent terms, other occupancy arrangements and global basin
coverage remain unresolved. Missing terms staynull. Prior stationary/soft-mode
failures and historical calculations remain intact.

The new machinery provides common physical frames with H symmetry, finite coupled
sampling, explicit importance/measure accounting, deterministic native checks and
actual comparison receipts. New protocol IDs:

- `native_r2scan3c_anchored_omol_finite_rigid_water_integral_v1`
- `native_r2scan3c_finite_water_basin_validation_v1`

There is no shared default dispatch change: sampling extends only its isolated
implementation snapshot.13 real-fixture tests pass, zero skips, including all16
actual native receipt/energy extractions. Parser/algebra tests are distinct from
the actually executed6254 MACE calls and16 nativeEnGrad calculations. One initial
preparation halted on a2e-19 floating-point weight difference; the exact chosen
sample and geometry were unchanged. We retained that partial directory and fixed
the numeric comparison without rerunning chemistry.

## Cost and recommendation

| Job | Work | Wall seconds | Allocated CPUs |
|---|---|---:|---:|
|1202081|all MACE samples/anchors|1171|16|
|1202085|1F6S Ca native 4|363|64|
|1202086|1F6S La native 4|393|64|
|1202087|6IP9 Ca native 4|401|64|
|1202088|6IP9 La native 4|423|64|

Total: 119856 allocated core-seconds = 33.2933 core-hours; 1171 GPU-seconds = 19.52 GPU-minutes.
GPU allocation: one H200 / 16 CPUs / 200000 MiB; native: four concurrent 16-rank tasks per
node. Peak reported CUDA allocation 1.92 GiB; peak worker host RSS 2.04 GiB. These are
development costs, not a production speed qualification. Login preparation/tests
were not metered. [Actual receipts](COSTS.json) preserve the measurements.

**Recommendation: retain MACE-assisted water preparation/local response, but do
not use this single-anchor wider-basin model for entropy or occupancy.** Any
continuation needs a better representation of the energy surface across the basin;
more importance samples alone would address the wrong failure. The main scanner
improvement track should proceed independently. No additional jobs were launched.

The older native orientation comparator1201825 was also collected without rerun:
2/4 qualified, one frozen-coordinate drift failure and one nonconverged optimization,
after 18490 s on 64 CPUs. It is separate from all costs/results above; see[status](STATUS.md).
The next runnable operation is the read-only final report replay in[commands](COMMANDS.md).

**Identity erratum (2026-09-19):** Earlier prose mistakenly called 6IP9 GGR.
Both source CIFs identify bovine alpha-lactalbumin. Numerical records and PDB
identifiers are unchanged; see [verified correction](IDENTITY_ERRATUM.md).
