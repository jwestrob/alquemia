# Hydration-network development — approved 2026-09-18

Jacob approved the proposed water-network/occupancy development with “proceed!”.
Standing discretionary authorization applies; no per-analysis permission gate.
Baseline and all completed experiments remain immutable. This is method
development on consumed alpha-lactalbumin structures, not blind validation.

## Scientific question

Can physically prepared, metal-specific water networks improve the electronic
contrast, and can different water counts be compared with consistent exchange
bookkeeping? First resolve water orientations in their actual local network.
A score-dependent choice of deleted water is not an occupancy model.

## Inputs and deterministic preparation

Use the two pinned 1F6S/6IP9 parents and the completed hydration-square results.
For each first-shell water, inventory deposited N/O/S atoms within 3.5 A of its
oxygen (candidate hydrogen-bond contacts, not asserted bonds). Include complete
amide units for contacted backbone atoms and complete supported sidechains for
contacted sidechain atoms; preserve source coordinates and the existing source
atom graph/capping policy. Include contacting deposited waters as frozen outer
waters. Do not iterate the shell outward. Use the union of protein fragment
selections across these identical-sequence structural replicates to keep protein
composition/charge common. Selection is based on geometry, not energy/labels.
Unsupported chemistry or missing covalent neighbors must fail explicitly.

Normalize every included water internally to the already pinned H2O geometry.
Variable occupancy positions are the original two/three first-shell waters;
outer waters are present and frozen in this conditional model. Keep all source
protein/cap coordinates and all oxygen/metal positions fixed in the first stage.

## Stage A: orientation and representation

Two initial orientation sets per structure, identical for both metals:
(1) normalized source orientations; (2) each variable water bisector points away
from the metal, with its plane determined by the nearest noncollinear retained
protein oxygen. No biological label enters this seed construction.

Eight native ORCA constrained optimizations: two structures x two metals x two
seeds. Only variable-water hydrogens move, with water OH lengths and HOH angle
constrained. All other atoms are Cartesian-fixed. Use analytic r2SCAN-3c/CPCM
water/DefGrid3 gradients, no numerical gradients or quantum Hessians. Native
Opt convergence, MaxIter 80 is a numerical optimizer safeguard, not an allocated
compute-time budget. Record every actual SCF/gradient evaluation and convergence;
unconverged searches cannot become optimized endpoints. Record effective SCF
settings and do not attribute unmatched old/new differences solely to geometry.

For the full-water state, compare both starts and both metals, check frozen
coordinates and intact water geometry, and retain all outcomes. Detect initial
or final source/cap clashes, absent gradient support and convergence failures.
Quantify orientation sensitivity and model cost before expanding searches.

## Stage B: occupancy table and bulk-water reference

Use all subsets of the original variable positions: four 1F6S and eight 6IP9
patterns, the same sets for both metals. Outer waters remain fixed. Complete
metal-specific orientation treatment using the same seeds/method; reuse matching
full states. Prepare finite manifests, reuse existing runners and count every
actual calculation, including failed attempts. No production rescore.

Construct an independently specified liquid-water reference at 298.15 K and
1 bar from native gas-phase H2O electronic energy plus traceable thermochemical
terms (gas entropy/enthalpy increment and saturation pressure). Any approximate
zero-point term must be explicit. Do not call isolated H2O+CPCM energy bulk mu.

For state s, Omega=G_s-n_s*mu_water. In the absence of validated bound-water
basin entropy and non-electrostatic solvation contributions, report electronic
state rankings and the additional free-energy terms required to change them;
do not invent G, occupancies, entropy, uncertainty or a calibrated score.
Parameterize missing corrections explicitly; no zero-filled successful result.
The scope permits a subsequent bounded water-position test if needed, with its
physical displacement domain and exact run list recorded before execution.

## Execution and interpretation

Existing isolated environment and ORCA task runner/MPI policy. No GPU planned.
Stage A requests 128 CPUs, eight concurrent 16-rank optimizations; actual cost
will exceed eight single points and is measured from the outputs. No allocated
CPU-time budget or requested wall-time cap. Baseline cost was about 88 s per
small-core endpoint, but expanded-core optimization cost is not established.

Keep old scalar references and decision bands unchanged. Evaluate affinity
claims only on qualified labels and require PQQ fidelity before any promotion.
The current task does not automatically promote a new score, fit missing water
terms to labels, run a broad benchmark, or introduce MACE without validation.
