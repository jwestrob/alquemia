# Native replay locates the large PedE solvent sensitivity

The qualified native ALPB replay attributes **+18.7467 of the +21.2566 kcal/mol**
source/metal contrast to additional solute electronic redistribution **including
its intrinsic-energy compensation**. The contribution evaluated on each
endpoint's vacuum charges is +2.5099 kcal/mol. This locates the sensitive part
of the existing model; it does not establish that this entire response is wrong.

## Exact accounting

For each of sample1/sample3 × Ca/La, the successful replay used the exact saved
geometry, state and eight-decimal vacuum/water charges. Native xTB6.7.1 provides
all ALPB components, including charge-dependent hydrogen bonding, SASA and the
original reference convention. No SCC, DFT or optimization ran.

`G(q) = Gelec(q) + Ghb(q) + Gsasa + Gshift`

`I = E_water − G(q_water) − E_vacuum`

`U = G(q_vacuum) + [G(q_water) − G(q_vacuum) + I]`.

The same signed contrast is applied to both terms:
`sample1_Ca − sample1_La − sample3_Ca + sample3_La`.

| Contribution | Contrast, kcal/mol |
|---|---:|
| Native G on each endpoint's vacuum charges | +2.5098749401 |
| Additional redistribution plus intrinsic compensation | +18.7467424167 |
| Sum from rounded-charge replay | +21.2566173568 |
| Original complete self-consistent correction | +21.2566283014 |
| Replay minus original | −0.0000109446 |

The small closure residual is within the predeclared propagated charge-rounding
budget. G(q_vacuum) still includes source geometry/cavity and vacuum-density
differences; it is not a pure cavity contribution. The bracket is not a pure
polarization measurement or an error estimate. Both come from the same GFN2
model, unlike the independent DFT comparison, which gives +10.0342 kcal/mol.
That nonzero DFT response cautions against dropping electronic response and
calling it a validated correction.

## Actual qualification and technical failures

Job4002 evaluated all eight components after all sixteen water-charge component
checks passed. Gelec residuals are within their analytic rounding bounds; HB,
SASA and shift pass the unchanged tolerances. No tolerance was relaxed.
Water/vacuum-charge evaluations have identical cavity and Born-matrix hashes
within each source/metal case. Exact raw residuals, bounds, charge/state/source
and binary pins are in [EVIDENCE.json](EVIDENCE.json) and its linked collection.

Two earlier technical attempts remain visible:

- 4000: CMake export integration failed before any component calculation.
- 4001: native build succeeded; four water diagonals failed tiny HB/SASA gates.
 Source inspection found the wrapper importing xTB's legacy input conversion,
 while the actual XYZ reader uses mctc_io_convert. Version3 changes only that
 input conversion import. Native solvent parameters were never globally changed.
 The four failed diagonal calls remain preserved; no vacuum calls ran in v2.

The actual final build/eight evaluations used65s on12 allocated CPUs
(780CPU-seconds). Including both preceding attempts, this branch used1668
allocated CPU-seconds, zero GPU, zero new SCC/DFT. Individual native component
calls took about0.07s, excluding building and all previously computed charges.
This is not a new-site classifier runtime or a claim about total screening cost.
Compiler-child peak memory is not relabeled as standalone replay memory.

## Next practical test

A separately named frozen-distribution candidate is defensible to test:
`E_candidate = E_OMOL,vacuum + G_ALPB(q_GFN2,vacuum)`.
It retains the continuum reaction field while omitting solute electronic
redistribution. This is a perturbative hybrid descriptor, not an exact solvent
transfer, corrected affinity or proof that the omitted response was artificial.
No fitted fraction or interpolation is justified by these two consumed sources.

Use existing complete common geometry pools, reselect each metal's geometry
under the new expression, and calibrate only on the designated canonical inputs
using the existing fixed rule. Then freeze those new bands before applying the
entire declared structural-transfer population, including Ca controls and all
fifteen two-Ca folds. Missing charges/cells remain unavailable. This tests actual
predictive usefulness rather than repeatedly decomposing two failures. At this
report's creation that candidate is being prepared, not calculated or promoted.

No PedE molecular/replay job remains active. Delayed wakes for4000/4001/4002
should be acknowledged without rerunning them. Production and old classifications
remain unchanged. Exact future scope and implementation live under the separate
frozen_distribution_candidate directory; root owns execution.
