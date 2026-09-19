# Exact water rotations proposed by MACE, adjudicated by native DFT

The16-geometry proposal check passed before this follow-up was specified. Keep
all prior checks/limitations. This is a contained continuation of the approved
hydration-modeling work, preserving production and the two live native DFT jobs.

## Fixed inputs and method

For1F6S/6IP9 and each metal, start from both original source/radial full-network
geometries. Use unchanged native vacuum MACE-OMOL0-100M float64 with physical
charge/spin. Optimize only the three rotation coordinates of each variable water
about its fixed oxygen. Rotate its intact reference H2O geometry analytically;
all other Cartesian coordinates remain exactly fixed. No water addition,
deletion, translation, donor change, Hessian/entropy model or force-field fitting.

Use scipy BFGS with analytic coordinate derivatives (SO(3) left Jacobian),
gradient infinity norm <=0.001eV/radian, max200iterations as a numerical optimizer
safeguard. Verify the coordinate chain rule on real fixtures before execution.
Do not impose an allocated compute/time budget. Record every model evaluation,
both outcomes, gradient norms and termination message. A nonconverged search
remains explicit; do not call its endpoint a minimum.

Four finite executor tasks (structure/metal), each running two starts with one
resident model: eight optimization searches. One standard GPU with16CPU and
64474MiB host memory. Choose the lower MACE energy independently for each metal
only if both searches meet the same convergence rule; no label-dependent choice.

## DFT adjudication

If all searches converge and exact geometry checks pass, prepare four new
native r2SCAN-3c/CPCM(Water)/DefGrid3/TightSCF analytic energy/gradient endpoints,
one per selected metal/structure proposal. Same expanded cores, charges,
protonation and water counts as the live DFT searches. These are not numerical
DFT derivatives. Compare proposed endpoint energies to both already-computed
DFT initial orientations and later native DFT optima, retaining all discrepancies.

Accept a proposal as improving initialization only when it lowers its own
metal's DFT energy below both original starts. Report DFT water-rotation gradients
directly: no claim of DFT stationarity merely from MACE convergence. A candidate
for substituting the expensive native search must additionally lie within
1kcal/mol per metal of the lowest geometry-qualified native DFT minimum and
within1kcal/mol of its Ca−La relaxation contrast. These provisional accuracy
checks are specified before the new predictions/DFT endpoints, not fitted to
biological labels. Failure does not authorize thresholds to be loosened.

Because native DFT optimization has small constraint drift, its final geometry
may fail the already-declared check. Such a comparator remains unqualified;
no MACE substitution claim follows until a matched exact-geometry DFT comparison
exists. No routine rescore or new benchmark classification is launched here.
Missing water basin entropy and occupancy free energies remain unavailable.
