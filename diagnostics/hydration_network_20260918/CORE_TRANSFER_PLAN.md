# Isolate water preparation in the original scoring cores

The expanded-network MACE proposals passed DFT adjudication: all four lower the
energies of both original starts, with orientation deltaR+21.69/+29.69kcal/mol.
Their nonzero target DFT rotational gradients (3.47–5.42kcal/mol/radian) mean
these are improved preparations, not claimed DFT minima.

The next contained test answers a narrower, scanner-relevant question: does
correcting only water H positions improve the original alpha/GGR electronic
ordering without adding protein atoms to the scored core or changing the
original SCF policy? This separates preparation utility from a new energy core.

## Fixed comparison before execution

Transfer the selected water H coordinates, by source-atom identity, from the
same expanded-context MACE proposals into the original40/43-atom amide-v3 alpha
cores. Retain every original protein/cap/metal/water-O coordinate, water count,
formal charge (La0/Ca−1) and atom order. The context chooses orientations; it is
not added as an environmental energy term. Restore the exact original native
r2SCAN-3c/CPCM/DefGrid3 single-point recipe, including its original SCF policy.
Four new quantum single points, two per structure. No optimization, Hessian,
new MACE call, altered occupancy or label-dependent seed choice.

Reuse the archived repaired52-atom GGR pair after verifying its original
coordinates/method/receipts. It has no waters, so this preparation operation is
an identity operation. Do not rerun it. Likewise no new PQQ calculation is
needed to establish that a water-only coordinate update leaves dry inputs
unchanged; preserve released PQQ inputs/results and classification policy.

Report both original and prepared alpha−GGR differences in R directly (the
common reference cancels). Existing experimental direction expects alpha above
GGR. Predeclare success of this development comparison as both structural
replicates having positive differences; no threshold fitting or best-replicate
selection. This is one consumed alpha biological group versus one GGR group,
not independent validation or a new universal decision threshold.

Also compare the prepared original-core result to the already executed
expanded-core result, retaining the different core/SCF policies explicitly.
A useful outcome is a preparation improvement at roughly original per-endpoint
quantum cost plus short MACE proposal cost. No default promotion or production
rescore, and no bound-water free energy or occupancy probability is assigned.
