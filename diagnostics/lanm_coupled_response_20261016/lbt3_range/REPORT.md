# More donor movement helps; substantial source dependence remains

The uniform range-only continuation is complete for both metals and both LBT3
sources. The paired source-local contrast spread falls from 57.895436 to
45.804712 kcal/mol. Including the original unaccommodated comparison, the sequence
is 72.966739 → 57.895436 → 45.804712 kcal/mol. This is improved consistency within the
same molecular model, not a reliable biological preference or affinity result.

## What physically changed

Only the source-relative all-atom displacement envelope changed 0.20 → 0.40 Å.
The original source remained the contact/bound baseline; target 0.20 Å, torsions
±0.20 rad, atom step 0.05 Å and all covalent/state/energy criteria were retained.
Every search began at its exact saved endpoint with saved forces and original
metric; no origin was recalculated. All final geometries pass admission.

|Source/metal|New accepted steps|Additional own work (kcal/mol)|Stop|
|---|---:|---:|---|
|7CCO/La|2|−0.266482|Retained torsion bound|
|7CCO/Dy|12|−0.687854|Declared accepted-step limit|
|7CCN/La|4|−15.277206|Expanded source-displacement bound|
|7CCN/Dy|5|−3.607853|Retained torsion bound|

The largest additional response occurs when La starts from the Lu-conditioned
7CCN structure. This supports a remaining mechanical limitation in that model;
it does not show that the deposited NMR structure is erroneous. Finite proposal
endpoints retain active forces and are not stationary minima.

## Fair paired comparison

All four new terminal points are distinct. The original six candidates/twelve
cells are retained, yielding ten candidates/twenty complete paired cells after
only four additional opposite-metal composites. Each metal sees the same menu.
The global menu still selects a 7CCO-derived La geometry and 7CCN-derived Dy geometry;
the global Dy-minus-La contrast changes −3.341370 kcal/mol. A common global menu
removes source choice by construction, so source-local menus remain the actual
robustness check. No hydrated-ion reference, population, threshold or Mex-site
experimental label is inferred.

## Execution and cost

Jobs 12544/12545/12546 all completed: 33/246/55 scheduler seconds, each on 12 CPUs
and one H200. This round used 28 new MACE/56 GFN evaluations, with no failed calls or
repeated starting endpoints. Allocation totals 4008 CPU-seconds/334 GPU-seconds.
Measured TotalCPU is 40.233/359.714/74.641 seconds respectively; these values
are distinct from resource allocation. Application and individual backend
costs, actual memory and exact pins are retained in RESULT/ACCOUNTING files.
Cumulative monomer branch through 12546: 9950 allocated CPU-s/767 allocated GPU-s,
including earlier preparation/support failures, excluding old Mex molecular work.

## Next discriminant

LBT3_COMMON_BASIS_PLAN.md declares a matched seven-versus-ten-coordinate test.
Both arms use the same range endpoints, budgets and feasible-direction fallback.
The candidate adds the already-mapped Asp3/Glu11 donor torsions and Ile10
backbone hinge selected from earlier saved-force evidence. This separates added
mechanical freedom from merely giving the old coordinate set more iterations.
No basis/control energy has run at this checkpoint. No new DFT, dimer campaign,
production change or PQQ work. The natural-protein profiler goal remains unmet.
