# Fresh prepared-context composite scanner — frozen 2026-09-20

Jacob said “pursue!” after the completed compact solvent comparison. Parent
assigned this contained execution/timing pilot, in parallel with independent
numerical qualification. No production/default change.

## Scientific scope

Exactly four existing consumed structures, in this order: 1H4I,
q9z4j7-pqq-la_model, 1F6S, 1GLG. Use their exact complete-context Ca/La coordinates,
charges, multiplicities and water states from the compact-solvation inventory.
No new geometry, water selection, protonation, DFT, model training or threshold.

Fresh calls: **8 native OMOL scalar endpoints +16 native GFN2 endpoints**
(four per site: each metal in vacuum and ALPBWater). No old result satisfies a
fresh timing task. No automatic retries. Preserve partial failures.

Native OMOL uses the exact 100M checkpoint and float64, no force calculation
requested. The existing captured energy-only worker computes the same scalar
Hamiltonian; this is not another scoring model. Compare every actual scalar to
its archived gradient-enabled source and exact checkpoint. Existing OMOL
endpoint repeat tolerance: 0.01 kcal/mol. Composite execution repeat target:
0.02 model-kcal/mol, small relative to the 5.076 class gap. Report all differences,
including exact decisions at fixed bands; tolerances never widen decision bands.
This allocation repeat does not replace the independent solver qualification.

GFN2 is unchanged installed native ORCA6.1.1 GFN2, native mixer, SmearTemp300,
NoAutostart, default parameter file exported, vacuum versus ALPB(Water). Native
parser verifies actual valence counts, atomic charges and total energies. No
inferred ALPB decomposition or invented solvent radii.

## Score and frozen decisions

E(M) = E_OMOL,vac(M) + E_GFN2,ALPB(M) - E_GFN2,vac(M).
R = E(Ca)-E(La), converting eV and Hartree separately exactly once. No aquo S.
Use existing published context composite bands exactly:
Ca-supported R <= -405464.18774828;
La-supported R >= -405459.1113819997; otherwise inconclusive.
Q9Z4J7 is the Ca calibration maximum; tiny floating differences can change its
literal band decision. Do not shift thresholds to conceal that sensitivity.
PQQ bands apply only to compatible PQQ preparation. Alpha/GGR is a relative
direction only, one biological comparison, not two absolute class assignments.

## Resources and measured scope

One H200 GPU node, 1GPU,32CPUs,200000MiB. Four sites sequentially. Each native
MACE pair runs sequentially using the existing worker, then four concurrent
8-rank native GFN2 tasks through the existing ORCA executor. This uses all32CPUs
for the solvent stage. Model loading occurs for each native endpoint; software
startup, validation and final collection are counted in total operational cost.
GPU reservation during CPU stages counts in allocated GPU-seconds.

No arbitrary project runtime cap. The manifest is finite. Record actual site
endpoint latency, process/model/inference breakdown, total score-path latency,
allocated core/GPU-seconds, host/GPU memory and Slurm receipt. This starts from
prepared context XYZ; folding, source preparation and calibration are excluded.
No matched DFT speedup can be claimed from this GPU-only pilot.

## Interface

The pilot uses the archived inventory for comparisons. A separate `prepare-pairs`
operation accepts explicit prepared pairs and the pinned native model/software,
context preparation provenance and frozen calibration. Previous MACE energies
are optional. This is prepared-input research use, not automatic raw-fold
scanning or global promotion. Only the four-site pilot is authorized here.
