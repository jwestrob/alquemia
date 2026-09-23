# One verified preparation failure is now scorable

The archived hydrogen repair recovers MMOL1770 La-conditioned sample4. All
three requested models call the repaired source La-supported under their own
existing references. Its9582 atom identities,4895 heavy coordinates, chemical
state and pH7 remain unchanged. Historical225 remains immutable; this is a
separate supplemental structural replicate and adds no independent biology.

| Conditional model | R, model kcal/mol | Own frozen La boundary | Decision |
|---|---:|---:|---|
| Released local static,160atoms | −405448.300670292 | −405459.111382000 | La-supported |
| Fixed union static,184atoms | −405449.184025901 | −405456.684468508 | La-supported |
| Fixed union adaptive,184atoms | −405447.048247691 | −405456.468817541 | La-supported |

R=ECa−ELa uses native MACE plus matched nativeGFN2 ALPB−vacuum. Absolute energy
scales are model-specific; these are electronic classifier contrasts, not
experimental affinities. No new aquo reference, fitted threshold, protonation
state, water or heavy-atom repair was introduced. Local bands come directly
from released RELEASE_PANEL.json frozen_bands; union/static and union/adaptive
use their separate pre-existing canonical references. All calls above are
explicit supplemental preparation-transfer checks.

## Actual execution and physical interpretation

All4 q0 MACE and8 q0 GFN2 endpoints succeeded. Both original-policy union
searches converged: Ca23iterations/24evaluations; La13/14. Their36 fresh MACE
calls reuse actual q0 energies/forces. Both own candidates plus origin were
cross-scored under both metals:2fresh cross-MACE and8GFN2, all successful.
Neither missing-state fallback nor rerun occurred. ZeroDFT or new preparation
force evaluations. Mathematical and operational pool minima agree.

Ca selects its own proposal (composite work−9.855899kcal/mol); La selects its
own (−11.991677), giving deltaR+2.135778. Native contributions are−11.267470 and
−12.662144, solvent-transfer changes+1.411572 and+0.670468 respectively. The
Ca proposal touches the declared displacement boundary; La does not. These
are bounded candidates, not both demonstrated unconstrained minima or free
energy populations. Adaptive scoring was unnecessary for the static call.

## Preparation scope and limitations

The original17 exclusions were inspected against the existing Sept20 diagnoses.
Sixteen actual ions lie25.016–25.232Å from completePQQ and remain unsupported;
no ions were moved and no structures replaced. The sole repairable case was
an independently diagnosed H-minimization failure, already repaired with the
captured same-potential OpenMM system in jobs1203794/1203809. Reuse retains
all4687 hydrogens and removes the sub-.45Å overlaps. No fresh protonation,
H normalization, missing-heavy-atom addition or donor substitution occurred.

The old H-placement potential is approximate and has unusually long H-parent covalent
bonds. Numerical minimization convergence does not validate its chemistry.
This adapter reuses one verified archived repair; it does not automatically
reprepare every incoming structure. Existing source/cap mapping and the frozen
union fragment membership reproduce exactly, and union state/protein signatures
match the canonical group. All original failure records remain accessible.

## Cost and records

Four completed allocations:1210669(16s×32CPU,1GPU),1210670(25s×64CPU),
1210819(41s×32CPU,1GPU),1210820(29s×64CPU). Total **5280 allocatedCPU-seconds,
57 requestedGPU-seconds**,42 new MACE calls and16nativeGFN2 calls. This includes
in-job startup, validation and collection. Local preparation/testing/reporting
is unmetered; earlier archived H repair and method-development cost are reused
and excluded. Peak measured GPU allocation5,254,282,240bytes. The pending
1210670 allocation was routed to the established H200 CPU-only host without
changing resources, priority, tasks or creating a duplicate job.

Exact matrices, components, proposal flags, raw receipts and hashes:
`workspaces/preparation_coverage_20260923/COMPARISON_v1.json`, `COSTS.json`,
`origins_v1/collection_final.json` and `pool_v1/collection_final.json`.
Seven real-fixture tests pass in7.636s (zero skips) and cover archived heavy/H identity, unchanged state/union,
paired coordinates, actual Kinematics replay, isolated optimizer policy, all16
completed native cells, actual q0 reuse, algebra and historical failure retention.
No scientific executable was mocked. Commands and immutable source snapshots
are documented in COMMANDS.md; do not resubmit the completed jobs.

Recommendation: retain this recovered source as supplemental coverage evidence;
leave the16 invalid ion placements and historical225 denominators explicit.
No production default or general preparation policy was changed.
