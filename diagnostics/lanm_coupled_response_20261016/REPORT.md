# The scaffold changes basin selection; fixed-source dependence remains

16 October 2026. Jacob recommended starting with monomers. We recovered Mex-LanM
as the natural-protein reference and the already-consumed LBT3 peptide as a small
opposite-preference control. LBT3 is not a second natural LanM. Exact primary
evidence, assay differences and source limitations are in MONOMER_TARGETS.md.

## What the new calculation says

Eight native MACE singlepoints compare the whole256-atom LBT3 peptide with a
48-atom donor-fragment control, using both deposited La7CCO and Lu7CCN model1
geometries and both La/Dy identities. Shared physical coordinates are exact.

At either fixed source, adding the scaffold changes the paired Dy-minus-La
contrast by only−0.033 or+1.608kcal/mol. The difference between source contrasts
remains72.049kcal/mol, versus73.690 for the pocket. Thus merely representing the
whole peptide does not resolve starting-geometry sensitivity.

However, the scaffold changes **selection between the two conformations**:

| Representation | La minimum | Dy minimum |
|---|---|---|
| Donor fragments | La-derived7CCO | La-derived7CCO |
| Whole peptide | La-derived7CCO | Lu-derived7CCN |

The whole-peptide common-pool comparison against the three compatible saved Mex
site cores has the qualitative expected relative direction in all three entries
(+19.623,+27.714,+19.598kcal/mol). The pocket comparison has mixed directions
(−3.810,+4.282,−3.835). These are balanced conditional exchange contrasts with
element offsets cancelled. They are **not three experimental site labels, a
measured affinity, a success rate or a population estimate**. No sites are summed.
Mex's protein-level qualitative affinity and LBT3 ITC use different conditions.

This makes scaffold-sensitive accommodation worth testing, but does not establish
a working LanM profiler. The dry vacuum model, assumed charged termini, NMR
coordination restraints and source dependence remain important limitations.
The whole-versus-core contrast combines scaffold, caps and the omitted D13 charge;
it does not isolate a pure mechanical contribution. Full matrices and exact pins
are in lbt3_static/RESULT.json and REPORT.md.

## Dy parameter repair now has real runtime evidence

The La-corrected nativeGFN2 profile still assigned three Dy arrays to the wrong
declared shells. We created a new full La+Dy profile that changes only those
explicit Dy arrays on top of the existing La correction; all other fields and
shell order are preserved. The old profile and failures remain unchanged.

Eight new48-atom real-core controls completed: two metals × vacuum/water × native
ORCA/standalone xTB. Every endpoint passed its original process, state, population
and convergence criteria. All four absolute energy/population comparisons passed
the existing1e−6Eh/1e−5e thresholds; both solvent-component comparisons also passed.
Residuals are reported individually in dy_parameter_qualification/RESULT.json.

This qualifies the tested parameter import/runtime comparison, not Dy affinity,
open-shell electronic response or the full protein. NativeGFN2's f-in-core
effective singlet is distinct from MACE's physical sextet input. No production
registry/default or old solvent score has been silently changed.

## Executed work and costs

| Work | Actual execution | Scheduler allocation |
|---|---|---|
| Hans coupled geometry/classical preparation |4 admissible proposals,8 classical cells; no QM|12460:13s,1CPU|
| LBT3 preparation |2 sources,8 endpoints;2 real missing-H helper calls|12469:3s,1CPU|
| Native MACE |8 fixed cells,all complete|12471:26s,12CPUs,1H200|
| Saved comparison |v1 pin-path parser failed;v2 completed, no molecular repeat|12478:3s;12484:1s,1CPU each|
| La/Dy parameter qualification |8 endpoints,all admitted and parity-passing|12488:19s,4CPUs|

Total for these listed jobs:408allocatedCPU-seconds and26allocatedGPU-seconds.
These are allocations, not measured utilization. Individual measuredCPU/memory
and exact scheduler records are in EXECUTED_ACCOUNTING.tsv. MACE application
time22.480s and qualification18.138s are separate from scheduler elapsed time.
Preparation's internal missing-H iteration/force-call count is unavailable;
it is not described as zero-energy work. Historical Mex endpoint costs are not
included in incremental work. No new DFT or experimental observation was generated.

## Solvent and response follow-up

The scout is now complete: four MACE analytic forces replay the earlier energies
exactly; all eight whole-peptide GFN2 energy/gradient endpoints pass admission.
Both selected conformations survive solvent, but the paired source gap remains
−72.966739 kcal/mol. Details and actual costs are in lbt3_solvent/REPORT.md.

A graph-defined catalogue covers backbone, sidechain and metal motions in the
intact peptide. Saved-force job12512 selects metal translation and four intact
Asp7/Asp3/Glu11 torsions by one common force-based rule across both sources and
metals. The nominated Asp7 carboxylate direction has opposite La/Dy energy slopes.
This is a local force diagnostic, not a finite relaxation or affinity result.

The actual two-step directional derivative check passed for both metals and
its paired response. Four subsequent own-metal searches and the complete shared
pool are now complete. The source-local Dy-minus-La contrast spread falls from
72.966739 to57.895436kcal/mol; La and Dy still select different source-derived
basins. No robust preference or experimental transfer is established.

All four searches reach the original0.20Å atom-displacement envelope at Asp7
oxygen atoms while separate covalent/contact checks pass. A uniformly declared
range-only continuation changes that envelope to0.40Å and reuses the exact saved
endpoints; original-source references, target/torsion limits and all model/state
assumptions stay fixed. See LBT3_RANGE_EXTENSION_PLAN.md. Added backbone freedom
is a separate possible test, not silently included in this continuation.

Read lbt3_accommodation/REPORT.md for the full paired matrix, force diagnosis
and actual costs. The completed response/pool round used32newMACE/64GFN calls,
385allocatedGPU-seconds and4621allocatedCPU-seconds. All results are preserved;
there is no new DFT campaign, production change or library rescore.

The overall LanM profiler goal remains active and unmet. A failed local branch
will change the model choice rather than be hidden or presented as completion.
