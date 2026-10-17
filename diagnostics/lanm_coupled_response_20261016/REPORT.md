# Coherent monomer response works; reliable La/Dy preference remains unproven

## Latest result — 17 October 2026 UTC

Mex-LanM is the natural monomer reference; LBT3 is a measured peptide control,
not a second natural LanM. The working energy combines local MACE chemistry,
a GFN protein-field correction, the retained protein mechanics and fixed-charge
ALPB solvent. Whole LBT lies inside the electronic region, so it has no separate
external-field or MM term. State, source and assay limitations remain explicit.

The Mex finite continuations preserve the checked covalent geometry and lower
each metal's own energy. Dy benefits more than La. The full paired comparison
retains the expected conditional Mex-La/LBT-Dy direction, but its shared-menu
margin is only1.019378kcal/mol; some other declared menus are about0.09. The two
LBT source-local contrasts still differ by28.262193kcal/mol. Positive signs here
do not establish a reliable preference classifier or an affinity scale.
See [the complete comparison](monomer_continued_comparison/REPORT.md) and
[Mex response](mex_continuation/REPORT.md).

The matched-solvent water rescore weakens the earlier water lead: its two
available conditional responses are-0.787/-0.740kcal/mol, within the existing
1kcal/mol interpretive resolution; four original geometry rejections remain.
See [matched water](lbt_matched_water/REPORT.md). No hydration expansion follows.
The latest molecular round used274GPU-seconds and3299allocatedCPU-seconds
including saved reporting; throughput is adequate for contained profiling tests.

The next declared experiment tests a real asymmetry: LBT's existing geometries
were generated under the former solvent objective and only rescored under the
current model. [Matched-response plan](LBT_MATCHED_RESPONSE_PLAN.md) retains all
four starts and unchanged physical limits. It is in preparation, not a result.
A separate [saved electronic-substitution diagnostic](ELECTRONIC_SUBSTITUTION_PLAN.md)
will distinguish sensitivity to MACE electronics from the shared mechanics/solvent.
The second Mex deposit requires new construct/occupancy/hydration assumptions;
see [its feasibility record](MEX_SECOND_SOURCE_FEASIBILITY.md). PQQ/defaults remain
unchanged. The active profiler goal is not complete.

## Historical progression under the earlier models

## Current result and next decision

**No reliable monomer preference profiler yet.** On the consumed LBT3 control,
added donor and backbone freedom reduces dependence on the starting structure
beyond equal extra search in the old coordinates. The remaining discrepancy
is still too large to claim a trustworthy La/Dy preference.

| Response model | Source-local Dy-minus-La spread, kcal/mol |
|---|---:|
| Original fixed sources | 72.967 |
| Original bounded response | 57.895 |
| Uniformly enlarged displacement envelope | 45.805 |
| Matched further search, original seven coordinates | 40.452 |
| Matched further search, three added donor/backbone coordinates | 29.509 |

These are computational source-sensitivity measures, not measured binding
energies. All original sources and failures are retained. The matched comparison
is in [lbt3_basis/REPORT.md](lbt3_basis/REPORT.md).

Two affordable components now work on actual fixtures: native fixed-charge
solvent forces pass their 26-call derivative/replay test; Mex's four local GFN2
endpoints converge with and without its actual surrounding protein field.
The latter field stabilizations are similar for La and Dy at the fixed origin.
See [solvent forces](frozen_solvent_force/REPORT.md) and
[Mex local field](mex_local_field/REPORT.md). Their scientific role is to enable
coupled protein response, not to supply a ready preference correction.

The finite explicit-water panel now supplies a concrete local response: moving
one modeled water inward favors Dy relative to La by about2.5kcal/mol in both
7CCO orientations. All four7CCN inner placements remain geometry-rejected, so
this is a conditional mechanism, not source-robust preference. See
[the completed water report](lbt3_water/REPORT.md). The fixed panel is closed.

The Mex whole-system assembly now works: two new localMACE and two globalALPB
calls completed in13s on oneH200/12CPUs, reusing existing localfield and MM data.
An independent saved comparison with matching archived PBE0 gradients agrees
on all six predeclared metal-differential response signs. Magnitudes differ,
and weak own-metal/common directions can reverse. See
[Mex assembly](mex_static/REPORT.md) and [reference response](mex_reference_response/REPORT.md).
The subsequent16-cell total-force test passes both external-site and cap-anchor
interfaces in69s on oneH200/12CPUs. The first coupled response12717 then passes
all four paired cells and covalent/contact guards in32s on the same resource
shape. La and Dy each select their own proposal, lowering their energies by
0.252489/0.354867kcal/mol and changing the paired contrast by−0.102379. The tiny
step mainly moves the target; residual scaffold/donor/metal loads remain. A
separate finite same-model continuation is being prepared. These results do not
validate affinity, electronic states or source robustness.
The matched fixed-FF solvent comparison now gives the expected conditional relative
ordering for **both** LBT3 starting sources: Mex–La plus LBT–Dy is favored, by
+2.533 or +68.101kcal/mol. The common two-source pool chooses7CCO for both metals
and gives+2.533. This is a useful directional result, but the large magnitude
sensitivity and different assays/partitions prohibit quantitative affinity or
source-robust classifier claims. Only four new native solvent calls were needed,
in5seconds/oneCPU; see [matched comparison](lbt_matched_model/REPORT.md).
Rescoring the complete existing LBT accommodation pool under that matched model
retains the expected direction in every predefined source/menu. Added donor and
backbone freedom reduces the source contrast gap from65.569 to28.262kcal/mol;
equal extra search in the original coordinates leaves38.435. Including Mex's
first response step retains every sign and gives+2.408 for the full shared menu.
This is a concrete accommodation benefit, with substantial residual magnitude
sensitivity. See [matched LBT pool](lbt_matched_pool/REPORT.md) and
[combined comparison](monomer_saved_comparison/REPORT.md). Only32 new ALPB calls
were needed,13seconds/oneCPU; the combination used saved data only.
No newDFT, production/PQQ change, dimer prerequisite or broad library campaign.
Earlier results below are historical evidence.


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
