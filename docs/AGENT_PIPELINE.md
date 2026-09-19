# Alquemia: current operating guide for agents

## Current: coupled water response complete — 2026-09-19

27native DFT checks and all MACE tasks completed. All11trial steps lower native
energy and pass energy checks;3/4endpoints meet water-coordinate minimum criteria.
1F6SCa remains nonstationary (projected max0.70245 versus0.4 limit) and retains
its failed soft-mode check.7/8initial direction checks pass; all4response paths pass.
This supports cheap water-relaxation proposals, not harmonic entropy or occupancy.
Thermal model extents exceed the validated local domain; missing free-energy terms
remain null. Baseline/default/PQQ unchanged.71real-fixture tests pass, zero skips.

Read [coupled report](../diagnostics/hydration_basin_20260919/REPORT.md),
[actual tables](../diagnostics/hydration_basin_20260919/export_final_v1/TABLES.md), and
[commands](../diagnostics/hydration_basin_20260919/COMMANDS.md). The finite continuation is
complete; do not launch a third recentering round or rerun completed jobs. One
qualified6IP9La endpoint was explicitly reused. Exact new allocation cost:
178496core-s/840GPU-s. Final report uses frozen analysis_implementation_v4.
Four initial collection equality failures were recovered without new chemistry.
The next research issue is wider water basins and complete state accounting,
not another identical local step or an invented entropy offset.

Older native orientation comparator1201825 for6IP9 remains live; preserve it and
collect it when finished.1201824 is complete with recorded geometry-tolerance
failures. These older jobs are separate from the completed coupled-response pilot.
Prior radial/occupancy/preparation results remain preserved; older summaries below
are historical checkpoints. Current work and vault note are recorded in SESSIONS.

**Updated 2026-09-18.** Start here for new work. Recover the latest user
instructions and [checkpoint](../diagnostics/mace_discriminator_goal_20260916/CURRENT.md)
for analysis authorization; a newer user instruction takes precedence over the
local AGENTS note. This guide supersedes older
operational/status prose; dated experiments and their numerical records remain
immutable. The baseline remains the default. Broad La/Ca affinity discrimination
has not been established. The [current MACE goal](../diagnostics/mace_pqq_utility_20260918/GOAL.md)
tests the frozen scorer on Jacob's labelled PQQ references. Latest steering:
accuracy and robustness take priority over speed. The [comparison is complete](../diagnostics/mace_pqq_utility_20260918/REPORT.md):
MACE retains 25/25 calls at 1.859x median speed, with a separate DFT reproduction
exception preventing the strict combined qualification. Earlier grouped-model trials and
further runtime variants remain paused; the approved water-state work is below. Matching calibration
calls preserves fidelity but does not demonstrate improved generalization.
Production remains unchanged.

**Current water-state experiment is complete:** all12 arrangements across1F6S
and6IP9,24 endpoints (20new/4reused),16 water-addition contrasts. All additions
shift the electronic contrast towardLa after contextual water-H preparation;
Ca/La fixed-count site rankings differ and water coupling is substantial.
This is a consumed mechanistic experiment, not a new accuracy or equilibrium
occupancy claim. Missing bound-water free-energy terms remain unavailable.
Production/PQQ unchanged;49fixture tests pass. Jobs1201908/1201910 completed.
[Result](../diagnostics/hydration_occupancy_20260918/REPORT.md),
[commands](../diagnostics/hydration_occupancy_20260918/COMMANDS.md),
[local-motion plan, now executed](../diagnostics/hydration_occupancy_20260918/NEXT_MOTION.md).
Native orientation1201824 is complete;1201825 remains live.

**Earlier classifier experiment:** Jacob approved the archived-feature classifier
comparison; it is now complete. Holding homolog groups together, DFT and
DFT+structure each give 25/25 PQQ calls, while adding the two MACE summaries
gives 24/25. The direct-affinity set contains only two biological groups, so
out-of-group validation is unavailable. No accuracy improvement established
and no default changed. The reusable inventory, feature, evaluate, export and
score operations are in [commands](../diagnostics/site_classifier_20260918/COMMANDS.md);
[results and limitations](../diagnostics/site_classifier_20260918/REPORT.md)
describe the six research models. Do not treat training fits as validation or
restart paused model campaigns from older notes below.

**Hydration mechanism pilot completed:** all 14 native DFT endpoints succeeded
(job 1201801; 333 s allocation wall, 21,312 allocated core-seconds). Normalizing
water H geometry shifts the two alpha contrasts +4.93/+2.48 kcal/mol toward La;
individual water effects range from −19.95 to +5.01. This identifies substantial
water-specific coupling, not a validated occupancy correction or accuracy gain.
[Results and operations](../diagnostics/hydration_square_20260918/REPORT.md).
That square pilot has no live jobs. Baseline and old inputs are unchanged. Jacob
removed the standing per-analysis approval gate again on 2026-09-18; historical
pauses above do not reinstate it. Recover the current checkpoint for next work.

**Hydration-network continuation is running:** jobs1201824/1201825 contain eight
native DFT water-orientation searches on source-defined expanded networks in
1F6S/6IP9. Missing hydrogen-bond partners are restored; protein/water oxygens
stay fixed. The independent liquid-water reference completed as1201831.
A preparation-only accuracy improvement is now available below; native DFT minima
and occupancy remain unestablished. Recover
[operations](../diagnostics/hydration_network_20260918/COMMANDS.md) and the
[Opt collection correction](../diagnostics/hydration_network_20260918/COLLECTION_POLICY.md)
before collecting; an absent final `.engrad` does not justify rerunning a
converged native optimization. Stage B remains prepared, unsubmitted, and needs
a fresh preparation with the corrected artifact expectation.

**A useful MACE water-proposal result:** native OMOL reproduced8/8 local DFT
energy-change directions and4/4 seed rankings, with median water-rotation gradient
cosine0.9913. All eight subsequent exact rigid-water optimizations converged.
Four selected configurations passed DFT adjudication as1201853; native
DFT searches remain live for comparison. [Proposal result](../diagnostics/hydration_network_20260918/MACE_PROPOSAL_REPORT.md)
and [follow-up scope](../diagnostics/hydration_network_20260918/MACE_OPTIMIZATION_PLAN.md).
This is local proposal utility, not an occupancy or broad accuracy claim.

**Prepared original-core scoring now repairs the alpha/GGR comparison:** only
water H positions change; the original40/43-atom cores and native DFT recipe stay
fixed. Alpha−GGR1GLG changes from−14.56/−17.54 to+10.92/+14.57kcal/mol. All six
comparisons to three archived GGR structures now have the expected ordering;
weakest margin0.655. One consumed biological comparison, not broad validation.
Practical route:96GPU-s for preparation plus118s on64CPUs for four DFT endpoints.
Production/default and PQQ remain unchanged;40fixture tests pass. Read the
[complete result](../diagnostics/hydration_network_20260918/REPORT.md).

**Earlier research result:** the complete frozen-response conductor correction
provides no overall improvement over the preceding GK hybrid:4/12 directional
comparisons remain correct. Parvalbumin EF improves, alpha-lactalbumin regresses.
All56 new solves complete and paired grid changes are small; earlier endpoint
and rotation qualification failures remain. Close this challenger without wider
panel expansion. [Results and replay](../diagnostics/mace_omol_20260917/CONDUCTOR_DISCRIMINATION_REPORT.md).

The paper challenge is complete: A0A7 and HEW5 agree with the qualified
experimental direction; RTX fails all three GGR comparisons (6/9 comparisons,
2/3 new domain groups). All22sites/44endpoints passed numerical checks.
[Paper results, caveats and figure](../diagnostics/mace_omol_20260917/KHOURY_BENCHMARK_REPORT.md).

The full3D metal-response pilot also completed: all four native GGR energy and
gradient checks pass, with0.176kcal/mol partition sensitivity. All alpha optima
exceed the frozen0.20A domain, so no alpha correction or discrimination gain is
released. [Response result](../diagnostics/mace_omol_20260917/METAL_RESPONSE_REPORT.md)
and [operations](../diagnostics/mace_omol_20260917/METAL_RESPONSE_COMMANDS.md).
The [bounded-response test](../diagnostics/mace_omol_20260917/BOUNDED_METAL_RESPONSE_REPORT.md)
completed: all28 energy/gradient checks pass. Both alpha-minus-GGR margins improve
(−14.92→−8.30 and−17.31→−6.07kcal/mol), but both directions still fail. This
supports a conditional response component, not broader classification accuracy.
Four new analytic DFT and eight short MACE calls took217CPU-allocation wall seconds
and81GPU seconds; GGR/native grids were reused. No free-minimum or entropy claim.
[Commands](../diagnostics/mace_omol_20260917/BOUNDED_METAL_RESPONSE_COMMANDS.md).
**Earlier development result: matched vacuum hybrid response improves raw alpha/GGR ordering from
2/4 to 4/4**, with actual margins +8.77 to +13.20 kcal-equivalent. All eight
native energy checks pass, but two radial signs and the final 2.024 kcal
partition difference fail qualification; qualified scores stay null. These
consumed comparisons represent two biological groups. Production remains default.
[Result and costs](../diagnostics/mace_omol_20260917/MATCHED_H_RESPONSE_REPORT.md)
and [operations](../diagnostics/mace_omol_20260917/MATCHED_H_RESPONSE_COMMANDS.md).
**The subsequent GGR structural transfer failed:** all eight new comparisons
remain misordered on2FW0/2FVY. Across all structures:4/12 raw directions versus
2/12 static; no new structural robustness. New partition changes are4.87–5.93
kcal-equivalent and qualification remains unavailable. Stop expanding/tuning
this metal-only candidate; preserve the production baseline and research tools.
[Final transfer result](../diagnostics/mace_omol_20260917/HYBRID_GGR_TRANSFER_REPORT.md)
and [operations](../diagnostics/mace_omol_20260917/HYBRID_GGR_TRANSFER_COMMANDS.md).
The [physical donor-coordinate diagnostic](../diagnostics/mace_site_response_20260918/REPORT.md)
now projects all saved matched-hybrid gradients onto actual donor torsions and
peptide motions. All geometry checks pass; donor sensitivity is substantial in
the failed GGR cases. Zero new energy calls; no response correction is available.
The [coupled donor/metal path test](../diagnostics/mace_site_response_20260918/COUPLED_PATH_REPORT.md)
is complete: 16 native DFT and 96 MACE calls retain 4/12 raw directions and no
qualified comparisons. All eight additional GGR-structure comparisons still
fail; coupled margins are worse than the earlier metal-only response. All
energy-prediction checks pass, while six gradient checks and all three
representation checks fail. Close this path as an accuracy improvement.
[Operations](../diagnostics/mace_site_response_20260918/COUPLED_PATH_COMMANDS.md).

The [whole-protein simultaneous-substitution test](../diagnostics/mace_collective_20260918/REPORT.md)
also completed: eight MACE calls, all 23 numerical/algebra checks pass, but the
author-domain challenge remains 6/9 and RTX worsens. Original site vectors and
evidence strata are retained. This static occupancy descriptor does not model
experimental folding or titration. Neither addition changes the production
baseline or supplies a broadly validated replacement. Both results have vault
notes, actual receipts and replay commands. Those two studies have no live jobs.
The [whole-protein charge-group POLAR trial](../diagnostics/mace_charge_groups_20260918/REPORT.md)
completed 24 MACE and 24 solvent calls. All numerical/grouping checks pass;
PQQ ordering improves from -17.63 to +19.96 model kcal, but all six alpha/GGR
directions still fail (1/7 total). Vacuum component directions are promising,
but its grouping shifts exceed the retained tolerance. No broad validation or
default change. [Completed operations and failures](../diagnostics/mace_charge_groups_20260918/COMMANDS.md).
The [carve-free typed-contact grouped-vacuum trial](../diagnostics/mace_group_transfer_20260918/REPORT.md)
completed55evaluations:17/22 raw directions,including all7original comparisons.
Khoury6/9 andparvalbumin4/6 retain their previous supporting-family failures.
All61numericchecks pass; threegrouping shifts2.51–3.28 failthe retained2gate,
althoughall22directional outcomes remainunchanged. No absolutecalibration or
defaultchange. [Completed operations](../diagnostics/mace_group_transfer_20260918/COMMANDS.md).
The [same-model PQQ calibration](../diagnostics/mace_group_canonical_20260918/REPORT.md)
completed50calls: class gap−213.316367 and106/154 pair orderings, no valid bands.
The short learned readout separates25/25; electrostatics correlates with protein
charge and destroys total-score separation. This is a diagnostic, not a new
classifier. The previous **whole-chain masked OMOL** classifier still separates
25/25 with its own9.108591 gap; do not describe it as a carved-core result.
The [large-checkpoint transfer](../diagnostics/mace_group_large_20260918/REPORT.md)
completed57calls:17/22, the same failures as medium, with about twice the observed
campaign GPU allocation. All69numerical checks pass; three grouping checks fail.
No larger canonical expansion or default promotion. Both jobs are complete.
Recover the [live checkpoint](../diagnostics/mace_discriminator_goal_20260916/CURRENT.md)
and inspect actual jobs before submitting research work.

## 1. Select the actual protocol

| Path | Protocol / implementation | Current use and interpretation |
|---|---|---|
| Existing automatic generic inbox | `generic_vertical_exchange_native_r2scan3c_v2`; `scripts/carve_generic.py` | Retained production/reference behavior. Its coordinating-backbone fragment defect is documented; it is not silently repaired in place. |
| Existing automatic PQQ inbox | `pqq_vertical_swap_r2scan3c_native_cpcm_typed31_v2`; `scripts/carve_with_pqq.py` | Buffered distance-selected PQQ core. This is a different protocol from the fixed-core calibration below. |
| Canonical PQQ benchmark | `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`; `diagnostics/pqq_pmdh_fixed_core_calibration_20260914/fixed_core_carver.py` | Exact homologous core, mapped catalytic partners, frozen PQQ state and no waters. Only compatible preparations inherit its released bands. It is not the inbox default. |
| Repaired generic benchmark | `generic_peptide_amide_vertical_native_r2scan3c_v3`; `scripts/affordable_peptide.py` | Explicit opt-in source-graph peptide-amide repair. Preserves source heavy coordinates and actual bonded amide N; merges overlaps. No calibrated absolute bands. |
| Global electrostatic research | `vacuum_r2scan3c_mbis_global_tabi_electrostatic_v1`; `scripts/global_electrostatic.py` | Completed physical pilot failed partition, mesh-refinement and rotation gates. Vacuum QM + direct protein Coulomb + whole-protein reaction field. No compatible aquo reference or calibrated decision; no default change. |
| Density/field diagnostic | `native_r2scan3c_permanent_field_density_diagnostic_v1`; `scripts/density_embedding.py` | Complete: exact-density coupling reduces the partition discrepancy by 8.06 kcal/mol; core response changes it by another −0.97, leaving 9.67 without solvent. Uniform native CHELPG improves potential/coupling agreement on all four consumed states. No solvent score or calibrated class. [Results](../diagnostics/density_embedding_20260916/REPORT.md). |
| Native interaction diagnostic | `native_r2scan3c_asp303_interaction_eda_v1`; `scripts/interaction_decomposition.py` | Native ghost-basis references already fail the frozen equivalence check. Retry 1199986 was running at checkpoint, with unstable Asp fragment SCFs; no complete decomposition claimed. An automatic collector writes its terminal result. [Status and completion location](../diagnostics/density_embedding_20260916/INTERACTION_STATUS.md). |
| Gradient/response research | `scripts/affordable_response.py`, `scripts/ggr_sensitivity.py`, `scripts/mace_output_audit.py` | GGR analytic-gradient checks completed. Saved full-protein MACE outputs now export source-mapped direct grad(E_Ca−E_La), charges and checkpoint comparisons. These are not hybrid gradients. The [MACE curvature screen](../diagnostics/mace_curvature_20260916/REPORT.md) passes on both checkpoints for the two archived GGR motions in two representations: DFT-anchored displacement errors below0.0045kcal/mol,40MACE+40GB calls,zero newDFT. This validates a narrow directional approximation; negative peptide curvature remains. The coupled scaffold test now passes local curvature and three Ca energy-change checks, but no paired La/Ca response score is supported. Relaxation/entropy remain `response_model_not_validated`; unavailable numerical corrections stay null. [Saved-output scope](../diagnostics/mace_large_20260916/OUTPUT_AUDIT_PLAN.md). |
| MACE hybrid research | `scripts/mace_hybrid.py`; distinct finite-medium, analytic-medium and analytic-large protocols in linked reports | Full 9,141-atom inference works on one A5000: analytic medium ~58 s/endpoint, 9.55 GiB allocated GPU; large ~121 s, 14.98 GiB allocated / 20.89 GiB reserved. Analytic dipoles fix the tested rotation defect; both checkpoints pass core/full rotation and charge checks. Hybrid partition shifts remain −5.023/−4.511 kcal/mol, failing the frozen 2 target. Large full raw contrast differs by +671.31 kcal/mol from medium, primarily in the recorded electrostatic component, with broad force/charge changes. No calibrated class or accuracy gain claimed; baseline unchanged. All 12 large calls completed (1200525). Read-only charge tracing (1200676/1200677) reproduces both models; no near-singular normalization denominator found. A uniform protein-H bond repair completed (1200681/1200682); global model disagreement increases to 923.917 kcal/mol, so the response problem persists. Frozen-monopole OBC-II solvent checks now pass (1200700): medium/large raw-contrast disagreement drops from 923.917 to 106.383 kcal/mol. Warm full-protein solvent evaluation takes about 0.09 s; no calibrated decision or combined gradient. [Solvent report](../diagnostics/mace_gb_20260916/REPORT.md). The five-structure whole-protein panel completed (24 MACE + 24 GB calls): both solvent models reverse all three predeclared relative-order tests. Numerical checks pass, but direct predictive utility fails; 1,370 allocated GPU seconds total. [Panel result](../diagnostics/mace_global_benchmark_20260916/REPORT.md). Local decomposition completed (20 MACE + 20 GB calls): the local PQQ ordering is correct and the full-system contribution reverses it; alpha is already misordered locally and worsens globally. [Decomposition](../diagnostics/mace_local_correction_20260916/REPORT.md). The saved short-range component gets PQQ right but misses both alpha/GGR comparisons in both checkpoints; no fitted rescue or calibrated decision. [Short component](../diagnostics/mace_short_range_20260916/REPORT.md). A separate [curvature pilot](../diagnostics/mace_curvature_20260916/REPORT.md) passes its declared local deformation-energy gates in both checkpoints; no production promotion. [Charge traces](../diagnostics/mace_response_trace_20260916/REPORT.md), [H preparation](../diagnostics/mace_hydrogen_20260916/REPORT.md). [Large results](../diagnostics/mace_large_20260916/REPORT.md), [runbook](../diagnostics/mace_large_20260916/RUNBOOK.md), [analytic medium](../diagnostics/mace_analytic_20260916/REPORT.md), [rotation audit](../diagnostics/mace_rotation_20260916/REPORT.md), [memory implementation](../diagnostics/mace_hybrid_20260916/MEMORY_RESULTS.md). |


The exact short-component engine now passes50 actual energy/rotation/physical
force checks (34calls, job1200736). It omits the global field blocks while
preserving their preceding learned local scalar exactly. This provides a tested
mechanical building block, not a successful standalone classifier; the earlier
alpha/GGR failures remain. [Engine report](../diagnostics/mace_short_engine_20260916/REPORT.md).

The [coupled response experiment](../diagnostics/mace_mechanics_20260916/REPORT.md)
is complete: local curvature, grid, analytic-gradient and SCF-bridge checks pass
for GGR in two representations and both alpha structures. Three eligible Ca
predictions agree with independent DFT+scaffold energy changes within0.001kcal/mol.
Every La endpoint is unstable in this quadratic model or outside its frozen
trust region, so all paired corrections remain null.39newDFT,116MACE-core,
112short-component and116GB calls cost2,560GPU-s and197,248allocatedcore-s.
This validates a narrow mechanical construction, not a discriminator improvement.
The [runbook](../diagnostics/mace_mechanics_20260916/COMMANDS.md) supplies exact
replay/collection commands; no mechanics job remains live. Production is unchanged.

The [canonical direct MACE trial](../diagnostics/mace_canonical_20260916/REPORT.md)
is complete and rejected: calibration classes overlap by 5.817 kcal/mol for
primary medium and 1.609 for large. All 25 calibration and three transfer scores
were computed, but no bands or transfer classifications can be released under
the frozen rule. Both crystal cases share calibration accessions. This is
retrospective PQQ functional-class evidence, not broad affinity validation.
108 new MACE + 108 GB calls, four endpoint reuses per method, zero DFT;
1,212 GPU-seconds and 19,392 allocated core-seconds. No canonical job remains
live. The [runbook](../diagnostics/mace_canonical_20260916/COMMANDS.md) and pinned
receipts support replay. No old reference, threshold or response term is inherited.

A subsequent [saved-output fixed-field screen](../diagnostics/mace_fixed_field_20260916/REPORT.md)
combines vacuum DFT + uniform CHELPG coupling + full-minus-core MACE local energy.
It reduces the consumed Asp303 partition jump from 10.24 to 6.60/5.56 kcal/mol
(medium/large), still failing the 2 target. All 38 accounting checks pass; zero
new quantum/model/solvent calls. This inherited sidechain test is not a complete
peptide-residue test. No solvent or affinity score is available and no expansion
of this frozen mixture is justified. The active goal continues.

The separate [MACE-OMOL trial](../diagnostics/mace_omol_20260917/REPORT.md)
passes all nine numerical checks and canonical PQQ transfer: 25 calibration
cases separate and all three consumed transfer structures classify correctly
with its own vacuum-descriptor bands. Broader affinity robustness fails:
alpha-minus-GGR is negative with the primary extended core and positive with
the connected core; GGR shifts by -26.19 kcal/mol. Do not select the favorable
representation or apply PQQ bands to generic sites. Protocol
`mace_omol_0_100m_vacuum_descriptor_v1` uses native MACE0.3.16, a pinned83-element
checkpoint, explicit charge/multiplicity and no density/solvent correction.
70 actual calls cost670GPU-s/10720allocatedcore-s; median pair inference1.454s.
Jobs1200797/1200799 are complete. [Commands](../diagnostics/mace_omol_20260917/COMMANDS.md).
Baseline remains default and broad-affinity validation remains incomplete.

The [GGR OMOL readout replay](../diagnostics/mace_omol_20260917/READOUT_REPORT.md)
is complete: four unchanged endpoints reproduce all energies/forces, with seven
replay/accounting checks passing. A geometry-independent linear charge/spin
readout contributes -36.188 kcal/mol to the -26.192 representation shift; the
remaining terms contribute +9.996. This identifies an extensive model term,
not an independently validated score correction. Cost: 50 GPU-seconds,
800 allocated core-seconds; zero DFT. A separate
[matched coordination candidate](../diagnostics/mace_omol_20260917/COORDINATION_PLAN.md)
compares bound and separated-metal geometries with identical atoms and charge.
Qualification1200803 and benchmark1200804 completed: all19+60 numerical checks
pass, but the [candidate fails its predictive gates](../diagnostics/mace_omol_20260917/COORDINATION_REPORT.md).
Calibration25/25 separates;1H4I/4MAE transfer correctly,1KB0 is inconclusive.
Both primary alpha-minus-GGRextended directions remain wrong; only1/4 total
non-PQQ comparisons passes. GGR sensitivity shrinks to-13.565kcal/mol but remains
large. Do not promote this refinement or inherit another scorer's bands.
70newcalls,64boundendpoint reuses,zeroDFT cost695GPU-s/11120allocatedcore-s;
median summed four-endpoint inference2.881s,excluding startup. The original
OMOL PQQ research candidate and production baseline remain unchanged.

The [intact-chain OMOL candidate](../diagnostics/mace_omol_20260917/INTACT_REPORT.md)
passes all three frozen development orderings: XoxF−MxaF +9.03 kcal/mol and
alpha−GGR +37.15/+28.09 on the two alpha structures. The primary core model
failed both alpha/GGR comparisons. This is one consumed biological affinity
comparison plus a separate PQQ functional-class comparison, with no absolute
bands or broad validation. Production remains unchanged.
Exact edge batching and allocator recovery fit the 9088-atom protein on A5000:
about 21.72 GB peak and 26 seconds per endpoint. Jobs 1200815/1200816 complete;
their 16 new successes and one OOM cost 1387 GPU-seconds / 22,192 allocated
core-seconds. Four qualified alpha endpoints are reused. Native CPU equivalence
and rigid/repeat checks pass. Reporting overhead remains substantial.
The [product batching qualification](../diagnostics/mace_omol_20260917/EXACT_PRODUCT_PLAN.md)
passes its core and full stages; all24 native CPU/product GPU comparisons agree
exactly. The [28-reference extension](../diagnostics/mace_omol_20260917/INTACT_CANONICAL_PLAN.md)
completed as1200819 with104 new endpoints and8 qualified reuses, but its
[calibration fails](../diagnostics/mace_omol_20260917/INTACT_CANONICAL_REPORT.md):
min(La)−max(Ca)=−306.006296kcal/mol. No new classification bands were issued.
Four reference charges exceed the training range and remain explicitly flagged.
An independent [connectivity audit](../diagnostics/mace_omol_20260917/INTACT_CANONICAL_INTEGRITY_ADDENDUM.md)
found two false peptide bonds across missing 1KB0 structure before its endpoints
ran. Its whole-chain score must remain unavailable, with the case retained in
the denominator. The new preparation policy rejects it before template matching;
27 other cases pass. Use the mandatory audited v2 reporter in the commands
below; the original reporter frozen with the run is superseded. Fixed-core
baseline1KB0 remains unchanged. Eight real-fixture tests pass.
H2001200809 remains queued. [Commands](../diagnostics/mace_omol_20260917/INTACT_COMMANDS.md).

The [saved-output locality audit](../diagnostics/mace_omol_20260917/INTACT_LOCALITY_REPORT.md)
reproduces all five scores and three contrasts; geometry-independent terms
cancel exactly and changes beyond18A are below4e-11kcal/mol per endpoint.
This supports intact local context, not full long-range electrostatics. Global
charge still conditions local nonlinear features. The separately declared
[disconnected sodium check](../diagnostics/mace_omol_20260917/INTACT_SPECTATOR_PLAN.md)
tested that ambiguity on20 new endpoints without changing the existing scores.
It [fails](../diagnostics/mace_omol_20260917/INTACT_SPECTATOR_REPORT.md): four
scores shift4–15kcal/mol and XoxF/MxaF reverses. The actual checkpoint has
identical charge features for−100..−6. This whole-chain version is rejected
for promotion; baseline and the original five-case observation remain intact.
The next [charge-feature ablation](../diagnostics/mace_omol_20260917/CHARGE_ABLATION_PLAN.md)
is a separately declared learned descriptor, with no molecular energy claim.
Recover the [active checkpoint](../diagnostics/mace_discriminator_goal_20260916/CURRENT.md)
before resuming implementation or execution.

Jacob reconfirmed blanket analysis/resource authorization on2026-09-17 and
explicitly requested removal of the old per-analysis check-in rule. Both project
AGENTS.md and `/home/jwestrob/.codex/AGENTS.md` now record this. Continue contained
experiments autonomously toward the active goal; baseline protection remains.

“Baseline” can refer to the electronic method or to a specific preparation.
Always name both. The baseline method is ORCA 6.1.1 native r2SCAN-3c,
CPCM(Water), DefGrid3, NoAutostart, using its native basis/ECP and composite
corrections. There is no custom La basis override. The new global experiment
deliberately uses vacuum endpoints under its separate protocol.

The historical point-charge embedding, CPCM/PB transfer challenger, and
whole-protein GFN2 attempt did not become production methods. Their archived
failures do not disprove global environmental models as a class. The completed
global v1 pilot also failed its own physical gates.

## 2. Interpret scores on the correct scale

With endpoint energies in Hartree:

```text
R_Ha = E_Ca - E_La
R_kcal = R_Ha * 627.509474
S_kcal = (R_Ha - A_Ha) * 627.509474
```

Convert once. Larger values are more La-like on the named protocol's scale.
A positive value alone is not a general affinity, occupancy, functional-use,
or biological “Ln-evolved” classification.

The exact fixed-core PQQ release defines Ca-supported S <=
**14.857129202922806**, La-supported S >= **23.460061205609236** kcal/mol,
and an indeterminate open interval between them. Read the authoritative
numbers from [result.json](../diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json),
`calibration.released_supported_bands`; do not refit or round before classifying.
Its primary released contrast is R; S uses the recorded additive reporting gauge.
That gauge's `registry_compatibility_claimed` is false: it is not a blanket
reference-registry compatibility claim for new preparations.

The finalized symmetric CN8 record is
`reference_inputs/aquo_cn8_native_r2scan3c_v2/aquo_reference.json`.
Its reporting value is A = -646.0775458314704 Hartree. Use it only through the
recorded compatible protocol or explicit reporting-gauge policy. Repaired v3
benchmarks use a shared gauge without inheriting the old zero or PQQ bands.
Matched differences between sites can be calculated directly from R, cancelling
the common offset. The global challenger has neither an S nor an absolute class.

## 3. What has actually completed

| Evidence | Actual outcome / authoritative record |
|---|---|
| Fixed-core PQQ calibration | 25/25 separated; gap 8.602932 kcal/mol. Motif, charge and composition already separate this panel, so it does not prove incremental DFT information. [Release](../diagnostics/pqq_pmdh_fixed_core_calibration_20260914/RESULT.md). |
| Frozen crystal transfer | 1H4I and 4MAE passed their released bands. [Result](../diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/result/HOLDOUT_RESULT.md). |
| Baseline/repair benchmark | All 26 endpoints complete. Fixed-core 1KB0 passes the Ca band. Original GGR passes its frozen direction test; repaired GGR changes sign without demonstrating improved prediction. [Results](../diagnostics/baseline_benchmark_20260915/RESULTS.md). |
| Hans-LanM / alpha-lactalbumin additions | All ten endpoints complete. Hans ranks above GGR; both alpha source geometries rank below GGR, conflicting with condition-qualified La-favoring evidence. [Results](../diagnostics/benchmark_set_20260915/SCORING_RESULTS_1199508.md). |
| GGR mechanism study | All 38 endpoints and 12 directional gradient checks complete. Representation and source geometry materially change scores; no validated mechanical correction. [Report](../diagnostics/ggr_mechanism_plan_20260915/REPORT.md). |
| Global electrostatic pilot | Four vacuum endpoints and four ESP checks complete. Primary partition shift is 18.08 kcal/mol against a 2 kcal/mol limit: failed. All 25 solver checks complete; surface/rotation tests also fail. The ten-endpoint accuracy stage did not run. The separate density/field diagnostic is complete; see its row above. [Final experiment report](../diagnostics/global_electrostatic_20260916/REPORT.md). |

These results are already consumed for development. Multiple chains, structures,
homologs, mutants and sites are not automatically independent observations.
Aequorin and Hans sites remain ordered vectors; parvalbumin evidence remains
supporting/cross-study. Keep canonical PQQ class transfer separate from direct
La/Ca affinity evidence. Do not pool all rows into one accuracy percentage.

The scored benchmark ledger is
`workspaces/benchmark_set_20260915/scored_release_1199508/benchmark_manifest.json`.
The GGR study exports a later derived ledger under
`workspaces/ggr_mechanism_20260915/report_v1/`; consult its report for the
named artifact. Construction-time `prepared_unscored` records are historical
snapshots, not the latest execution status. The newer
[challenge curation](../diagnostics/accuracy_strategy_20260915/CHALLENGE_PANEL_CURATION.md)
records additional evidence/construct constraints and proposals, not new scores.

## 4. Environment and safe read-only entry points

On biotite:

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ALQUEMIA_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ALQUEMIA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

git status --short
tail -n 80 SESSIONS.md
squeue -u jwestrob

"$ALQUEMIA_PY" scripts/affordable_workflow.py dry-run \
  --manifest workspaces/baseline_benchmark_20260915/run_v1/manifest.json
"$ALQUEMIA_PY" scripts/affordable_workflow.py dry-run \
  --manifest workspaces/benchmark_set_20260915/ready_tasks_v4/manifest.json
```

ORCA is pinned at
`/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca`.
Use the manifest's executable and implementation pins, rather than an arbitrary
PATH executable. The older PDBFixer entry point is
`/home/jwestrob/miniconda3/envs/fep/bin/python`; exact benchmark reproduction
should replay the archived protonated coordinates instead of regenerating H atoms.

Recollect completed baseline results without rerunning quantum calculations:

```bash
"$ALQUEMIA_PY" scripts/affordable_benchmark.py collect \
  --manifest workspaces/benchmark_set_20260915/ready_tasks_v4/manifest.json \
  --output workspaces/benchmark_set_20260915/collection_agent_recheck_01.json
"$ALQUEMIA_PY" scripts/affordable_benchmark.py report \
  --collection workspaces/benchmark_set_20260915/collection_agent_recheck_01.json \
  --output workspaces/benchmark_set_20260915/report_agent_recheck_01.md
```

Writers reject existing output paths. Preserve prior versions; choose a new
explicit output filename for another recollection. Collectors verify normal
termination, convergence, input/coordinate/receipt hashes and protocol grouping.
Missing or invalid results remain visible and unscored.

## 5. Preparation and execution

For an approved generic benchmark preparation, this real example replays the
exact alpha-lactalbumin source preparation and writes a fresh workspace:

```bash
"$ALQUEMIA_PY" scripts/affordable_benchmark_set.py prepare-generic \
  --root "$ALQUEMIA_ROOT" \
  --config diagnostics/benchmark_set_20260915/preparation_configs/1F6S.json \
  --protonation-report workspaces/benchmark_set_20260915/prepared/alacta_1f6s_v1/preparation_report.json \
  --output workspaces/benchmark_set_20260915/prepared/alacta_1f6s_agent_replay_01
```

This prepares inputs; it submits nothing. For other targets, use an explicit
reviewed configuration, exact metal/chain/site and source-state inventory.
[Benchmark commands](../diagnostics/benchmark_set_20260915/COMMANDS.md) document
task packaging. The generic repair's lower-level CLI requires the original
v2 manifest **and** a pinned molecular topology; residue-number arithmetic is
not a replacement for connectivity. Unsupported chemistry must remain explicit.

The existing inbox driver is an automatic **submission** path, not a dry-run.
It resolves one exact metal site, validates PQQ routing before and after
protonation, then runs the existing v2 carver. Fold-catalogue qualification is
CN >= 6 at 3.1 A; automatic ORCA/carver qualification is CN >= 7, with <= 2
direct N donors. Its separate 3.3-A inclusion buffer does not widen that gate.
It does not invoke fixed-core PQQ or peptide v3 automatically. Do not feed new
work into `inbox/`, run `process_inbox.sh`, restart watchers, or resubmit a
completed historical batch as a documentation verification step.

Approved manifested quantum jobs use `affordable_workflow.py execute` and the
existing `run_orca_task_manifest.py` / `render_orca_runtime_input.py` machinery.
The recent 64-CPU batches use four concurrent 16-rank ORCA endpoints, one
OpenMP thread per rank. Runtime `%pal` insertion and MPI launch fixes are
intentional; the old “NEVER MPI” instruction is obsolete. Retain the recorded
runtime inputs and receipts. Valid completions are cached; partial quantum
attempts require an explicit fresh retry directory. Never remove another
executor's lock or overwrite an attempt.

Use the batch script associated with the approved manifest. Jobs1199299 and
1199508 are complete; their launch commands are historical/recovery recipes,
not pending tasks. The general inbox produces auxiliary apo/water inputs;
the recent benchmark manifests use only the required La/Ca endpoint pair.

## 6. Current research jobs and handoff discipline

The current global experiment has its own [agreement](../diagnostics/global_electrostatic_20260916/AGREEMENT.md),
[frozen settings](../diagnostics/global_electrostatic_20260916/NUMERICS.md),
[report](../diagnostics/global_electrostatic_20260916/REPORT.md), and
[exact runbook](../diagnostics/global_electrostatic_20260916/RUNBOOK.md).
Read them before touching its tasks. Collect the original 25-task master
manifest even though execution is split across jobs. A subset success cannot
trigger Stage 2. Do not use an old baseline cache entry as a challenger result.

Candidate products belong under `workspaces/`; compact plans/inventories/results
belong under `diagnostics/`. Preserve unrelated dirty changes and immutable
experiments. Record scientific scope and consequential settings in project notes.
For the active MACE discriminator goal, Jacob explicitly authorized autonomous
analyses and method development and removed the per-analysis approval gate;
see [AGENTS.md](../AGENTS.md) and the durable goal authorization. Do not restore
that gate after compaction. The user has
removed project CPU/time stopping budgets; record actual costs and task counts
without inventing new spending limits. Existing scheduler policies still apply.

For cost reports, distinguish allocated Slurm CPU-seconds, measured CPU use,
rank-times, wall time and physical-core topology. Shared allocations can expose
SMT siblings; CPU utilization alone does not establish clock rate or throughput.
Current [resource observations](../diagnostics/global_electrostatic_20260916/PERFORMANCE_OBSERVATION.md)
document one same-input recovery, with all interrupted cost retained.

Update SESSIONS.md and the relevant experiment report at handoff. Commit only
your own scoped files/hunks; do not push, change defaults, or launch a production
rescore through this guide. The [August campaign record](../CANONICAL_METHOD_OPERATIONS_RESULTS_2026-08-30.md)
and older reports remain provenance for those campaigns.


### Charge-feature ablation development result (2026-09-17)

The separate whole-chain MACE descriptor passes its five-case relative-order
and numerical tests after masking only the raw global charge feature. All
physical inputs and learned parameters are retained. Its outputs are explicitly
model units, not quantum energies or binding free energies. Native whole-chain
canonical/spectator failures remain recorded. No default changes.

Job1200828 completed42calls in798GPU-allocation seconds. The conditional
25-reference calibration is submitted as1200830 (100newcalls,8exactreuses);
1KB0 remains unsupported. Recover [CURRENT.md](../diagnostics/mace_discriminator_goal_20260916/CURRENT.md)
and use the [frozen commands](../diagnostics/mace_omol_20260917/CHARGE_ABLATION_COMMANDS.md).
Do not issue candidate bands before the declared full-panel criterion passes.


### Masked MACE canonical milestone (2026-09-17)

Job1200830 completed100/100newforwards on oneA5000 in59m06s. The masked
whole-chain descriptor separates all25designated PQQ calibration proteins
(gap9.108591modelkcal), and correctly transfers to1H4I/4MAE.1KB0 remains
unsupported; retain the full3-transfer denominator and false all-three gate.
These consumed classes do not establish broad affinity or incremental value
beyond composition. [Canonical result](../diagnostics/mace_omol_20260917/CHARGE_ABLATION_CANONICAL_REPORT.md).

All189panel factorization checks pass. Optional two-call research calibration
is implemented in mace_omol_prepared.py report --calibration, with strict
numeric/model/preparation compatibility and no PQQ-band transfer to non-PQQ
sites. Use the existing OpenMM driver; workers use their separately pinned
MACE environment. [Calibration](../diagnostics/mace_omol_20260917/MASKED_CALIBRATION_REPORT.md),
[commands](../diagnostics/mace_omol_20260917/PREPARED_INPUT_COMMANDS.md).
Production/default and all old references remain unchanged. Recover CURRENT.md
for actual integration test status and outstanding research.


### Expanded GGR robustness and multisite development (2026-09-17)

The masked model's five-case development pass does not survive all existing
GGR structures: only2/6alpha-minus-GGR structural comparisons pass. All four
new forwards passed numerical accounting; this is a predictive robustness
failure. The PQQ calibration remains unchanged.
[GGR result](../diagnostics/mace_omol_20260917/GGR_STRUCTURE_ROBUSTNESS_REPORT.md).

New source-backed preparation paths check raw chemistry before filtering.
Multisite support retains other calcium ions, the fixed union of site waters
and actual N-terminal acetyl chemistry where present. It is a separate policy;
old single-metal inputs are unchanged. Recover CURRENT.md for the five-site
4CPV/1SL8 experiment and live jobs; do not borrow PQQ bands for these proteins.

The five-site multisite trial is now complete:10/10calls and all-Ca permutation
checks pass; the parvalbumin supporting direction gate fails4/6contrasts.
Aequorin remains an ordered vector with unresolved site labels.
[Multisite report](../diagnostics/mace_omol_20260917/MULTISITE_REPORT.md).
Exact masked-descriptor gradients now work on intact GGR:23/23 core checks and
11/11 whole-protein checks pass. The source-mapped4698-atom gradient costs about
41seconds per endpoint and12.0GB peak GPU memory on oneA5000. A separately
versioned checkpoint/scatter adapter preserves the scalar; production energy-only
adapters remain unchanged. This is numerical derivative qualification, not a
physical relaxation correction. [Gradient report and commands](../diagnostics/mace_omol_20260917/MASKED_GRADIENT_REPORT.md).
The completed [response screen](../diagnostics/mace_omol_20260917/MASKED_RESPONSE_REPORT.md)
used40 real-core evaluations and zero new DFT. Numerical derivatives pass24/24;
direct responses pass6/24 and curvature19/24. No mechanical correction is
supported. The separate static DFT+masked(full-core) candidate fails its cached
partition test by20.16kcal versus the2kcal target; its six conditional full calls
must not be launched. [Partition result](../diagnostics/mace_omol_20260917/MASKED_SUBTRACTIVE_CONTEXT_REPORT.md).
The shared learned-neutral feature now passes16/16 numerical checks on eight
real cores, but its subtractive partition shift is10.91kcal versus the2kcal limit.
Its conditional six whole calculations must not be launched. No charge-category
sweep or production change. [Result and replay](../diagnostics/mace_omol_20260917/SHARED_NEUTRAL_FEATURE_REPORT.md).


The [matched vacuum diagnostic](../diagnostics/mace_omol_20260917/MATCHED_VACUUM_REPORT.md)
completed eight native-ORCA centers. Solvent mismatch accounts for18.33 of the
20.16kcal GGR discrepancy; vacuum DFT minus masked-core residual is1.83kcal,
within the declared2 target. Shared-neutral residual still fails. Four tests
pass. No aqueous score or baseline change. The separate six-call whole-protein
[vacuum hybrid](../diagnostics/mace_omol_20260917/VACUUM_HYBRID_REPORT.md) now
passes11/11 numerical and partition checks, but only1/4 ordering comparisons.
It has not earned an accuracy claim. The [matched H-preparation result](../diagnostics/mace_omol_20260917/MATCHED_H_REPORT.md)
now passes13/13 numerical checks and reduces GGR partition sensitivity to0.1185kcal.
Eight DFT and eight small MACE cores completed; six whole results reused.
Only2/4 ordering comparisons pass: both1F6S comparisons succeed, both6IP9 fail.
The accuracy gate remains false. Five new tests plus six native regressions pass.
Heavy atoms/caps/chemical inventories stay fixed; no aqueous score or relaxation.
[Commands](../diagnostics/mace_omol_20260917/MATCHED_H_COMMANDS.md). Prior failures remain immutable.
[Vacuum replay commands](../diagnostics/mace_omol_20260917/MATCHED_VACUUM_COMMANDS.md).

The [normalized charge check](../diagnostics/mace_omol_20260917/NORMALIZED_CHARGE_REPORT.md)
completed eight native CHELPG/eight exact-potential utilities. Endpoint and paired
quality gates pass; projected Ca/La field errors remain below0.6%relativeRMS.
This does not establish solvent-energy accuracy. Four real-fixture tests pass.
The [full-boundary solvent plan](../diagnostics/mace_omol_20260917/FULL_BOUNDARY_GB_PLAN.md)
is declared but not yet implemented or submitted. It preserves all physical
atoms and uses a local charge ledger, no coreCPCM or added bare Coulomb term.

The [full-boundary GB test](../diagnostics/mace_omol_20260917/FULL_BOUNDARY_GB_REPORT.md)
completed48calls:46numerical checks pass, partition4.48kcal fails2 and all four
ordering contrasts reverse. Saved-state coupling shows a large direct/solvent
cancellation that the entangled OMOL context cannot be assumed to reproduce.
Do not promote or add arbitrary bare Coulomb to that model. A distinct
[explicit-field/short-readout plan](../diagnostics/mace_omol_20260917/EXPLICIT_FIELD_SHORT_PLAN.md)
is declared, not yet implemented:8short-core+8nativepotential calls, existing
DFT/GB/whole-short reuse. Earlier short-component failures remain recorded.

## Latest MACE hybrid development: explicit field result

The [explicit field + short MACE candidate](../diagnostics/mace_omol_20260917/EXPLICIT_FIELD_SHORT_REPORT.md)
completes with all near-boundary charge-representation checks passing, but
fails the partition gate (-4.139kcal) and passes only1/4 consumed alpha/GGR
ordering comparisons. Baseline remains default; no class/reference or combined
gradient exists for this candidate. Eightshort/eightpotential calls cost63GPU-s
and1240allocatedcore-s, reusing quantum/whole/GB sources.
The next [responsive-field plan](../diagnostics/mace_omol_20260917/RESPONSIVE_FIELD_PLAN.md)
is declared, not yet implemented or run. Recover CURRENT.md before launch.

The [responsive-field follow-up](../diagnostics/mace_omol_20260917/RESPONSIVE_FIELD_REPORT.md)
now completes8nativeDFT,16charge/potential utilities and44GB solves. All component
checks pass; updated density/solvent improves ordering to3/4, but GGR boundary
shift-4.438kcal still fails.21distinct tests pass. No reference/class/combined
gradient or production change. Cost45529core-s/66GPU-s includes recovery.
Next is the declared [AMOEBA capability/accounting investigation](../diagnostics/mace_omol_20260917/AMOEBA_CAPABILITY_PLAN.md),
with3real parameterization preps and no energy/force evaluation in that stage.


## Research component update — 2026-09-17

Native coupled polarization now accepts supplied per-site fields and reproduces
its original energies/dipoles on all six real controls. Energy accounting and
rigid checks pass; this is still an ion-excluded framework test, with no new
metal prediction or change to the production default. See the
[native field report](../diagnostics/mace_omol_20260917/NATIVE_FIELD_INPUT_REPORT.md)
and [current goal checkpoint](../diagnostics/mace_discriminator_goal_20260916/CURRENT.md).
Compact charge/dipole source fits failed their full spatial field screens;
they were not promoted. The next declared step supplies exact density coupling
to permanent protein multipoles before a full hybrid score.

### Same-density charge sampling completed

Forty native NoIter property replays and40density queries pass all electronic
identity/ESP-quality checks. All8finer-to-finest source-self sensitivity checks
pass, but5/8fitting-extent checks fail the frozen0.5kcal criterion. The GGR
source-self discrepancy remains14.97–18.13kcal across settings: sampling matters
but does not explain most of it. No default change or full hybrid rescore.
[Result and costs](../diagnostics/mace_omol_20260917/CHELPG_SAMPLING_STABILITY_REPORT.md).

An isolated ddX0.9.0 build imports successfully. The subsequent real solver
pilot is described below. Exact-density coupling needs
both surface potentials and a density integral; see the
[interface assessment](../diagnostics/mace_omol_20260917/DDX_CAPABILITY_NOTE.md).


### Resolved protein-boundary solver development, 2026-09-18

The isolated ddX backend now executes real source-only protein PCM calculations.
The initial native-source FMM precision gate failed; a separately versioned
exact Coulomb source evaluation preserves the physical inputs and passes that
check. The full grid/rotation pilot1201279 is running. One unchanged La state
needed337 dielectric iterations, exceeding the original300; logged replay1201280
converges at the same1e-10 tolerance in264.183solve seconds. Numerical credibility
across the inventory and predictive value remain unestablished.

The Python native Model retains an error after nonconvergence. The new recovery
runner creates fresh Model objects, preserves failed attempts, and reuses actual
successful coefficients. [Recovery scope](../diagnostics/mace_omol_20260917/DDX_ITERATION_RECOVERY_PLAN.md).
This is a frozen projected-charge solvent-component diagnostic: no complete
hybrid rescore, exact-density claim or production/default change. Recover
CURRENT.md for current job and manifest status before executing anything.


### Affordable conductor component,2026-09-18

The distinct conductor-like protein-boundary pilot completed20/20real solves
in285allocated wall seconds(64CPUs). It reduces the same-charge GGR source-self
structural discrepancy from18.129kcal(GK) to1.89–2.10across its three resolutions.
Its frozen numerical gate still fails420/430checks; refinement/rotation errors
remain. This is model sensitivity, not demonstrated biological improvement.
[Actual result](../diagnostics/mace_omol_20260917/DDX_CPCM_REPORT.md).

Higher-resolution job1201302 retains all physical parameters and tolerances:
18/974versus24/2030,16newroles plus4exactreuses. Native energies, conductor
prefactor and unrounded scaled values remain separate. No baseline/default,
aquo reference, threshold, fullhybrid score or combinedgradient change.
[Scope](../diagnostics/mace_omol_20260917/DDX_CPCM_REFINEMENT_PLAN.md).
The separate GMRES qualification converges but fails strict equivalence checks
and is slower; its cost split motivated the conductor approximation, not a new
class label. [GMRES result](../diagnostics/mace_omol_20260917/DDX_KRYLOV_REPORT.md).


### Polarization accounting and conductor refinement, 2026-09-18

The full frozen-response reaction-energy expression passes 28 native checks
on four real GGR endpoints. It retains both AMOEBA dipole states and all
permanent/induced cross terms. A corrected isolated wrapper restores native
cutoff initialization; the initial failed execution remains recorded.
[Derivation and result](../diagnostics/mace_omol_20260917/FROZEN_RESPONSE_REPORT.md).

Conductor refinement completes but fails six of 450 numerical checks.
The between-structure contrast stabilizes, while individual energies still
need refinement. Job 1201312 separates angular-basis and integration errors;
[scope](../diagnostics/mace_omol_20260917/DDX_RESOLUTION_PLAN.md).
Finite-dielectric PCM recovery has finished unsuccessfully; no further PCM
campaign is running. Recover CURRENT.md before execution. No new hybrid score,
calibration, baseline integration or production/default change is available.
