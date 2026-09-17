# Active MACE discriminator checkpoint — 2026-09-17

**Goal active; candidate milestone reached.** Jacob authorized autonomous
analyses and available resources. Home/project AGENTS now remove the old
per-analysis permission gate. Preserve baseline/default, immutable studies,
other agents' edits/jobs. No push/promotion. No project CPU/time/token budget.
Recover GOAL.md and inspect live jobs before resuming.

## Latest execution checkpoint — intact gradients qualified

This supersedes the earlier active-gradient notes. Full job1200886 COMPLETED;
all six calls and11/11 numerical checks pass. WholeGGR has4698 mapped atoms;
Ca/La gradient evaluations take40.699184/40.789196seconds and12,011,144,704bytes
peak GPU (11.19GiB). Both center energies equal the archived scalar. Along the
fixed metal-to-GLN140/O direction, grad(R)=+8.836946204modelkcal/A; actual odd
0.01A change0.088339618 versus predicted0.088369462. No physical-force,
relaxation, entropy or biological-improvement claim. Production unchanged.

Final adapter `omol_exact_checkpointed_edge_product_autograd_v3` requalified
on the exact73-atom1H4I core (job1200885,23/23checks). The actual fullV2 CUDA
OOM was resolved by replacing index_add with equivalent scatter_add: saved
edge-message tensors no longer consume GPU memory outside checkpoints.
Core/native max gradient error9.02e-14eV/A. Earlier coordinate/report roundoff,
TorchScript early-stop, collector-serialization and OOM attempts remain visible.
No coordinates, scientific tolerances or criteria changed during recovery.

Authoritative artifacts under workspaces/mace_omol_20260917:
- core: masked_gradient_core_v4/manifest.json, core_report_v5 is named in full:
  `masked_gradient_core_report_v5/result.json` (all checks pass).
- whole: `masked_gradient_full_v3/manifest.json`,
  `masked_gradient_full_report_v1/result.json`, plus all-atom TSV/paired array.
- costs: `masked_gradient_cost_v1.json`, seven terminal jobs1200863/64/65,
  1200878/84/85/86.27successful calls,2failed model calls,410GPUallocation-s,
  6560allocatedcore-s,470.555reportedCPU-s. Two preflight failures made0calls.
- cumulative: `intact_engineering_status_v12.json`:389successes,4failedcalls,
  12022GPU-s,255584allocatedcore-s,30231.464reportedCPU-s.
- seven final real-fixture tests pass:6in12.883s plus mapped-export test1.531s;
  none skipped. Earlier native regressions and all resource receipts retained.

[Gradient report](../mace_omol_20260917/MASKED_GRADIENT_REPORT.md) includes replay
commands and immutable attempt details. Source files: mace_omol_gradients.py,
mace_omol_gradient_worker.py, mace_omol_gradient_run.py; minimal dispatch in
mace_omol.py. No production score/energy-only adapter changes.

## Latest scientific results after the gradient milestone

Response job1200888 COMPLETED40/40calls, zero failures. Exact archived DFT inputs
and physical Jacobians were used. Own analytic derivatives pass24/24; direct
DFT deformation-energy comparisons pass6/24; curvature passes19/24; DFT-gradient
anchored predictions pass24/24 under their separate0.02 absolute floor. The
anchored pass does not override five curvature failures. Maximum direct error
0.320106; anchored0.019487kcal-scale. No relaxed score is supported. Report:
`masked_response_report_v1/result.json`; compact MASKED_RESPONSE_REPORT.md.

Cost320GPUallocation-s,5120allocatedcore-s,416.572reportedCPU-s;20.148532summed
model-s,1,031,137,792bytespeakGPU. Three distinct real-fixture tests pass; actual
integration was explicitly skipped before outputs then passed after completion.
CumulativeengineeringV13:429successful calls,4failedcalls,12342GPU-s,
260704allocatedcore-s,30648.036reportedCPU-s. No newDFT/solver/training.

A separate static candidate was then declared: DFT/CPCMcore + masked(full-core).
Its cached partition prerequisite FAILS: connected-minus-extended GGR shifts
DFT−7.343500873, maskedcore−27.503809106, hybrid+20.160308232 versus2kcal gate.
All original-H shared source/cap mappings pass (max4.99e-11A). Zero new inference;
DO NOT run its conditional six whole-chain calls. Full result:
`masked_context_partition_v1/result.json`; compact
MASKED_SUBTRACTIVE_CONTEXT_REPORT.md. Two actual-fixture tests pass4.384s.
No adaptive change, learned multiplier or omitted failure. Baseline unchanged.

**Next declared development:** [SHARED_NEUTRAL_FEATURE_PLAN.md](../mace_omol_20260917/SHARED_NEUTRAL_FEATURE_PLAN.md).
One fixed alternative to the raw-zero vector: use the checkpoint's learned
charge-zero feature for both endpoints/all atoms, preserving actual physical
charges and spin. This remains an empirical descriptor, not physical
neutralization. No charge-category sweep. First eight real core-center
energy/gradient calls on the same GGR/alpha set; verify four actual La charge0
states against exact native archived results where compatible. Require native
numerical agreement and the same2kcal cached partition target before any
conditional six whole-system calls. No optimization or newDFT. Plan declared;
new adapter/reference audit/manifest NOT IMPLEMENTED OR LAUNCHED yet.

Source scripts for completed work: mace_omol_response.py, mace_omol_context.py,
minimal mace_omol.py dispatch, real-fixture tests. Gradient milestone committed
c22a881; response/context work is the next scoped commit. Do not repeat passed
qualifications or rerun the failed fixed context candidate.

All new research GPU jobs are terminal. Native H2001200809 remains pending at
last check; preserve other agents'PLM1200794–1200796 and inspect live state.
Goal active. Full autonomy, baseline/reference protection and no push remain.

## Actual result: masked MACE candidate works on the supported PQQ panel

Protocol `mace_omol_intact_charge_feature_ablation_descriptor_v1`.
The same pinned OMOL checkpoint and whole-chain physical preparation are used,
with only the raw total-charge embedding zeroed before the joint projection.
Actual physical charges, spin and weights remain recorded and unchanged.
This is a learned descriptor, not a quantum electronic-state energy or binding
free energy. It lacks responsive long-range electrostatics; whole proteins are
outside published model training sizes.

- Development1200828:42/42newforwards, all numerical/three-order gates pass.
  XoxF−MxaF44.329179modelkcal; alpha1F6S−GGR12.741688;
  alpha6IP9−GGR4.112398. The two alpha structures are one qualified affinity
  comparison. Sodium score change−1.07e-8 is architectural consistency.
- Canonical1200830: **COMPLETE100/100newforwards**,8exact crystal reuses.
  All25calibration proteins separate, gap9.108591modelkcal. MxaF1H4I17.579956
  is Ca; XoxF4MAE61.909135 is La. **2/3transfers valid and correct;1KB0
  unsupported. The declared all-three transfer gate remains false.** All cases
  consumed; composition already separates the canonical panel. No broad
  affinity accuracy or incremental value beyond composition established.
  [Canonical report](../mace_omol_20260917/CHARGE_ABLATION_CANONICAL_REPORT.md).

Actual canonical cost3546GPU/allocation-s(59m06s)oneA5000,56736allocatedcore-s,
3925reportedCPU-s(seconds precision),2657.578770summed model-evaluation-s,
13751386112bytes peakGPU,1769204KiBpeakRSS. Local work separately measured.
Existing bound pairs measured49.85685–60.772997forward-s(median53.068669),
53.744777–65.495802worker-wall-s(median57.267661), excluding prep/controller/report.

## Two-call reference and runnable research classifier

The exact descriptor factorization is
R_mask=[T_bound,Ca−T_bound,La−(C_Ca−C_La)]*23.06054783061903.
Fixed1H4I native detached-node-plus-embedding terms are Ca−18430.794927644074
and La−850.2720512362149 model eV. These are not aquo/quantum-ion energies.
Factorization development60checks and full-panel189checks pass; full-panel
maxerror2.357415596e-8 vs declared0.01tolerance. Zero additional inference.

Use `scripts/mace_omol_prepared.py`: audit,prepare(optional--factorization),
existing mace_hybrid dry-run/execute/collect,report(optional--calibration).
Source-backed whole-chain preparation required; arbitraryXYZ is unsupported.
Strict peptide connectivity, atom/source replay, charge/spin/paired coordinates,
cofactor/water state and actual endpoint receipts enforced.

IMPORTANT driver for prepared-input operations:
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
(OpenMM). GPU workers use the separately pinned existing MACE venv. Never
install into either to hide a mismatch. MACEvenv works for pure panel reports.

Matching two-call PQQ bands:Ca<=42.346663801320496;La>=51.455254788976355;
otherwiseinconclusive. PQQ functional association only; no non-PQQ affinity band.
Original four-call bands stay immutable and numerically separately typed.
Reference: workspaces/mace_omol_20260917/masked_calibration_v1/reference.json
SHA217a127f7ba8f6cccb59d1c63b2d389e9022449a5feb09f58eba2e3b04cba54d.
Frozen current reporter: masked_calibration_source_v1/implementation/.
[Commands](../mace_omol_20260917/PREPARED_INPUT_COMMANDS.md),
[reference report](../mace_omol_20260917/MASKED_CALIBRATION_REPORT.md).

Executed no-new-inference interface reports:
prepared_interface_pqq_calibrated_report_v1:17.57995558550866,Ca;
prepared_interface_ggr_calibrated_report_v1:23.98742145552895,classificationnull,
outside_canonical_PQQ_calibration_scope. Reference and direct algebra replay.
All five actual-artifact calibration guard tests PASSED in675.322seconds,
none skipped; masked_calibration_validation_v1.json pins tests/log/resources.
Standalone inspected PDF/SVG/PNG figure inmasked_calibration_figure_v1. Earlier4prepared,3factorization,
8ablation/interface and6native regressions passed. No synthetic science fixtures.

## Artifacts and live work

All candidate products under workspaces/mace_omol_20260917/ unless specified.
- charge_ablation_canonical_v1/manifest.json SHA
  2415b4dacd755530dec888017d76afe8895708213a36ef09497c7264e341476f.
- collection_job_1200830.json SHA
  144d7e312716672345358b319823da8cd4b61a02782a849ea1b013c921b660de.
- charge_ablation_canonical_report_v1/result.json SHA
  9c12a55ac09aab4c68f2399a1ff784ebd078d6f6a8ea020b85fccb577d017f25.
- factorization_report_v2/result.json SHA
  35ba8bb53b209f650b8392ff03f0c426d0f95455b39b18c082e25aabb5434fbe.
  V1auxiliary NumPy reduction differences preserved; v2uses math.fsum.
- prepared_interface_pqq_two_call_v2/manifest.json has0newtasks/2realreuses.
- prepared_interface_ggr_two_call_v1/manifest.json has0newtasks/2realreuses.
- prepared_interface_ggr_fresh_v1 contains4UNEXECUTEDtasks. Do not duplicateGGR.
- Engineering statusV8:348successfulforwards,2OOMs,11282GPUallocation-s,
  243744allocatedcore-s,29420.954reportedCPU-s. Includes4invalid native1KB0
  diagnostic calls and prior failures; local work separate.
- Native H200job1200809 remains pending unchanged14-taskmanifest. Inspect live
  state; do not cancel/duplicate/alter its manifest or any executor lock.

Earlier own commits2fcf2f3,287d3d0,30cb9b4,410390a. This checkpoint accompanies
the scoped canonical-calibration implementation/results commit; recover its
hash fromgitlog. Many unrelated dirty files are concurrent/historical.
Never blanket-stage/reset/stash. SESSIONS own append
must be staged separately from another agent's working additions.
Email accepted by local relay20260917T145316Z_masked_canonical_pass. No delivery
confirmation claimed. Vault note now includes the canonical milestone, costs, limitations and source finding.

## Important source finding: 1KB0 is not a simple terminal-cap repair

[Raw-source diagnosis](../mace_omol_20260917/ONE_KB0_RAW_CHEMISTRY.md): deposited
TRO512 is present, but fixed-core-oriented normalization removed it and
whole-chain parsing connectedLEU511directlytoGLU513(4.750089A). A separate
missing574–578loop produces573to579(18.969263A). Raw HEC802 is covalently joined
toCYS604/CYS607 and coordinatesHIS608/MET647. Heme was explicitly excluded from
the old normalized source. Do not fabricate a whole-chain pass by converting
TROtoTRP, omittingheme, guessingironstates, or bridgingthegap. New full-system
support needs actual modified-residue/heme/loop chemistry under a new policy.
Old fixed-core1KB0 remains valid under its own recorded scope and unchanged.

## Earlier failures remain results

Native whole-chain OMOL canonical gap−306.006296; no bands. DisconnectedNa
shifted four scores4–15kcal/mol and reversedXoxF/MxaF. Checkpoint charge features
alias−100..−6and+11..+100; verified behavior, no inferredtraininghistory.
Masked representation is a separate development response, not erased native
history. Read INTACT_CANONICAL/SPECTATOR/LOCALITY reports, chargeembeddingaudit.
Prior POLAR+GB reversed allthree directions; local/hybrid partition and paired
mechanics failed. Relaxation/entropy remain unavailable. Do not rediscover them.

## Next useful work

Calibration tests, scoped implementation/results and vault update are complete. The operational
candidate milestone merits reporting, but goal remainsactive: all-three transfer
support and broad/composition-challenging predictive usefulness remain open.
Existing CHALLENGE_PANEL_CURATION.md records eligible versus unresolved affinity
controls. AQUALYSIN_EVIDENCE_FOLLOWUP.md records read-only follow-up: adding/-char/en
recovers the primary2002PDF. The2019Ca-1assignment remains an inference from
structural analogy, not direct assay-site mapping. No new label/preparation/score.
Do not spend effort merely recounting provenance: prioritize chemically sound
new inputs and tests that can reveal predictive weaknesses. No new user approval
is needed within this goal. No automatic production/default promotion.


## Expanded GGR robustness: completed failure, 2026-09-17

Jobs1200845/1200846 COMPLETE four forwards, all numerical checks pass.
GGR ordered [1GLG,2FW0,2FVY] = [23.987421458,44.884127357,45.975500410]
model kcal. Both new structures reverse both alpha-minus-GGR comparisons:
only2/6margins pass; all-case robustness FALSE. No label/threshold/input rescue.
The separate PQQ result stands; broad affinity improvement is unestablished.
See ../mace_omol_20260917/GGR_STRUCTURE_ROBUSTNESS_REPORT.md.

Actual cost114GPUallocation-s,1824allocatedcore-s,118.469reportedCPU-s,
52.124891summed model-s,6108277248bytespeakGPU. Cumulative engineeringV9:
352successfulforwards,2OOMs,11396GPU-s,245568allocatedcore-s,29539.423CPU-s.
All costs include prior failures; local work separate. Four source/actual-result
regressions pass6.713s,none skipped. No newDFT/solver/training/forces.

Source bridge: scripts/mace_omol_source_prepare.py with explicit raw inventory,
exclusions, source metal/paired-coordinate and peptide checks; original1GLG
physical coordinates exactly reproduced, raw1KB0TRO512 explicitly rejected.
Policy omol_source_backed_chain_A_ff19sb_H_explicit_exclusions_v1.
New actualpreps ggr_source_bridge_{2fw0,2fvy}_v1/source_preparation.json;
source-row/evidence/exclusions JSONs in ggr_structure_sources_v1/.
Source tests and initial failed comment-byte assertion preserved; no scientific
coordinate discrepancy. Full result ggr_structure_report_v1/result.json.
Saved readout diagnosis places changes across metal and several donor residues,
not a unique causal term. Same donor identities; several distance shifts0.1–0.2A.

## Multisite family coverage: complete, supporting gate fails

PlanMULTISITE_PANEL_PLAN.md; resultMULTISITE_REPORT.md/MULTISITE_RESULT.json.
Jobs1200851–1200855 COMPLETE tenforwards, allnumericchecks pass. Identical
all-Caenergies acrossselectedsiteatomorderings (reportederror0).
Orderedparvalbumin[CD,EF]=[27.431183345093064,64.34764681174364]modelkcal.
CDexceedsonly1GLG;EFexceedsallthreeGGR. Supportingcross-studycontrasts4/6pass,
all-casegateFALSE. No sub-kcalaffinitysiteorderclaim. Aequorin[EF1,EF3,EF4]=
[13.04767291865186,33.17334887989371,51.5641437611785]. Siteunresolvedassay,
no bestsiteselection,physicalzero,PQQband orbinaryassayreproductionclaim.

New scripts/mace_omol_multisite.py, multisite_report.py. Policy
omol_intact_multisite_fixed_background_Ca_source_acetyl_ff19sb_v1.
Rawsourceheavyatoms/covalentlinksreplay;4CPVactualACE0cap,H,bondretained.
4CPV1611atoms,1backgroundCa,1water,QCa−3/QLa−2.1SL82866atoms,
2backgroundCa,3waters,QCa−4/QLa−3;its10missingN-terminalresiduesexplicit.
Onlyselectedmetalchanges.Background-Caindicesmustbeexplicit;default
single-metalcheckandallpreviouspreparationpathsremainunchanged.
V1prepfailedbeforeinferenceatoldsingle-metalguard;V2fixedthespecifictechnical
limitation.Noalteredscientificselection.Unknownrawchemistry/gapsstillfail.

Artifactsunderworkspaces/mace_omol_20260917/:
- multisite_recipes_v1/{PARV_4CPV,AEQ_1SL8}.json
- multisite_prepared_v2/CASE/preparation.json (all5prepared)
- multisite_scoring_v1/CASE/manifest.json and actualcollections/receipts
- multisite_reports_v1/CASE/result.json; multisite_comparison_v1/result.json
- multisite_cost_v2.json, multisite_sacct_v1.tsv. CostV1hadoneauxiliarytest
  resourcehashcapturedinflight;V2correctsit,scores/allocationcostunchanged.
- Sevennew/source/resulttestspass14.973s,fourlegacyPQQ/GGRtestspass28.313s,
  onenewactualmultisitereporttestpasses9.923s;12distincttests,noneskipped.

Actualmultisitecost216GPUallocation-s,3456allocatedcore-s,221.486reported
actualCPU-s;63.758972model-s,80.695251workerwall-s;3936891392bytespeakGPU.
Localpreparation/dryrun/report/testreceiptsseparate. CumulativeV11through1200855:
362successfulforwards,2OOMs,11612GPU-s,249024allocatedcore-s,29760.909CPU-s.
NoDFT/solver/training/forces/relaxationcalls. Goalactive;productionunchanged.

## Next declared work: analytic descriptor gradients

MASKED_GRADIENT_PLAN.md declaredbeforeimplementation/inference. Add a separate
exactbackward-capableedge/productadapterwithcheckpointing;currentqualified
energy-onlyadaptersmuststayunchanged. Sameweights,mask,spin,float64andalgebra.
StageA10forwards onreal73-atom1H4Icores from ablationdevelopmentV2:2nativegrad,
2batchedgrad,2rotatedgrad,4signedmetaldisplacements. StageBconditional6forwards
onreal4698-atomGGR1GLG:2gradcenters,4signedmetaldisplacements. Frozen numerical
criteria and exactsources inplan. No gradientimplementation/callsdone yet.
Outputsarederivativesofmaskeddescriptor,notvalidatedphysicalforces. Do not
addrelaxation/entropy merelybecauseagradientorcurvaturematrixcanbecomputed.

GGRcommit1038b55;nextownscopedcommitcoversmultisitecode/result,nextplan.
Vaultandagentguideupdated. LatestemailacceptedrelaynotifiedGGRfailure;do not
spamrepeats. NativeH2001200809stillpendingunchanged;inspectlivejobs. Allnew
A5000jobsabovearecomplete.Noautomaticpromotion/push/per-analysispermissiongate.
Preserveotheragents'PLMjobs,watchers,dirtyfilesandlocks.


## Gradient implementation history

The original in-flight notes are superseded by the latest checkpoint above.
All recovery versions and failure receipts remain under workspaces; see the
completed gradient report for the full sequence. The next scientific task is
the declared response screen, not another rerun of passed engineering checks.
