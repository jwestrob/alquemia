# Active MACE discriminator checkpoint — 2026-09-17

**Goal active; candidate milestone reached.** Jacob authorized autonomous
analyses and available resources. Home/project AGENTS now remove the old
per-analysis permission gate. Preserve baseline/default, immutable studies,
other agents' edits/jobs. No push/promotion. No project CPU/time/token budget.
Recover GOAL.md and inspect live jobs before resuming.

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
