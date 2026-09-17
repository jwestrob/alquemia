# Active MACE discriminator checkpoint — 2026-09-17

**Goal active.** Jacob authorized autonomous analyses and available resources.
The per-analysis permission gate has been removed from home/project AGENTS.
Preserve baseline/default, immutable results and other agents' jobs. No automatic
promotion or push. Recover the active goal and this checkpoint after compaction.

## Current candidate: charge-feature ablation descriptor

Protocol `mace_omol_intact_charge_feature_ablation_descriptor_v1`.
Same pinned OMOL100M weights, whole-chain coordinates, protonation, spin,
physical charges and water inventory; raw total-charge embedding is zeroed
before the native joint projection. This is explicitly a **learned descriptor,
not a quantum energy for those physical charges**. No baseline band or affinity
zero applies. See [plan](../mace_omol_20260917/CHARGE_ABLATION_PLAN.md).

Job **1200828 completed**, 42/42 forwards, no failures. All numerical and
three-direction development gates pass. XoxF−MxaF +44.329179 model kcal;
alpha1F6S−GGR +12.741688; alpha6IP9−GGR +4.112398. GGR sodium perturbation
changes its score by -1.07384e-8 (tolerance0.1). The latter invariance is imposed
by the architecture; it does not demonstrate physical electrostatics.
All cases consumed; two alpha structures are one qualified affinity comparison,
and PQQ functional class is a separate stratum. Broad validation remains open.
[Report](../mace_omol_20260917/CHARGE_ABLATION_REPORT.md).

Actual job cost798 GPU/wall seconds,12768 allocated core-seconds,894.358
reported CPU-seconds;409.63256 summed evaluation seconds;11536478720bytes
peak GPU allocation,2876296KiB peak processRSS. Local receipts separate.
Eight real-fixture/actual-forward guard tests pass57.308s.

## Live execution: conditional canonical test

**Job1200830 running** after all development gates passed; inspect live state.
Recent check74/100completed, no failures; do not infer final calibration yet.
100 new forwards for all25canonical calibration proteins,8exact modified
crystal reuses, all28evidence rows retained. 1KB0 remains unsupported.
OneA5000/16CPUs/64474MiB, existing runner, no project compute/time budget.

Workspace: `workspaces/mace_omol_20260917/charge_ablation_canonical_v1/`
ManifestSHA:2415b4dacd755530dec888017d76afe8895708213a36ef09497c7264e341476f.
Preparation and frozen dry-run pass. Exact submission argv in submission.json.
The frozen wrapper will collect into collection_job_1200830.json.
Do not duplicate the job, modify its manifest or delete executor locks.

Use [commands](../mace_omol_20260917/CHARGE_ABLATION_COMMANDS.md).
After collection, use its frozen implementation's mace_omol_ablation_panel.py
report command with --collection and a new --output directory. The report
retains all denominators and model units. All25calibration cases must be valid
and min(La)−max(Ca)>0.02 before new candidate-specific bands. Unsupported1KB0
cannot count as a correct transfer. Do not retune on inconvenient outputs.

Source development workspace: charge_ablation_development_v2/;
manifestSHA78c781671295658656111f94f9b7cd6625b5bb696ab85e726fe70ea5de2d30c5.
Full report: charge_ablation_report_v1/result.json (SHA
0d7b7f325bf8032b1f11c1811c994684f028d76e0810fa335069bd830659a947).
The unexecuted development_v1 was superseded to fix receipt row counting:
MACE expands graph charge to atom rows before embedding. No scientific change.
Component_v3 checks201charge categories and exact zero-charge/original-spin
projection; no learned parameter changed. V1import collision/V2wrong reference
feature order remain recorded technical failures, without molecular calls.

## Earlier native whole-chain candidate: rejected for promotion

`mace_omol_intact_chain_matched_coordination_v1` passed the original five-case
three-direction development panel, then failed the larger calibration and
charged-spectator consistency test. Original results remain immutable.

- Job1200819:104newforwards+8reuses, all25calibration scores, separation gap
  -306.006296kcal/mol. No bands. Two crystal raw scores, unsupported1KB0.
  [Report](../mace_omol_20260917/INTACT_CANONICAL_REPORT.md).
- Job1200823:20forwards, four scores shift4–15kcal/mol under a disconnected
  Na+; XoxF/MxaF reverses. [Report](../mace_omol_20260917/INTACT_SPECTATOR_REPORT.md).
- Actual100Mcheckpoint has identical charge features for -100..-6 and another
  group +11..+100. MxaF's invariance reflects that aliasing. This is observed
  checkpoint behavior, not an inferred training history. [Audit](../mace_omol_20260917/CHARGE_EMBEDDING_AUDIT_REPORT.md).
- Saved-output locality audit:18closure checks pass; beyond18Å changes are
  <=3.66e-11kcal/mol. Native matched subtraction cancels distant additive terms,
  but not nonlinear global charge conditioning. [Audit](../mace_omol_20260917/INTACT_LOCALITY_REPORT.md).
- A common separately evaluated ionic reference cannot fix the raw bound
  XoxF−MxaF contrast (-623.299274kcal/mol). No new reference energy was computed.
  [Screen](../mace_omol_20260917/SEPARATED_REFERENCE_SCREEN_REPORT.md).

## Preparation and resource safeguards

1KB0 whole-chain input has false peptide connections at4.750089Å and18.969263Å.
Its original four native outputs are diagnostics only. Mandatory native report
is intact_panel_reporting_source_v2, with its geometry audit. Original frozen
reporterV1is superseded by REPORTING_SUPERSEDED.json; do not use it to classify1KB0.
Strict new policy `omol_canonical_chain_A_ff19sb_H_peptide_connectivity_v2` gives
27supported/1unsupported; all1718bonds in the original five-case preparations
pass independently. Old fixed-core baseline1KB0 remains unchanged.

Exact1024edge/product batching is qualified against native CPU energies.
Original native H200 job1200809 remains pending on its unchanged14-task manifest;
no cancellation/duplicate. Completed engineering summaryV7 totals248successful
forwards,2OOMs,7736GPU-s,187008allocated core-s,25495.954reported actual CPU-s.
This includes4invalid-preparation diagnostic forwards; local receipts separate.
New canonical job cost is not yet included. No newDFT/solver/training/forces.

## Recovery and remaining work

Read the actual canonical result next. Numerical credibility, usefulness and
cost are separate judgments. Canonical composition alone already separates
classes, and independent non-PQQ labels remain scarce; consult the existing
challenge curation without inventing new site labels. Earlier POLAR/hybrid/
mechanical failures remain recorded; a failed pilot does not complete the goal.

Email accepted by local mailer20260917T134908Z_charge_ablation_pass; vault note
updated with actual results and current job. Prior own commits19cede4,ed2e65c,
2fcf2f3. Conditional canonical support/report changes are scoped separately.
Many other repository edits are concurrent/historical; never blanket-stage.


## Prepared-input interface now works

New scripts/mace_omol_prepared.py audits source-backed whole-chain preparations,
creates four descriptor tasks or explicit actual reuses, and reports raw values.
Exact protein/cofactor/water replay and peptide connectivity required. Single
metal-bearing chainA only; recorded source omissions remain visible. No new
protonation/model/inference was used to validate this interface.

IMPORTANT driver is the existing OpenMM environment:
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python.
Its executable/versions are pinned. The existing executor launches actual MACE
workers in the separate software-recorded MACE venv. Do not install OpenMM into
that working venv; a failed first audit/import is preserved and recovered with
the correct existing driver.

Under workspaces/mace_omol_20260917/: prepared_interface_pqq_v1 has zero new
tasks and four actual reuses; prepared_interface_pqq_report_v1 returns the exact
17.579955566129197 descriptor with classification unavailable. Four real tests
pass28.325s, including actual1KB0 rejection and corrupted-coordinate rejection.
prepared_interface_ggr_fresh_v1 contains four unexecuted tasks; frozen dry-run
passes. DO NOT submit duplicate GGR calculations; this is only an interface
example. See [commands](../mace_omol_20260917/PREPARED_INPUT_COMMANDS.md).

Current interface deliberately has no classification/reference backend until
the running canonical test supplies a compatible passing calibration. Finishing
that conditional integration and evaluating broader evidence remain ahead.
Latest completed own commit287d3d0; this interface work is scoped separately.


## Latest factorization work (no new inference)

[Proof](../mace_omol_20260917/DISCONNECTED_FACTORIZATION_REPORT.md) passes60checks
with maxerror2.36e-8model kcal. The exact masked descriptor can use two bound
forwards and fixed disconnected-metal node+embedding terms. Optional
--factorization is implemented in the prepared-input interface; four-call path
retained. Authoritative factorization_report_v2/result.json SHA
35ba8bb53b209f650b8392ff03f0c426d0f95455b39b18c082e25aabb5434fbe.
V1had only cross-NumPy auxiliary sum discrepancies; v2uses math.fsum. Constants,
scores and criteria unchanged. Original records/failure timing retained.

prepared_interface_pqq_two_call_v2 reuses two actual bound receipts; its report
is complete. Three factorization tests pass44.009s. No new model call. Use the
OpenMM driver for prepared-input operations; MACE workers remain separately
pinned. DO NOT infer detached components that the two-call path does not compute.

Canonical1200830 still runs the original100tasks. Before its final calibration
was read, FACTORIZATION_NUMERICAL_BANDS_ADDENDUM.md declared separate numeric
band records for two/four-call evaluation, using the same25calibration cases,
only after full-panel equivalence. No changed labels, coefficient or epsilon.
Next implement/execute that full-panel proof and compatible reference integration
once its actual collection/report is complete. No calibration integration exists
yet. Latest own committed interface work30cb9b4; factorization changes separate.
