# Alquemia: current operating guide for agents

**Updated 2026-09-16.** Start here for new work. Jacob has approved contained
pilots autonomously; see root AGENTS.md for the instruction superseding older
per-pilot approval language. This guide supersedes older
operational/status prose; dated experiments and their numerical records remain
immutable. The baseline remains the default. Broad La/Ca affinity discrimination
has not been established. The [active MACE goal](../diagnostics/mace_discriminator_goal_20260916/GOAL.md) continues through intermediate failed pilots; production remains unchanged.

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
experiments. Record the agreed scientific scope in project notes and follow
the user's AGENTS.md analysis policy: routine implementation/recovery proceeds
within that scope, while materially new analyses need agreement. The user has
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
