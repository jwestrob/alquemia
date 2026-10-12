# Native GFN2 parameter lineage — 12 October 2026

## Finding

The completed corrected rescore invalidates the original accommodated canonical25 calibration. Its 2,796/2,796 fresh standalone GFN2 points are valid, but the reference gap changes from +7.333342 to −4.556602 model kcal/mol. Corrected numeric scores remain available; corrected accommodated class calls remain unavailable. The historical 207 correct / 0 wrong / 1 inconclusive / 17 unavailable result cannot serve as corrected accuracy evidence.

The released static PQQ mode and opt-in three-source envelope also consume the affected native La solvent term. Their reference scopes differ from strict208: the strict208 overlap cannot be transplanted into either reference or repaired with a universal +22.60 offset. The known +22.60006 kcal/mol shift belongs to the motivating 222-atom control only.

This bounded audit read the dispatch after Jacob confirmed vault synchronization. It wrote only this report and [LEDGER.json](LEDGER.json). It launched no molecular jobs, Slurm jobs, or tests; read no reserved labels; changed no production, shared code, historical artifacts, or vault documents. Root owns runtime enforcement and integration. No old LanM continuation was launched.

## Evidence checked directly

Eighteen selected snapshot artifacts match their entries in `/home/jwestrob/jwestrob/obsidian-vault/projects/lanthanide-binding/orca_parameter_audit_2026-10-12/MANIFEST.json`. These include the full La report, causal audit, strict208 calibration/validation/summary and component tables, and the separate transition-metal audit. This verifies the selected copied bytes, not the entire original compute archive.

Four actual historical La/vacuum endpoints were inspected: released static q9z4j7, strict canonical origin, strict alternate-fold transfer, and the PLM envelope candidate. For each, the ledger pins the manifest, geometry, template, runtime input, stdout, execution receipt and parameter export. Every sampled export has SHA256 `f3dac846a59238947d35c23ba5bdb40c19a34482fe9860e8e11bf9c00d719107`, exactly the original uncorrected parameter identity. The recorded input contains `ReadXTBParam false` and `WriteXTBParam true`; La retains shells `[5d, 6s, 6p]` with the affected original assignments. This establishes actual consumption along these paths rather than inferring it only from source names.

All four receipts record ORCA 6.1.1, driver SHA256 `38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83`, at `/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca`. Commands invoke that driver with `endpoint.runtime.inp` and `--bind-to none`. Geometry hashes, formal charge/multiplicity, solver input, command and receipt completion are recorded per sample. The receipts do not pin all solver helpers/shared libraries. Current runtime equivalence therefore remains a separate check; driver equality alone cannot establish it.

The tested corrected full JSON is pinned by [the October 9 record](../../native_gfn2_la_shellfix_20261009/PARAMETER_RECORD.json), SHA256 `a4f83ca72d58d6e80f634e61bb8e99b7af0f26c38cc4a9295e85215ecf72d05c`. Its only semantic changes are La `kcn`, `shpoly`, and `lgam`. Both source and staged copies match their recorded digest. This La-only artifact does not correct the separately diagnosed transition metals.

## Affected dependency paths

| Scope | Concrete dependencies | Consequence |
|---|---|---|
| Released static PQQ | `scripts/affordable_workflow.py standard` → `scripts/pqq_standard.py` → `params/pqq_fast_v1.json` → `compact_solvation_scanner.py` / `compact_solvation.py` → `run_orca_task_manifest.py` | Composite solvent term and native-derived static bands require parameter-consistent qualification. The static release is separate from the accommodated calibration. |
| Original strict accommodation | `workspaces/strict_native_pool_20260923/run_v1/REFERENCES.json` → `workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json` | Exact completed corrected canonical25 overlap makes the original accommodated classification unavailable. Reuse the completed recompute rather than resubmit its points. |
| Opt-in PLM envelope | `params/pqq_plm_envelope_v1.json` → `pqq_three_source_envelope.py` / `pqq_three_source_envelope_execution.py` → `motion_envelope_scalar.py` → shared native runner | Candidate selection, composite medians, ranges and reference decisions depend on the native La term. Its distinct `workspaces/motion_envelope_20260923/REFERENCE_PINNED_v1.json` must retain its own identity. |
| External PQQ production lineage | Motivating native-control receipts link to `EastRiver_PLM/revision_analysis/2026-09-23_Nikasha_PQQ/preparation_revision_4A/mapping_repair/...` | Direct affected source is established for the motivating control. This audit did not exhaustively inventory the concurrent external campaign or change its outputs. |

The static release explicitly marks `fast_PQQ_OMOL_GFN2_ALPB_v1` promoted for compatible explicit PQQ requests. Its calibration pin points to `workspaces/compact_solvation_20260920/full_v1/comparison_v1.json`, not the strict canonical25 reference. All current files named by its implementation pins match those pins in this inspection. That is historical implementation consistency; it does not qualify the scientifically incorrect La assignment or authorize continued classification under corrected parameters.

The PLM operations package separately pins envelope calibration, strict numerical qualification, rank qualification and MaxIter qualification. Matching solver/implementation pins cannot substitute for matching parameter identity. Original numerical qualifications and original observations remain historical evidence.

## Downstream claims requiring corrected qualification

- `diagnostics/strict_native_transfer_20260923/REPORT.md` and `diagnostics/overnight_discrimination_20260923/DELIVERY.md` report the historical 207/0/1/17 outcome and 199/1/8 → 207/0/1 accommodation comparison. The completed corrected audit supersedes using those as current accommodated accuracy claims.
- `diagnostics/nikasha_manuscript_20260922/PLM_ACCOMMODATION_DRAFT_20260923.md` carries those outcomes, claims an isolated useful accommodation contribution, and reports native-composite PLM medians/range reductions. Retain raw MACE structural observations, but energy-selected/composite conclusions need corrected provenance and calibration.
- `diagnostics/pqq_fast_release_20260920/REPORT.md` recommends the static native-GFN2 route and reports its 25/25 + 3/3 native-calibrated fidelity. This is an affected historical release statement; strict208's corrected overlap is not a measured replacement result for this distinct static release.
- `docs/PLM_PQQ_SOP.md`, `diagnostics/motion_envelope_transfer_20260923/REPORT.md` and `diagnostics/plm_candidate_overlay_20260923/REPORT.md` consume the envelope route or its scores. Preserve original rows and label new corrected results explicitly. No broad new biological labels were established by the unlabeled PLM cases.

No manuscript, SOP, result table or production parameter package was edited by this subtask. These paths are integration targets for root, with their inspected digests in the ledger.

## Components that remain distinct

**Standalone xTB.** The historical September standalone lane has its own manifest at `workspaces/standalone_xtb_reference_20260923/run_v1/manifest.json`: xTB 6.7.1 driver pin `481a1e7d8151fc08e0341b8e06ab60ac09cf1c2628738eeb83346271274bf8f4`; explicit parameter pin `f6f05c6c5264cb2f75c54f8ac6838044541c1e493b16f1a4d4ee4fa64de5f545`; package metadata, temperature, accuracy, charge and restart conventions retained. The parameter copies and package record were hash-checked here; the executable pin was read from provenance, not requalified. This branch is not an ORCA-native assignment consumer. Tests that transferred its scores onto native bands remain incompatible-reference diagnostics. The October strict208 rescore uses a separately recorded arm64 xTB build; it must not be relabeled as this cluster executable.

**MACE.** The stored MACE-OMOL checkpoint pin is `9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`, float64. The correction changes no MACE weights, raw MACE energies or proposal forces. It can change which geometry is selected by the composite. At q9z4j7 the stored accommodation MACE shift is +8.517684 kcal/mol; the corrected solvent shift is +3.568463 instead of −5.005229. Thus the old error partly compensated context/accommodation sensitivity. A successful parameter intervention does not resolve that remaining protocol problem.

**Explicit DFT.** The native-GFN2 JSON is not consumed by the explicit r2SCAN-3c/CPCM reference or prepared PBE0-D4 frozen-f FNR continuation. The latter's manifest and both exact inputs were inspected: PBE0-D4, RIJCOSX, DefGrid3, VeryTightSCF, EnGrad/TRAH/MORead; La ECP core 46 and Dy ECP core 55. Sharing an ORCA driver does not make these calculations GFN2. No DFT claim is newly validated by this separation, and numerical/state limitations remain independent.

**Geometry and biological metadata.** The synced completed summary retains the frozen `original_24start` comparator on original raw poses at 208/208. It was not rerun here, and differs from geometry evaluated on energy-selected poses. The 208 alternate folds inherit labels from 25 parent proteins. Original sequence/expression metadata and archived DFT columns remain separate from a corrected native-GFN2 score.

## Completed results versus remaining controls

The strict208 recomputation is complete: 25 canonical + 208 alternate pools, each three geometries × two metals × two media; 2,796 valid points, no fresh failures. All 225 declared folds remain represented, including 17 original unavailable preparations. Three selection profiles have complete numeric alternate scores and null corrected bands. Fourteen alternate selections and two canonical selections change under reselection. No additional calculation is required to establish that failure of calibration.

The original-context corrected origin-only ablation has its own +8.867782 gap and 203 correct / 3 wrong / 2 inconclusive result. The expanded accommodated panel has its own +3.472171 gap. They do not repair the original accommodated protocol. Context composition and reaccommodated coordinates change together; exact atom/cap/protonation and candidate mapping must precede the planned context × accommodation factorial.

Keep numerical exclusions explicit. The synced strict208 validation has 300 historical canonical SCF records (56 maximum-density overruns) and 2,496 alternate native SCF/population records unavailable locally. Eight direct corrected native boundary endpoints closely match energies, but four exceed density criteria; strictly accepted paired boundary scores remain null. The full original archive and all 2,796 output files were not rehashed by this audit.

Cu/Mn/Fe/Co/Ni four-array corrections, including `refocc`, are established for the pinned macOS build and finite states in the separate audit. All 14 tight corrected aquo parity controls pass; unchanged imports reproduce default energies/populations, while four original-setting cases remain unaccepted. Cluster reproduction remains required. Reference populations are not actual SCF occupations; requested multiplicity is not proof of physical spin energetics in spin-independent GFN2.

## Smallest safe enforcement integration

Use one parameter-identity helper, called by existing preparation/dry-run, launch and cache admission. Root owns implementation.

1. `run_orca_task_manifest.py`: preserve policy metadata during task normalization; check approved method/build, full task-local parameter digest and exact import directives before launch; recheck after execution and before cached acceptance. Record actual import evidence and executable/runtime identity. Gate by native GFN2 so DFT stays outside this contract.
2. `compact_solvation.py`: add a separately versioned import profile and full parameter digest to scientific keys. Imported-file diagnostics must replace the old assumption that `endpoint.runtime.xtb.json` is always emitted. Preserve historical recipe validation in immutable snapshots.
3. `motion_envelope_scalar.py` and analogous new native preparation: propagate parameter identity into cache keys, both solvent legs, candidate pools and reference/calibration identities. No old-band fallback or mixed original/corrected endpoint reuse.
4. `render_orca_runtime_input.py`: keep scientific templates unchanged; validate/record the contract while inserting allocation controls. Do not silently rewrite the Hamiltonian during rendering.
5. `affordable_workflow.py::dry_run/execute`: run the same preflight before submission/execution and repeat at launch. An explicit historical-control profile can retain original parameters for declared causal controls; missing identity should not silently choose installed defaults.

Full pinned JSON import, unchanged shell declarations, an exact permitted semantic diff, and scope-specific supported elements/builds are essential. New-metal work follows the dispatch's primary standalone xTB 6.7.1 + tblite 0.7.0 per-input-check lane until exact cluster-native reproduction qualifies a corresponding policy. No blanket claim that all orbital parameters or builds are correct follows from the La correction.
