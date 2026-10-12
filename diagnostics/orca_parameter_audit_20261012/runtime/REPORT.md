# Linux native-GFN2 parameter runtime qualification — 12 October 2026

All **60 declared native controls completed** on the pinned Linux ORCA 6.1.1 runtime. **52 passed every original state, process, population and printed SCF criterion**. Eight default/roundtrip endpoints retain their density failures. No retry, tolerance change, production edit or independent-solver call was made.

## Results

| Control | Actual result |
|---|---|
| Original Ca/La eight-cell replication | 8/8 accepted; energies and all 222 atomic charges reproduce archived Linux controls |
| Corrected full La JSON imported into Ca, both media | 2/2 accepted; energies and all 222 charges exactly equal current defaults |
| Cu1/Cu2/Mn2/Fe2/Fe3/Co2/Ni2 corrected aquo, both media | 14/14 accepted against both archived xTB 6.7.1 and tblite 0.7.0 |
| Unchanged aquo export/import roundtrip | 14/14 raw numerical matches; 10/14 pairs accepted under all density gates |
| Ca2/Mg2/Sr2/Zn2 default aquo | All 8 native endpoints admitted; 7/8 upstream parity passes; Zn water charge parity fails |

Corrected aquo maximum residuals are 5.27820134e-08 Eh and 7.27960941e-07 e, below the frozen 1e-6 Eh/1e-5 e thresholds. Unchanged roundtrip thresholds remain 1e-8 Eh/1e-7 e. Comparisons reuse accepted archived upstream outputs; no Linux xTB/tblite execution is implied.

Default and unchanged-roundtrip Cu2 vacuum/water, Fe2 water and Fe3 water fail printed density limits despite normal termination. Their exact roundtrip agreement does not admit them. Corrected_tight runs pass the original density criteria. Zinc water retains a maximum all-atom charge residual 8.74284557e-5 e; its small energy residual does not rescue parity.

The original La correction is exactly three arrays with declared shells and all other fields unchanged. Fresh corrected La reproduces the archived result: the motivating Ca−La solvent correction changes by **+22.600061589 kcal/mol**. This example-specific correction does not restore old calibration or establish protein prediction.

## Protocol, identity and evidence

The immutable original plan declared 58 calls. Before execution, PLAN_EXTENSION_Ca2.json declared two unchanged-Ca full-file import controls and the new ceiling 60. Stage1 used the exact existing eight fixtures. The aquo 50 geometry and native input byte hashes match the synced source artifacts. The source inputs establish TolE 1e-10 for default/roundtrip and 1e-12 for corrected_tight, overriding the task prose ambiguity before stage2. Other solver/state/temperature settings and all gates were retained.

RUNTIME_BINDING.json pins the actual driver, startup/SCF helpers, shared libraries and environment. OBSERVED_RUNTIME_FINAL.json adds actually observed orca_guess and dash, for 12 runtime files. No preregistration of those two additional observed helpers is claimed. All initial pins still match; the separate guard binds the final full identity. APPROVAL_INPUTS.json supplies each real corrected import, exact full-file SHA, medium and printed-state convention. La/Ca use full-file SHA a4f83ca72d58d6e80f634e61bb8e99b7af0f26c38cc4a9295e85215ecf72d05c.

RESULTS.json preserves raw versus accepted comparisons, all 60 receipt pins, original failures, source tables, finite manifests and costs. The actual full-precision property populations are compared, including charge sums and atom indices. Each imported file is byte checked before/after and stdout records the import. Ca2 receipts retain harmless copied null parameter_source/empty allowed_input_changes metadata; actual hashes, import lines and the separate declared extension establish the real change. No historical receipt was rewritten.

## Cost and limits

Actual summed endpoint wall time: 753.687 s. Measured reaped child CPU time for the 50 aquo plus 2 Ca controls: 55.951 s. CPU time for the first 8 was not separately captured; their measured wall times remain in receipts. Peak owned process-group RSS: 298131456 bytes. At most 4 concurrent single-thread numerical children; every observed process had 1 thread. Zero Slurm jobs/GPU/new independent-solver calls; no owned processes remain.

These controls qualify this runtime and explicit parameter/state combinations for further guarded research. Distinct target-metal full JSON profiles do not qualify mixed-metal unions. Spin-independent GFN2 does not establish physical spin-state energetics. Protein affinities, discrimination, accommodation/calibration and untested builds remain separate questions.

## Independent review and new guarded entry point

The independent review in `../runtime_review/` checked all 60 raw stdout/property pairs and 720 artifact hashes, including exact Ca import equivalence. It found no numerical discrepancy. The separately implemented `scripts/guarded_native_gfn2_dispatch.py` now requires shared parameter staging, preflight and postflight plus separate scientific admission. Final validation records seven passing real-input/rejection/cleanup checks in GUARDED_DISPATCH_VALIDATION_V2.json, with zero added molecular calls. Its execution branch remains unexecuted on chemistry; see GUARDED_DISPATCH.md. The final reviewed dispatcher SHA is `581af16b3184e883a2d46c3c0dea9396c28016e304086e41b6fab23f9c7ae7ff`.
