# Vacuum hybrid: partition passes; accuracy test fails

All six whole-protein MACE calls completed. The separate vacuum hybrid passes
11/11 numerical checks and its2kcal partition diagnostic, but only **1/4**
predeclared alpha-minus-GGR ordering comparisons passes. Retain the baseline;
this fixed candidate is not a validated improvement.

| Alpha structure | GGR representation | Alpha minus GGR, kcal-scale | Gate |
|---|---|---:|---|
| 1F6S | Extended | +0.558714448 | Pass |
| 1F6S | Connected | -1.271091793 | Fail |
| 6IP9 | Extended | -0.265848708 | Fail |
| 6IP9 | Connected | -2.095654950 | Fail |

The fixed criterion is each difference>0.02; the all-four gate fails. Both alpha
structures remain one condition-qualified biological group, and both GGR cores
one direct-affinity group. These are consumed development cases with different
assay conditions. No broad or blind accuracy denominator is implied.

## What this separates

The matched vacuum diagnostic established that solvent mismatch accounted for
18.3305 of the earlier20.1603kcal partition discrepancy. The explicit full-system
test confirms the remaining shift is **1.829806242**, because the same whole
GGR term cancels. Correcting that mismatch improves representation consistency;
it does not by itself deliver the desired biological ordering. Some ordering
margins are smaller than this remaining representation sensitivity.

Protocol `masked_omol_subtractive_context_DFT_vacuum_v1` uses
`H_M = E_DFT,vacuum(core,M) + T_mask(full,M) - T_mask(core,M)`.
Unrounded component energies and baseline/CPCM comparisons are retained in
the result. No fitted weight, changed label, favorable-core selection,
solvent term, physical zero or old threshold was introduced. Aqueous and
calibrated scores, relaxation and entropy remain unavailable.

All core/full coordinates share the original source H positions, charges,
protonation, waters and assembly. Known H-bond-length defects remain a limitation:
the earlier normalized-H whole-protein scores cannot substitute for these
outputs while reusing original-H DFT. Any matched normalization experiment must
be separately versioned and preserve this failed result.

## Execution and verification

Job1200924: six successful energy-only masked-OMOL forwards on oneA5000,
16CPUs/64474MiB host memory; no failures. All eight DFT and eight core-MACE
centers were reused exactly. No new DFT, solver, gradient, optimization or training.
Existing scalar worker/adapters are unchanged; a new comparison module and
small existing-runner dispatch handle this protocol.

**98 GPU allocation-seconds,1,568 allocated core-seconds,117.848 reported
CPU-seconds**. Model evaluations total47.057306305seconds; peak GPU allocation
6,186,040,832bytes. The preceding vacuum DFT diagnostic cost39,552core-seconds
and35,369CPU-seconds additionally. Local preparation/report/test receipts are
separate. No production speedup is claimed.

Three distinct new real-fixture tests and six existing native regressions pass.
The new tests cover exact original-H preparation,
rejection of corrupted actual charge/cache metadata, and actual executed result
and gate replay. Actual integration was skipped before inference, then passed
in3.895seconds. An initial dry-run CLI rejected an execute-only argument before
any calculation; the corrected frozen dry-run passed without source changes.

## Artifacts and replay

Under `workspaces/mace_omol_20260917/`: `vacuum_hybrid_v1/manifest.json`
SHA256 `8d054667b104384a49f3cbaa7b97ddd712b8c620c4660f2907f79afa4bc72d7e`,
actual execution receipts, `vacuum_hybrid_report_v1/result.json`,
`vacuum_hybrid_cost_v2.json` and `intact_engineering_status_v15.json`.
[Compact result](VACUUM_HYBRID_RESULT.json) includes all components and gates.

From the repository root, replay into a fresh output directory:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/vacuum_hybrid_v1/implementation/mace_omol_vacuum_hybrid.py report --manifest workspaces/mace_omol_20260917/vacuum_hybrid_v1/manifest.json --output workspaces/mace_omol_20260917/vacuum_hybrid_report_replay_v1
```

The MACE goal remains active. The solvent finding is useful and reproducible;
this unchanged vacuum hybrid has not earned an accuracy claim or promotion.
