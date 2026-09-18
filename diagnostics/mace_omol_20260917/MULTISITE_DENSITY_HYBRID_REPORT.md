# Frozen MACE hybrid fails the parvalbumin extension

**All five sites completed; all39 numerical checks pass; all six supporting parvalbumin ordering comparisons fail.** Baseline/default remains unchanged. This frozen candidate has not earned calibration or promotion.

| Supporting comparison | Difference kcal/mol | >0.02 |
|---|---:|---|
| PARV_4CPV_CD minus GGR_extended | -30.860274666 | fail |
| PARV_4CPV_CD minus GGR_2FW0 | -55.751496119 | fail |
| PARV_4CPV_CD minus GGR_2FVY | -32.059759458 | fail |
| PARV_4CPV_EF minus GGR_extended | -3.371289637 | fail |
| PARV_4CPV_EF minus GGR_2FW0 | -28.262511091 | fail |
| PARV_4CPV_EF minus GGR_2FVY | -4.570774429 | fail |

Parvalbumin contributes one biological group, two sites; the three GGR crystals remain one group. Its La/Ca direction is supporting cross-study evidence, not a direct same-assay gold standard. The failed alpha/2FW0 test remains a separate failure. No structure, threshold, donor, charge, water or site was selected to improve the result.

## Complete ordered vectors

| Site | Raw R, kcal/mol |
|---|---:|
| AEQ_1SL8_EF1 | -405434.172796625 |
| AEQ_1SL8_EF3 | -405422.617731380 |
| AEQ_1SL8_EF4 | -405431.161142844 |
| PARV_4CPV_CD | -405438.895340047 |
| PARV_4CPV_EF | -405411.406355019 |

The large common metal-energy offset is retained. Zero is not a threshold, and the baseline bands do not apply. Aequorin remains [EF1,EF3,EF4]; its protein-level weak Ca direction does not assign individual site labels. One-at-a-time substitution with the other sites fixed Ca is not a cooperative titration. These are two additional consumed development groups, with no blind validation.

## What the components show

Relative to GGR1GLGextended, the intrinsic quantum contrast shifts +52.95 kcal/mol for parvalbumin CD and +81.57 for EF. Environmental terms reverse that ordering; the MACE short-context differences are under 1 kcal/mol. Both parvalbumin sites therefore fail, though their margins differ. This locates the dominant terms in this model; it does not establish a unique cause or validate the supporting label across assay conditions. Exact component differences are retained in RESULT.json.

## Executed work and physical state

Ten native r2SCAN-3c responsive endpoints,20 MACE short full/core forwards,10 CHELPG fits,10 saved-density field/ESP queries,20 source/reference boundary initializations,63 static energies,60 field queries and75 response solves completed. No DFT retry or substituted output. The five actual AMOEBA/native parameter receipts were imported with exact physical-state and backend checks. Preparation preserves the actual ACE0–ALA1 bond, all background Ca ions, water unions, paired coordinates and endpoint charge ledgers.

All39 numerical/identity/rigid/refinement/radius checks pass. Largest response-refinement effect is3.411e-7 kcal/mol; rigid component effect6.367e-12; supporting-comparison radius effect0.001185. Aequorin raw R changes under the prescribed radius variants are explicitly retained (largest0.005017 kcal/mol); they are not binary affinity tests. No new partition geometry was tested. The earlier GGR partition check remains separately recorded.

The generating density responds to ff19SB plus the explicitly pinned amber19/tip3p Ca/water templates. Every retained water charge was verified equal to the original TIP3P template. That generating field is distinct from AMOEBA/GK evaluation, so the trial density is not self-consistent with the current hybrid functional. Projected CHELPG is still a GK source proxy. No combined gradient, mechanical correction, entropy, aquo reference or absolute decision is available.

## Measured cost

| Phase | Wall s | Allocated core-s | Actual CPU s | GPU s |
|---|---:|---:|---|---:|
| quantum, 1201164 | 426 | 27264 | 05:47:55; seconds precision | 0 |
| short, 1201165 | 208 | 3328 | 238.153 | 208 |
| density, 1201167 | 326 | 2608 | 2072.832 | 0 |
| native, 1201168 | 502 | 32128 | 2737.174 | 0 |

Total recorded allocation: **65,328 core-seconds and208 GPU-seconds**; summed job wall1462seconds. Quantum CPU was recorded at whole-second precision. This includes development controls on five sites, not a measured ordinary-production latency or matched baseline overhead. Local preparation, reused parameterization and earlier reference jobs are additional. Their receipts remain in the linked manifests; unmetered debugging is not claimed as zero.

Native kernels sum79.107079 wall seconds; native child processes173.166642 wall/2427.394123 CPU seconds. Native batch peakRSS1,854,748KiB. The rest of its502-second allocation includes launch, validation, parsing and collection. No H200 allocation or training was used.

## Implementation and validation

Scientific model remains `saved_responsive_trial_density_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1`. New comparison runner: `declared_source_graph_responsive_density_GK_POLAR_panel_v2`. The opt-in multisite input protocol is `source_graph_multisite_normalized_responsive_density_inputs_v2`; framework import is `source_multisite_AMOEBA2018_GK_verified_framework_import_v1`. Legacy defaults and archived execution snapshots remain intact.

New input4, runner3, grouped-comparison3 and native-panel3 tests pass (13distinct); actual native integration now runs with no final skip. Legacy input5, expansion6 and native-panel2 regressions pass. Two local preflights failed before scientific execution (missing Ca template and one-Ca validator); both partial preparations remain visible. An initial test compared paths of equal boundary copies; exact-content comparison fixed it without changing inputs.

**Judgments:** numerical consistency passes the prescribed tests; predictive transfer fails; execution is practical for this pilot, while matched production overhead remains unmeasured. Retain baseline and investigate the environmental representation before another accuracy trial.

## Runnable replay

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/multisite_density_native_v1/implementation/mace_density_panel.py dry-run --manifest workspaces/mace_omol_20260917/multisite_density_native_v1/manifest.json

/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_multisite_density_report.py --collection workspaces/mace_omol_20260917/multisite_density_native_v1/collection_job_1201168.json --output workspaces/mace_omol_20260917/multisite_density_comparison_review_v1
```

Use a fresh output directory for report replay. The completed native manifest must not be resubmitted. Input/config/manifest/collection hashes, states, full numerical checks and costs are in [RESULT](MULTISITE_DENSITY_HYBRID_RESULT.json). Exact original stage configurations are under `workspaces/mace_omol_20260917/multisite_density_*`; existing prepare/dry-run/execute/collect commands accept explicit paths.
