# Multisite masked-MACE test: physical preparation works; supporting gate fails

Ten of ten real endpoint calls completed; all numerical checks pass. All-Ca energies are identical across the different site/atom orderings within each protein (reported difference exactly0). Parvalbumin passes4/6supporting contrasts against the threeGGR structures; the predeclared all-case gate is **false**. Aequorin has no site-resolved direction test.

| Protein/site, frozen order | Descriptor, model kcal |
|---|---:|
| PARV_4CPV_CD | 27.431183345 |
| PARV_4CPV_EF | 64.347646812 |
| AEQ_1SL8_EF1 | 13.047672919 |
| AEQ_1SL8_EF3 | 33.173348880 |
| AEQ_1SL8_EF4 | 51.564143761 |

Parvalbumin EF exceeds all threeGGR scores; CD exceeds only1GLG. These are supporting cross-study comparisons, not same-assay gold labels. The two sites remain one protein observation; no sub-kcal site-order claim is tested. Aequorin is reported as [EF1,EF3,EF4], without selecting a favorable site or interpreting the vector as reproduction of a global/cooperative binding assay. No PQQ band or universal affinity zero is used.

## What was implemented

New policy `omol_intact_multisite_fixed_background_Ca_source_acetyl_ff19sb_v1` retains all other calcium ions, a fixed union of archived site waters and the actual4CPVACE0acetyl group. Raw heavy-atom mapping and covalent connections replay exactly; ff19SB matching, charge/electron counts and paired coordinates pass. No heavy atom optimization or reconstruction was done. The10missingN-terminal1SL8residues remain explicitly absent.

4CPV has1611atoms, onebackgroundCa, onewater, physicalQCa−3/QLa−2. 1SL8has2866atoms, twobackgroundCa, threewaters, QCa−4/QLa−3. Only the named metal changes. The baseline and prior single-metal paths remain unchanged; background-calcium support requires explicit indices.

Initial physical-preparationV1stopped at the old single-metal state guard for all five cases, before inference. V2added explicit background inventory to that guard and prepared all cases. Both attempts are retained. Seven source/preparation/result tests passed14.973s; fourlegacyPQQ/GGRinterface tests passed28.313s; the actual multisite-result regression passed9.923s. None skipped. These parser/preparation tests are distinct from the ten executed GPU calls.

## Actual cost and scope

Jobs1200851–1200855 completed on A5000s: **216GPUallocation-seconds**,3456allocatedcore-seconds,221.486reportedactualCPU-seconds. Ten forwards took63.758972model-evaluation seconds total; summedworkerwall80.695251seconds. PeakallocatedGPUmemory3,936,891,392bytes; peakworkerRSS1,415,676KiB. Localpreparation/dryrun/report/testcosts are separately retained. CostV2corrects an auxiliary test-resource hash initially captured while that test was running; scores and measured allocation totals are unchanged.

Candidate products are under `workspaces/mace_omol_20260917/`: `multisite_prepared_v2/`, `multisite_scoring_v1/`, `multisite_reports_v1/`, and `multisite_comparison_v1/result.json`. Full manifests, endpoint energies, failed preparation attempts, source mappings and execution receipts remain pinned. No newDFT, solver, training, gradients or relaxation.

## Judgment and next operation

The preparation/scoring implementation is numerically credible and affordable on these inputs. Predictive generality remains unestablished: the expanded GGR and supporting parvalbumin gates fail. Retain the production baseline and the separate opt-in PQQ research classifier. Do not promote this as a broadly validated affinity discriminator.

The next declared development experiment is MASKED_GRADIENT_PLAN.md: obtain exact derivatives of the same masked descriptor with bounded real-coordinate checks. This may help design a structural-response score; it supplies no automatic mechanical correction. No gradient calculations have run yet.

An existing result can be audited without new inference:

```bash
MACE_DRIVER=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$MACE_DRIVER" scripts/mace_omol_prepared.py audit \
  --preparation "$PWD/workspaces/mace_omol_20260917/multisite_prepared_v2/PARV_4CPV_CD/preparation.json"
```

New preparations take explicit `mace_omol_multisite.py --recipe --site-key --agreement --output` arguments. Existing exact recipes are `multisite_recipes_v1/PARV_4CPV.json` and `AEQ_1SL8.json`; outputs must be new. Use the established prepared-input prepare/dry-run/execute/collect/report operations, with factorizationV2 and no PQQ calibration for these cases.
