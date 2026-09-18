# Complete conductor transfer: discrimination comparison

Protocol: `trial_density_MACE_AMOEBA_frozen_response_ddCPCM_v1`. Production baseline unchanged. This compares the prior
GK hybrid with its frozen-response conductor correction, not two calibrated classifiers.
All structures are consumed development cases; crystal replicates and sites are grouped.

| Comparison | Prior GK | Conductor primary | Conductor refined | Refinement change |
|---|---:|---:|---:|---:|
| ALPHA_1F6S − GGR_extended | +10.127062 | +3.769137 | +3.754399 | -0.014738 |
| ALPHA_1F6S − GGR_2FW0 | -14.764159 | -7.050328 | -7.061140 | -0.010812 |
| ALPHA_1F6S − GGR_2FVY | +8.927578 | -2.325348 | -2.346042 | -0.020695 |
| ALPHA_6IP9 − GGR_extended | +9.252283 | -2.225645 | -2.246291 | -0.020646 |
| ALPHA_6IP9 − GGR_2FW0 | -15.638938 | -13.045111 | -13.061830 | -0.016720 |
| ALPHA_6IP9 − GGR_2FVY | +8.052798 | -8.320130 | -8.346733 | -0.026603 |
| PARV_4CPV_CD − GGR_extended | -30.860275 | -8.838376 | -8.828939 | +0.009437 |
| PARV_4CPV_CD − GGR_2FW0 | -55.751496 | -19.657841 | -19.644478 | +0.013363 |
| PARV_4CPV_CD − GGR_2FVY | -32.059759 | -14.932860 | -14.929380 | +0.003480 |
| PARV_4CPV_EF − GGR_extended | -3.371290 | +16.621807 | +16.599991 | -0.021816 |
| PARV_4CPV_EF − GGR_2FW0 | -28.262511 | +5.802342 | +5.784452 | -0.017890 |
| PARV_4CPV_EF − GGR_2FVY | -4.570774 | +10.527322 | +10.499549 | -0.027773 |

Values are kcal/mol differences in R = E_Ca − E_La. The predeclared directional
screen is >0.02. These twelve comparisons are not twelve independent biological tests.
GGR has direct same-assay Ca-favoring evidence; alpha is a qualified strong-site
comparison and parvalbumin is supporting cross-study evidence. No absolute bands.

**Outcome: no overall discrimination gain.** The prior and new methods each pass
4/12 directional comparisons. Gains on parvalbumin EF are offset by losses on
alpha-lactalbumin. Retain the baseline and close this conductor challenger; no
additional panel expansion is planned. The small paired refinement changes
(maximum margin change0.027773 kcal/mol) do not explain the large failed margins.

## Directional summary

- alpha: prior 4/6; primary 1/6; refined 1/6 (6 available).
- parvalbumin: prior 0/6; primary 3/6; refined 3/6 (6 available).

## Numerical qualification and cost

Conductor endpoint/rotation qualification previously failed; not waived.
New native-functional and paired-refinement checks: 105/105 passing.
Source representation, solved-state and contraction checks remain in each endpoint receipt.
A small paired grid change does not erase earlier endpoint or rigid-transform failures.
No new DFT or MACE inference. Native GK evaluations reuse four exact prior endpoints.
Scheduler allocation cost is recorded separately from worker timings; cache preparation
and prior development costs are not implied to be free.

Executed56/56 new conductor solves, zero failures. Seven jobs used3436 summed
allocation-wall seconds,219904 allocated core-seconds and209515 reported CPU-seconds;
zero GPU-seconds. Longest job810s, largest worker RSS about2.6GiB. These are
development costs, not a matched production endpoint-pair timing. Four real-fixture
parser/algebra/cache/output tests pass (1.681s); all105 new functional and paired
refinement checks pass. Earlier endpoint/rotation failures remain.

## Reproduction

From the repository root (choose a fresh output filename):

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_conductor_discrimination.py collect --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/conductor_discrimination_v1/manifest.json --output workspaces/mace_omol_20260917/conductor_discrimination_v1/recollected.json
```

Collection: `/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/conductor_discrimination_v1/collection_final.json`.
Manifest SHA256: `6577e76ea5f43d3742471080fa230d46073fa039e982d6d8b137cf14c565ab36`. Unrounded values and components are retained.
