# Shared learned-neutral feature: numerically correct, partition test fails

The eight core calculations completed (job 1200901). All 16 numerical checks
pass, but the GGR hybrid partition shift is **+10.912283925 kcal-scale**, above
the declared 2 limit. Do not launch this candidate's six conditional whole
calculations or sweep charge categories. Baseline and earlier results remain
unchanged; the wider MACE discriminator goal remains active.

Protocol: `omol_shared_neutral_charge_feature_descriptor_v1`. Both endpoints
use the checkpoint's actual learned charge-zero vector (table index 100), while
physical charges, spin, coordinates, parameters and source/cap mappings remain
unchanged. This is a distinct learned descriptor, without an affinity reference.

| Connected minus extended GGR contribution | kcal-scale |
|---|---:|
| DFT/CPCM core contrast | -7.343500873 |
| Shared-neutral learned core contrast | -18.255784798 |
| DFT core minus learned core | **+10.912283925** |

The identical full-system term cancels in this partition difference. No new
whole-system inference is required to reject the declared expression. The
CPCM-core/vacuum-learned mismatch remains a possible contributor, not a proven
unique cause. No solvation energy, mechanical correction or classification is
inferred from this result.

All four physical-charge-zero La cores reproduce exact native archived energies;
maximum gradient component error is 9.77e-15 eV/Angstrom. Eight native-readout
closures also pass. Projected center gradients versus DFT are retained in the
result, but do not validate curvature or relaxation.

## Execution and checks

V1 job 1200893 failed preflight before any model call: replaying an existing
native qualification on the compute host changed an auxiliary floating error
by 8.88e-16. V2 retains exact raw receipts, arrays, pins, fixed thresholds and
gate decisions while allowing host-derived floating summaries to differ. All
eight scientific tasks are unchanged; implementation/cache hashes are new.
The independent neutral adapter also has explicit dispatch in the existing
runner. Previous energy and gradient adapters remain unchanged.

Three new real-artifact tests and six existing native regressions pass. The
actual-result test was explicitly skipped before execution and then passed
against the real result (0.139 seconds). No synthetic scientific output was used.

Including preflight failure: **82 GPU allocation-seconds, 1,312 allocated
core-seconds, 93.244 reported CPU-seconds**. Successful model evaluations total
6.595379632 seconds, with 1,031,137,792 bytes peak GPU allocation. Local prepare,
report and test resources are recorded separately. No new DFT, solver, training
or optimization calls. Cumulative engineering V14 records 437 successful model
calls, four failed model calls, 12,424 GPU allocation-seconds and 262,016 allocated
core-seconds; preflight failures are also retained.

## Artifacts and replay

Under `workspaces/mace_omol_20260917/`: `shared_neutral_core_v2/manifest.json`
(SHA256 `4b3420ef72fc742b25cc407b1992073e3a55da2d28f25f15c9055781bb7db814`),
`shared_neutral_core_report_v1/result.json`, `shared_neutral_cost_v1.json`,
`shared_neutral_reference_recovery_v2.json` and `intact_engineering_status_v14.json`.
[Compact result](SHARED_NEUTRAL_FEATURE_RESULT.json) includes source pins,
unrounded components and all checks.

Replay into a fresh output directory from the repository root:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/shared_neutral_core_v2/implementation/mace_omol_neutral_run.py report --manifest workspaces/mace_omol_20260917/shared_neutral_core_v2/manifest.json --output workspaces/mace_omol_20260917/shared_neutral_core_report_replay_v1
```

Recommendation: retain the production baseline. This shared-neutral subtractive
candidate did not pass physical consistency; no predictive improvement is claimed.
