# Full canonical check: union context × adaptive accommodation

The combined method retains **25/25 canonical classifications and3/3 consumed
crystal transfers**, with an independently frozen canonical-only gap of
**7.240239674 model kcal/mol**. Replaying the already consumed common8 pilot
after that freeze gives8/8 correct, including all four noncanonical hard-source
folds. This preserves the pilot benefit under the method’s own reference. It
does not establish broad biological validation or justify automatic promotion.

## Frozen reference

Exact protocol: `nikasha_union_adaptive_minimal_common_geometry_native_OMOL_GFN2_ALPB_v1`.
Reference: `Nikasha_union_adaptive_minimal_canonical25_v1`. Both mathematical
and operational variants retain all25, exclude crystals and noncanonical folds
from calibration, and use the unchanged0.02kcal minimum-gap/extrema rule. No
threshold padding, fitted coefficient or favorable source substitution.

| Variant | Ca supported maximum R | La supported minimum R | Gap |
|---|---:|---:|---:|
| mathematical | -405463.70905721450 | -405456.46881754074 | 7.240239674 |
| operational | -405463.70905721450 | -405456.46881754074 | 7.240239674 |

Previous canonical gaps were5.076366 for released static,2.546690 for original
adaptive, and11.066622 for static union. Thus the combination improves separation
relative to the original adaptive model but has a narrower training gap than
static union. Training gap alone is not evidence of improved prediction; static
union had added abstentions in the prior fold challenge. The four noncanonical
pilot cases are encouraging development checks, not fresh blind validation.

## All28 source outcomes

R = E_Ca−E_La; the table adds the same+405000 display offset to all raw values.
Machine artifacts retain the exact unrounded components and R. Larger values
are La-like; this offset is not an aquo reference or binding-free-energy zero.

| Source | Role / known class | Operational R +405000 | Own-reference outcome |
|---|---|---:|---|
| a0a3f2yly8-pqq-la_model | canonical / La | -444.233240 | correct |
| a0acd6b9f2-pqq-la_model | canonical / La | -444.092074 | correct |
| c5atj3-pqq-la_model | canonical / La | -420.953729 | correct |
| c5axv8-pqq-la_model | canonical / La | -456.468818 | correct |
| c5b120-pqq-la_model | canonical / La | -440.544556 | correct |
| i0jwn7-pqq-la_model | canonical / La | -429.476450 | correct |
| mmol_1770-pqq-la_model | canonical / La | -440.711102 | correct |
| mmol_2048-pqq-la_model | canonical / La | -438.417922 | correct |
| q88jh0-pqq-la_model | canonical / La | -449.587814 | correct |
| q89gy2-pqq-la_model | canonical / La | -441.912653 | correct |
| q92wy9-pqq-la_model | canonical / La | -431.762979 | correct |
| a8r3s4-pqq-la_model | canonical / Ca | -471.678579 | correct |
| atq70401.1-pqq-la_model | canonical / Ca | -485.096244 | correct |
| bbl57595.1-pqq-la_model | canonical / Ca | -485.074192 | correct |
| o24759-pqq-la_model | canonical / Ca | -473.773643 | correct |
| p12293-pqq-la_model | canonical / Ca | -486.677138 | correct |
| p15279-pqq-la_model | canonical / Ca | -486.805165 | correct |
| p16027-pqq-la_model | canonical / Ca | -487.891880 | correct |
| p38539-pqq-la_model | canonical / Ca | -473.436320 | correct |
| q4w6g0-pqq-la_model | canonical / Ca | -471.545681 | correct |
| q60ar6-pqq-la_model | canonical / Ca | -473.799941 | correct |
| q88jh5-pqq-la_model | canonical / Ca | -463.709057 | correct |
| q8gr64-pqq-la_model | canonical / Ca | -479.073805 | correct |
| q9l935-pqq-la_model | canonical / Ca | -486.490339 | correct |
| q9z4j7-pqq-la_model | canonical / Ca | -464.235878 | correct |
| 1H4I | crystal / Ca | -505.420107 | correct |
| 4MAE | crystal / La | -444.304425 | correct |
| 1KB0 | crystal / Ca | -482.166486 | correct |

The1KB0 mathematical R is−405482.115532438 and operational R is
−405482.166486368, reflecting the pre-existing0.1kcal origin rule. Both are
Ca-supported. No score variant was selected by classification.

## Existing hard-source replay after the freeze

| Noncanonical source | Operational R +405000 | Frozen-own-reference call |
|---|---:|---|
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1 | -450.716602 | La-supported |
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-3 | -455.470172 | La-supported |
| a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4 | -454.844132 | La-supported |
| a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-4 | -439.808385 | La-supported |

These are the same four previously scored noncanonical pilot sources. No new
fold energies or225-fold comparison occurred here. The pilot raw minimum La−Ca
gap remains8.238885 versus−4.290011 under original adaptive; A0A3 pair spread
narrows8.201646→4.753570, while A0AC stays15.035748. See [pilot report](REPORT.md)
for all source-pair, component-work and transferred-band results.

## Actual execution, coverage and cost

All48 bounded searches were admitted, all48 new cross-MACE evaluations and
all192 GFN2 single points completed. The remaining four sources were reused
exactly from the pilot. No q0 force calculations, DFT, retries, new chemistry or
fallback scores were introduced. All28 stay available.

| Stage | Job | Wall s | CPUs | New work |
|---|---:|---:|---:|---|
| Native adaptive searches |1210465|657|32|2275 MACE evaluations|
| Cross-metal scoring |1210505|36|32|48 MACE single points|
| Vacuum/ALPB GFN2 |1210506|501|64|192 single points|
| Collection |1210507|16|1|No chemistry|

New canonical-stage allocation: **54256 core-seconds / 693 requested GPU-seconds**.
Both GPU stages used oneH200 and200000MiB host memory; solvent used64CPUs/128GiB
withoutGPU. Including the entire earlier eight-case pilot (which also includes
four noncanonical folds), this development branch used69,331core-seconds and
820requestedGPU-seconds. Historical reused origins, local preparation and tests
are additional/unmetered. These are not end-to-end production scanner timings.

## Runtime limitation: numerical boundary settling

Q89GY2 La converged after134 SLSQP iterations /1304 function evaluations, taking
364.093s. The last30 saved trials had exactly equal printed MACE energy and
angular differences below7e−7rad, at a0.800000000000008Å heavy-displacement bound.
Bond/cap/overlap checks passed. This is consistent with numerical line-search
settling at a constraint boundary, not evidence for a different chemical basin.
No tolerance or acceptance criterion changed. It remains a significant cost
outlier; C5B120 La/Ca took39.34/26.45s, while most endpoints were much shorter.
Nineteen of48 new endpoints reached a trust-domain boundary; these are finite
constrained proposals, not proven full/composite minima.

Any future optimizer repair should align numerical convergence with the already
frozen physical admissibility and meaningful energy precision uniformly, under
a new engineering version. This run’s actual attempts/costs must remain intact.

## Scope and recommendation

Pursue the frozen-reference fold challenge next, after parent coordination;
no225 expansion has been launched. Preserve production/defaults. Membership is
the source-backed union across ten saved folds per protein, not a validated
future three-fold preparation policy. Changed membership affects composition,
cavity, native response and selected coordinates together; no unique second-shell
or hydrogen-bond mechanism is established. Canonical labels also correlate with
composition, so this panel alone cannot prove broad added discriminatory value.

The implemented energy remains E_OMOL,vac + E_GFN2,ALPB − E_GFN2,vac. Both metals
score the same {origin,adaptive_Ca,adaptive_La} physical pool, with unchanged
source states/caps/waters/protonation. No solvent-force optimization, scaffoldFF
energy, entropy, extra correction or whole-protein inference is included.

## Artifacts and operations

- [Frozen continuation plan](CANONICAL_PLAN.md).
- [Runnable validators/reference replay](CANONICAL_COMMANDS.md).
- [Compact result](CANONICAL_RESULT.json) and [actual allocations](CANONICAL_SACCT_v1.txt).
- Full frozen reference: `workspaces/union_adaptive_20260923/CANONICAL_REFERENCE_v1.json`.
- Physical/optimizer audit: `workspaces/union_adaptive_20260923/CANONICAL_OPTIMIZER_AUDIT_v1.json`.
- Actual matrix: `workspaces/union_adaptive_20260923/canonical_pool_v1/collection_final.json`.
- [Final real-fixture tests](CANONICAL_TESTS_v2.txt), including actual extrema and deleted-energy rejection.
