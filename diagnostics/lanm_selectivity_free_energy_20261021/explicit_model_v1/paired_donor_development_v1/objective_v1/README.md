# Saved-only donor objective

`OBJECTIVE.py` implements the unchanged `../DESIGN.md` objective. Root accepted the existing28 consistent ECP-PBE0-D4 rows as provisional development targets, preserving the narrow central hydrated-amide translation qualification and method/state limitations. This is not universal reference-force accuracy, a bound water candidate, a fitter, or an execution release.

## Interface

Call `evaluate(data, model, body, water_baseline)`:

- `data['rows']`: the exact28 admitted training rows, including references, geometry, covalent graph and rigid bodies. Existing training authority is `matched_dry_v1/training_v1/RESULT.json`, SHA `99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb`.
- `model`: exactly those28 IDs, each with `energy` in kcal/mol and `forces` in kcal/mol/Å. Forces have shape atom-count×3. No missing-cell averaging or extra outcomes are accepted.
- `body`: the existing `matched_dry_v1/OBJECTIVE.py` helper, SHA `bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819`. Its category/body/translation/rotation weighting and centroid torque/1Å are reused unchanged.
- `water_baseline`: explicit E/forces for `water_central_{La,Dy}` and `aquo8_central_{La,Dy}` from the future frozen water candidate. The module never picks the historical baseline implicitly. Exact equality is required: donor changes cannot alter these controls or the aquo anchor. The additional four radial-water records and full eight-water numerical/development qualification remain external authorities.

The caller must admit full source, native and candidate provenance before calling. The objective itself performs no file access, native parsing, molecular calls or parameter changes. It validates complete counts, finite arrays, chemical states, paired geometry/body consistency and supplied water invariance.

## Arithmetic and output

Four loss blocks each receive weight0.25: energy, La, Dy and Dy−La rigid response. Each force block averages all12 donor pairs. The energy block uses0.50 mean square of12 aquo-referenced exchanges,0.30 mean square of3 wet-minus-matched-dry double differences and0.20 mean square of12 individual radial works, equally including both metals. Six paired radial-work errors are controls outside the loss.

Energy errors are model−reference; metal contrasts are Dy−La; radial works are endpoint−center. Each fixture must satisfy rigid RMS≤2 and maximum body component≤5 kcal/mol/Å; every reported donor energy error must have magnitude≤1 kcal/mol. Panel averages cannot hide a failed fixture. Atomic/internal-force residuals are retained but excluded from fitting.

`donor_development_checks_passed` covers the24 donor rows. `development_checks_passed` additionally requires the supplied four water-center controls to pass; its explicit scope remains this28-row panel. It does not assert the external eight-water qualification, transfer, physical accuracy or production readiness. Exact four-center water invariance is required before a score is returned.

## Actual saved validation

`VALIDATE.py` admitted the pinned original28 data and body helper and passed72 checks in0.997 seconds on one local thread, with zero native programs or molecular calls. Checks include independent signed cycle arithmetic, element-energy-zero and common-composition cancellations, analytical weighting cases, paired-force cancellation without hiding individual errors, missing/extra/nonfinite cells, paired geometry mismatch and changed-water rejection.

The historical baseline score is29.555283241651008. This is only an arithmetic replay of the original baseline; it is not the future passing water candidate or a newly fitted result. No reserved validation outcome was accessed. An initial repository-path error occurred before data admission; its source and failure receipt remain in `validation_v1`. The corrected actual result is in `validation_v2`.

| Artifact | SHA256 |
|---|---|
| `OBJECTIVE.py` | `86e85fcb0537bdf5b459874eeda0c8bb68aead9758c712014b18b2db94757974` |
| `VALIDATE.py` | `bdcb6efc0e644d0b9a3c445bb74fe255b45d38e1a87b36ce09c5ff403b3cca00` |
| Workspace `objective_v1/validation_v2/RESULT.json` | `86efece9bd9aa4ac633f0b79c6c1350b46154fe9f2d43d3243da64e66c535694` |
| Original `../DESIGN.md`, unchanged | `343efe3eed4a124e6dea4fc64d37821d5b6beb78ce1705f4ae7a11beb0425a31` |

Diagnostic prefix: `diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_development_v1/`. Workspace prefix uses `workspaces` in place of `diagnostics`. No runnable molecular configuration or candidate binding exists in this bundle.
