# All PQQ reference mappings are ready

2026-09-20. **28/28 contexts supported**, covering all 25 canonical references and consumed crystal controls 1H4I, 4MAE and 1KB0. This prepares a common accommodation rule for later testing across both reference classes. No energies, optimization, folds, protonation or submissions ran; the baseline and ongoing DFT campaign are unchanged.

## Common physical rule

Every case activates its original anchor-Glu chi3. Twelve cases also activate the original extra-Asp chi2; the other sixteen have ALA (10), SER (3) or THR (3) at that homologous role and keep it fixed. This follows actual chemistry, without consulting labels or energy scores. All other coordinates are zero. Source/cap mappings use the existing torsion-control Kinematics and the first nonlinear plan unchanged.

| Case | Anchor-Glu mode | Extra homolog | Extra active mode |
|---|---|---|---|
| a0a3f2yly8-pqq-la_model | A/197/chi3 | ASP | A/329/chi2 |
| a0acd6b9f2-pqq-la_model | A/204/chi3 | ASP | A/347/chi2 |
| c5atj3-pqq-la_model | A/191/chi3 | ASP | A/319/chi2 |
| c5axv8-pqq-la_model | A/198/chi3 | ASP | A/319/chi2 |
| c5b120-pqq-la_model | A/192/chi3 | ASP | A/320/chi2 |
| i0jwn7-pqq-la_model | A/206/chi3 | ASP | A/335/chi2 |
| mmol_1770-pqq-la_model | A/200/chi3 | ASP | A/362/chi2 |
| mmol_2048-pqq-la_model | A/195/chi3 | ASP | A/357/chi2 |
| q88jh0-pqq-la_model | A/199/chi3 | ASP | A/325/chi2 |
| q89gy2-pqq-la_model | A/194/chi3 | ASP | A/322/chi2 |
| q92wy9-pqq-la_model | A/193/chi3 | ASP | A/321/chi2 |
| a8r3s4-pqq-la_model | A/213/chi3 | SER | fixed |
| atq70401.1-pqq-la_model | A/204/chi3 | ALA | fixed |
| bbl57595.1-pqq-la_model | A/206/chi3 | ALA | fixed |
| o24759-pqq-la_model | A/211/chi3 | ALA | fixed |
| p12293-pqq-la_model | A/209/chi3 | ALA | fixed |
| p15279-pqq-la_model | A/204/chi3 | ALA | fixed |
| p16027-pqq-la_model | A/204/chi3 | ALA | fixed |
| p38539-pqq-la_model | A/173/chi3 | ALA | fixed |
| q4w6g0-pqq-la_model | A/207/chi3 | THR | fixed |
| q60ar6-pqq-la_model | A/205/chi3 | ALA | fixed |
| q88jh5-pqq-la_model | A/221/chi3 | SER | fixed |
| q8gr64-pqq-la_model | A/195/chi3 | THR | fixed |
| q9l935-pqq-la_model | A/199/chi3 | ALA | fixed |
| q9z4j7-pqq-la_model | A/213/chi3 | SER | fixed |
| 1H4I | A/177/chi3 | ALA | fixed |
| 4MAE | A/172/chi3 | ASP | A/301/chi2 |
| 1KB0 | A/185/chi3 | THR | fixed |

## Actual geometry/state checks

- All 56 endpoint origins retain exact archived XYZ pins and unchanged context composition, charge and singlet multiplicity. Original-core and context electron parity checks pass. Ca/La nonmetal coordinates and physical mappings are identical within each case; endpoint charges differ by exactly one.
- All 112 q=0 core/context replays pass the existing 1e-10 Angstrom criterion. In each representation 24/56 are bitwise identical; the remaining differences are at most 4.440892098500626e-16 Angstrom. Original XYZ files are not rewritten or rounded.
- All active terminal groups retain both carboxylate oxygen atoms. Existing analytic source/cap Jacobian checks pass: maximum derivative discrepancy 5.980210104894468e-9 Angstrom/unit, maximum bond discrepancy 2.4424906541753444e-15 Angstrom, fixed atoms unchanged.
- Active modes have the same IDs and index mappings for both metals. All unused modes, including metal translations and other donor torsions, remain zero. No geometry-only result implies a stationary energy minimum or better classification.

## Artifacts and operation

`workspaces/accommodation_nonlinear_20260920/reference_preparation_v1/READY.json` pins the supported output. The compatible torsion-style request is **`reference_preparation_v1/ready_v2/design.json`**, with `cases`, `source`, `preparation`, `modes`, `maps`, `origins` and model/software/ORCA/agreement pins. Per-case source/cap JSON and checks are under that same directory; `result.json` preserves the full denominator. Original-core parity checks are in `original_core_state_checks.json`.

The first attempt retained all three supported crystal maps but failed the 25 canonical cases because their old output paths were relative to the parent manifest. Its results and builder source remain at `reference_preparation_v1/`. The corrected adapter resolves those paths against each unchanged parent manifest and verifies the original SHA256; it changes no scientific input. Corrected preparation took 44.761992918 seconds, first attempt 12.254693497 seconds, with zero molecular calls.

Replay the geometry preparation into a new directory from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/accommodation_nonlinear_20260920/prepare_reference_maps.py \
  --root . \
  --output workspaces/accommodation_nonlinear_20260920/reference_preparation_replay_v1 \
  --agreement diagnostics/accommodation_nonlinear_20260920/REFERENCE_PREPARATION_PLAN.md \
  --parent-plan diagnostics/accommodation_nonlinear_20260920/PLAN.md
```

Root owns review and the decision to run a reference energy experiment. This preparation neither submits nor authorizes those calculations. The common rule still needs actual energy and discriminatory-utility testing; all references are already consumed evidence.
