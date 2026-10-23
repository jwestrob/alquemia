# Independent review: eight hydrated-acetamide single points

**No launch blocker identified.** The saved-input review passed 815 checks and 48 unique file pins in 2.08 seconds on one local thread. It made zero native calls and zero submissions.

- The eight actual XYZ/input files preserve the original PBE0-D4 scientific recipe, removing only EnGrad and changing allocation controls. Metal and seven waters stay fixed; the nine ligand atoms translate rigidly by the declared ±0.010/±0.005 Å. Maximum translation error is 1.56×10⁻¹³ Å; maximum internal-distance change is 8.88×10⁻¹⁶ Å. La/Dy geometries match exactly.
- Independent inventories confirm 31 real atoms, 110 explicit electrons, 502 orbital basis functions, charge +3 and effective singlet with La46/Dy55 frozen cores. Physical La singlet/Dy sextet labels remain distinct.
- Both saved central outputs replay under the strict native parser. Independent ligand-force projections give La 14.852538269 and Dy 31.729427236 kcal/mol/Å. Positive displacement is outward and the derivative is minus the summed ligand force projection.
- Signed polynomial controls verify electronic/D4/final derivatives and the Dy-minus-La sign. Source review confirms separate individual and paired two-step/analytic 0.1 kcal/mol/Å checks, conditional smaller-step summaries and explicit energy-stage change flags. Missing results retain all eight unavailable slots without a derivative conclusion.
- The four-rank, one-thread-per-rank launcher uses standard_expand, default memory/time and no requeue. Native admission, source/input replay and actual allocation-derived MaxCore remain enforced.

The old coarse-screen failure remains preserved. This review does not establish derivative agreement before the new results arrive, physical reference accuracy, or eligibility of other donor states.

## Exact evidence

| Artifact | SHA256 |
|---|---|
| Independent RESULT.json | `3e685421cf9836ebe2d8ee3e37f6bd25fc0da592c69d4633e62e29da8f3eca09` |
| Independent REVIEW.py | `0cbc0a63d1f477ec8eaff1f029b24cbf340c9b99419db5cf27d9aee9d09b3455` |
| Author REFERENCE.py | `7de9c17aa9198ba07c66c3f934390fc2a67d94a77860eca8ba4c0378b25ad677` |
| SPEC.md | `f5de75cfa8b006a07e93267f6ebcdeef1d02e37a31d544872656698bc92112d4` |
| prepared_v2/PANEL.json | `38b42c1020c47b12512076d59d66e299d4c8eb7f14754f8cfe5f71bec502e4cc` |
| READY.json | `a88e0f63cfa69ab19ab986b69af05c076bf84fb9829a0b027fe8d92a9157dff3` |

Result location: `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_amide_force_independent_review_v1/RESULT.json`.

Executed command:

```bash
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_amide_force_independent_review_v1/REVIEW.py
```
