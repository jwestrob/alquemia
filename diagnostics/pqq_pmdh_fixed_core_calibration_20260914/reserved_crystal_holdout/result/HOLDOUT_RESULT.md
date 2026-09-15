# Reserved PQQ-MDH crystal transfer result

**Protocol:** `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`
**Primary verdict:** **PASS**

The primary rule was frozen before these energies were read: 1H4I must have S <= 14.857129203 kcal/mol and 4MAE must have S >= 23.460061206 kcal/mol. The open interval is indeterminate and fails the primary transfer test.

| PDB | Role | Expected | S, kcal/mol | Frozen-band call | Primary pass |
|---|---|---|---:|---|---:|
| 1H4I | primary | Ca-supported | 7.517843566 | Ca-supported | True |
| 4MAE | primary | Ln-supported | 38.089733042 | Ln-supported | True |

6OC6 was not run; it remains a nonindependent optional secondary geometry check.
4MAE was evaluated after explicit removal of coordinating 15P603/OXT without replacement; it is a dry fixed-coordinate transfer test with a ligand vacancy.
No threshold or band was fit, shifted, or widened using holdout results.
