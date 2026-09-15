# Fixed-core PQQ-MDH calibration result

**Protocol:** `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`  
**Verdict:** **CALIBRATABLE**

## Scientific interpretation

The D+2 Asp motif-only baseline classifies all 25/25 controls, and that motif also determines the fixed-core atom count and net charge. Therefore this panel cannot support any claim that DFT adds discrimination beyond the motif/charge confound, even if every preregistered calibration gate passes.

## Frozen gates

| Gate | Result |
|---|---:|
| Valid unchanged pairs | 25/25 |
| Strict separation (L > U) | True |
| Gap G >= 5.0 kcal/mol | True |
| Strict leave-one-out midpoint | 25/25 |
| AUROC | 1.000000 |

U (largest Ca R): -405404.923819 kcal/mol  
L (smallest La R): -405396.320887 kcal/mol  
G = L-U: 8.602932 kcal/mol  
Released midpoint T_R (only if every gate passes): -405400.62235271255

Centered U_S (largest Ca S): 14.857129 kcal/mol  
Centered L_S (smallest La S): 23.460061 kcal/mol  
Released centered midpoint T_S: 19.158595204295125  
The gap and every classification decision are identical in R and S; S is only a constant-shifted display scale.

## Scores

| Panel ID | Class | D+2 Asp | E(La), Eh | E(Ca), Eh | R, kcal/mol | S (aquo gauge), kcal/mol |
|---|---|---:|---:|---:|---:|---:|
| a0a3f2yly8-pqq-la_model | La | True | -2536.124449548277 | -3182.164609392479 | -405396.320887 | 23.460061 |
| a0acd6b9f2-pqq-la_model | La | True | -2536.132448594577 | -3182.147533525803 | -405380.586141 | 39.194807 |
| c5atj3-pqq-la_model | La | True | -2535.681705566840 | -3181.698676297841 | -405381.769498 | 38.011449 |
| c5axv8-pqq-la_model | La | True | -2536.039032979030 | -3182.068409048596 | -405389.553966 | 30.226982 |
| c5b120-pqq-la_model | La | True | -2536.014569267596 | -3182.039960411339 | -405387.053387 | 32.727561 |
| i0jwn7-pqq-la_model | La | True | -2535.998285328168 | -3182.014865191739 | -405381.524225 | 38.256722 |
| mmol_1770-pqq-la_model | La | True | -2536.113625900349 | -3182.140435153774 | -405387.943265 | 31.837683 |
| mmol_2048-pqq-la_model | La | True | -2535.652060183824 | -3181.674695049542 | -405385.323797 | 34.457151 |
| q88jh0-pqq-la_model | La | True | -2535.670813216289 | -3181.697457169118 | -405387.839537 | 31.941411 |
| q89gy2-pqq-la_model | La | True | -2536.099441956949 | -3182.122583121985 | -405385.641504 | 34.139444 |
| q92wy9-pqq-la_model | La | True | -2536.047162111444 | -3182.069151556201 | -405384.918789 | 34.862159 |
| a8r3s4-pqq-la_model | Ca | False | -2307.221739442743 | -2953.289004339939 | -405413.329564 | 6.451384 |
| atq70401.1-pqq-la_model | Ca | False | -2307.459959775945 | -2953.532355126346 | -405416.548972 | 3.231976 |
| bbl57595.1-pqq-la_model | Ca | False | -2307.405597819419 | -2953.472230162447 | -405412.932631 | 6.848317 |
| o24759-pqq-la_model | Ca | False | -2307.342380513589 | -2953.415601226655 | -405417.066895 | 2.714053 |
| p12293-pqq-la_model | Ca | False | -2307.448060803240 | -2953.517484612976 | -405414.684302 | 5.096646 |
| p15279-pqq-la_model | Ca | False | -2307.488894192972 | -2953.560615376699 | -405416.125926 | 3.655022 |
| p16027-pqq-la_model | Ca | False | -2307.436394131634 | -2953.506681693916 | -405415.226315 | 4.554633 |
| p38539-pqq-la_model | Ca | False | -2307.473491720060 | -2953.528662530320 | -405405.740410 | 14.040538 |
| q4w6g0-pqq-la_model | Ca | False | -2197.973912007076 | -2844.027781496087 | -405404.923819 | 14.857129 |
| q60ar6-pqq-la_model | Ca | False | -2307.214693465895 | -2953.277332483763 | -405410.426781 | 9.354167 |
| q88jh5-pqq-la_model | Ca | False | -2307.088362633949 | -2953.147587295830 | -405408.284240 | 11.496707 |
| q8gr64-pqq-la_model | Ca | False | -2198.040985286417 | -2844.117766486728 | -405419.301135 | 0.479813 |
| q9l935-pqq-la_model | Ca | False | -2307.461311336739 | -2953.523479200126 | -405410.131127 | 9.649821 |
| q9z4j7-pqq-la_model | Ca | False | -2307.346992498117 | -2953.411529480632 | -405411.617772 | 8.163176 |

R is the primary preregistered contrast. S differs only by the pinned constant aquo gauge; it was not used for any gate.
