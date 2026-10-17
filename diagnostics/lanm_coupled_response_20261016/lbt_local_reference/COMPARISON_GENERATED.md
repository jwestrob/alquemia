# Local LBT48 electronic reference comparison

Status: complete.

Same-metal source strain S = E(7CCN) − E(7CCO), kcal/mol. No new molecular calls.

| Method | La S | Dy S | S(Dy) − S(La) |
|---|---:|---:|---:|
| MACE | 77.640037 | 3.950156 | -73.689881 |
| GFN2 | -14.563185 | -22.508206 | -7.945021 |
| PBE0 | 69.911132 | 3.584798 | -66.326335 |

| Method minus PBE0 | La strain error | Dy strain error | Differential error |
|---|---:|---:|---:|
| MACE | 7.728904 | 0.365358 | -7.363546 |
| GFN2 | -84.474317 | -26.093003 | 58.381314 |

| Reference state | Availability | All printed residuals pass |
|---|---|---|
| 7CCO_La | available | False |
| 7CCO_Dy | available | False |
| 7CCN_La | available | False |
| 7CCN_Dy | available | False |

PBE0 admits native convergence plus the existing state/gradient gates. Printed residual failures remain visible in RESULT.json; producing a gradient does not qualify its derivatives or electronic stability.

MACE uses physical La1/Dy6 metadata. GFN2 and the frozen-f PBE0 Dy model use effective closed valence states; spin–orbit physics is unqualified. This finite comparison does not change a default, establish affinity, or select a new search.

Exact parents, all four availability records, raw receipts and continuous method-minus-PBE0 differences are retained in RESULT.json.
