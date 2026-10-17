# Original four-state hybrid response

Three searches are admitted but nonconverged. The fourth is unavailable after
an embedded-GFN SCC failure. No complete paired response comparison is available.

| Source / metal | Steps / evaluations | Work (kcal/mol) | Max force (kcal/mol/Å) | Result |
| --- | ---: | ---: | ---: | --- |
| 7CCO / La | 80 / 80 | −28.244859 | 3.317400 | Admitted, step limit |
| 7CCO / Dy | 80 / 80 | −28.397116 | 1.788255 | Admitted, step limit |
| 7CCN / La | 2 / 3 | unavailable | unavailable | Embedded SCC failure |
| 7CCN / Dy | 80 / 83 | −31.300723 | 7.717793 | Admitted, step limit |

The failed state's partial history is preserved in SUMMARY.json; its last
successful point is not substituted for the required endpoint. The La/7CCO
80-step count includes the 12-step scout. Original budgets and convergence gates
are unchanged. Terminal work is predominantly MM lowering and is not evidence
of improved metal discrimination.

Jobs13017/13018/13019/13020 used315/361/87/385 scheduler seconds, each12CPUs
and oneH200. Scout13016 used95seconds on the same resources separately. Actual
application/component calls and raw result hashes are in SUMMARY.json.
