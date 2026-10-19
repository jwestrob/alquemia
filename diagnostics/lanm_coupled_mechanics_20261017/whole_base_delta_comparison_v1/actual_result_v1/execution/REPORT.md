# Whole-base regional-delta: consumed 7CCO comparison

Eight declared states; four development cells and four conditional 7CCN challenge slots.
Work is after minus before. Energy units are kcal/mol; force units are kcal/mol/Angstrom.

## Deformation work

| Quantity | La | Dy |
|---|---:|---:|
| GFN_full | 5.016670 | 4.845637 |
| MACE_region | -10.371217 | -11.618549 |
| negative_GFN_region | 9.843322 | 10.024922 |
| fixed_FF_ALPB | 5.382042 | 5.948114 |
| Candidate | 9.870817 | 9.200125 |
| Base only | 10.398712 | 10.793752 |
| Whole PBE0 + same ALPB | 12.279593 | 11.722656 |
| Whole MACE + same ALPB | 11.290076 | 11.082401 |

## Work discrepancies from the whole-PBE0 reference

| Model minus reference | La | Dy |
|---|---:|---:|
| candidate | -2.408776 | -2.522531 |
| base_only | -1.880882 | -0.928904 |
| whole_MACE | -0.989517 | -0.640255 |

**Development triage: failed.** 7CCN candidate challenge closed.
The frozen rule requires absolute candidate work error <=1 kcal/mol for both metals. It is not a reference error bar or force-accuracy criterion.
Base-only results are reported independently and do not control this candidate gate.

## Force offsets and changes

Errors are model minus PBE0 plus the same ALPB. Change is after minus before; all256 atoms remain represented.

| Metal / model / field | Full256 vector RMS | Phe1-force projection | Chord projection |
|---|---:|---:|---:|
| La / candidate / error_before | 5.428814 | 5.007225 | -2.791638 |
| La / candidate / error_after | 5.291432 | 9.831135 | 3.230741 |
| La / candidate / error_change | 2.028263 | 4.823910 | 6.022379 |
| La / base_only / error_before | 7.249080 | 5.007225 | -2.587036 |
| La / base_only / error_after | 7.163349 | 9.831135 | 2.931978 |
| La / base_only / error_change | 2.258968 | 4.823910 | 5.519015 |
| La / whole_MACE_composite / error_before | 3.511987 | -2.542480 | -1.172226 |
| La / whole_MACE_composite / error_after | 3.516438 | -3.612843 | 1.496902 |
| La / whole_MACE_composite / error_change | 0.530091 | -1.070363 | 2.669128 |
| La / hybrid_composite / error_before | 10.141916 | 47.721979 | 16.289791 |
| La / hybrid_composite / error_after | 8.850041 | 38.955928 | 2.243337 |
| La / hybrid_composite / error_change | 4.463637 | -8.766050 | -14.046454 |
| Dy / candidate / error_before | 5.471512 | 4.995297 | -2.765024 |
| Dy / candidate / error_after | 5.335039 | 9.825720 | 3.252820 |
| Dy / candidate / error_change | 2.032134 | 4.830423 | 6.017844 |
| Dy / base_only / error_before | 7.643386 | 4.995297 | -2.949072 |
| Dy / base_only / error_after | 7.536375 | 9.825720 | 2.847353 |
| Dy / base_only / error_change | 2.306745 | 4.830423 | 5.796425 |
| Dy / whole_MACE_composite / error_before | unavailable | — | — |
| Dy / whole_MACE_composite / error_after | 3.531992 | -4.263123 | 1.397318 |
| Dy / whole_MACE_composite / error_change | unavailable | — | — |
| Dy / hybrid_composite / error_before | 10.162397 | 47.719232 | 16.320952 |
| Dy / hybrid_composite / error_after | unavailable | — | — |
| Dy / hybrid_composite / error_change | unavailable | — | — |

Full vectors, original named terminal-atom summaries and every component force change are preserved in RESULT.json and FORCE_TRANSITIONS.json.

## Denominator and limits

7CCN remains pending and conditionally gated here; no pending reference outcome is read or used.
Before-Dy whole-MACE forces and after-Dy old-hybrid forces/work remain unavailable. No missing array is replaced with zero.
The candidate uses a mechanically embedded vacuum regional correction; it has no mutual high/low polarization. No MM or embedded-GFN term is added.
Conditional terminal chemistry, vacuum-dianion/frozen-f/SCF limitations and nonstationary endpoints persist. No independent candidate derivative or force-accuracy qualification, optimizer, affinity or preference claim follows.
Two endpoints do not determine a Hessian or uniquely separate equilibrium-position and stiffness errors.

This process makes zero molecular calls. Parent allocation/call records remain separate from this saved analysis.
