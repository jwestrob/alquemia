# GGR native OMOL readout diagnostic

All replay/accounting checks pass: **True**.
Connected minus extended R: **-26.191960550 kcal/mol**. No score has changed.

| Native readout group | Extended R | Connected R | Contribution to shift |
|---|---:|---:|---:|
| metal | -405436.519091179 | -405438.794685744 | -2.275594565 |
| shared_source_atoms | 18.605936995 | 24.246435687 | 5.640498692 |
| unique_source_atoms | 0.000000000 | 6.964273698 | 6.964273698 |
| synthetic_caps | 3.165736990 | 2.832504833 | -0.333232157 |
| embedding_graph_term | -11.143976457 | -47.331882674 | -36.187906218 |

Values are kcal/mol. Each graph term is retained separately; source hydrogens remain source atoms, not caps.

The explicit linear embedding readout has the form sum(a_element) + N*b(charge, spin). Its nonmetal Ca-minus-La term is -0.682790683351 kcal/mol per atom. The extra 53 atoms predict -36.187906218 kcal/mol of embedding shift; observed minus predicted is -7.11e-15. This term does not depend on geometry. The remaining readouts still depend on charge-conditioned features and local context.

Atomic readouts are learned bookkeeping, not a unique physical energy partition. Identifying this extensive term does not establish that simply deleting it gives correct affinity. Both original scores remain in use; no corrected score is produced.
