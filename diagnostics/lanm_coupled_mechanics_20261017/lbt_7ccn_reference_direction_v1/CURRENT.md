# Directional reference check complete; strict gate failed

21October2026UTC. Both declared singlepoints15847/15848 are complete and pass
unchanged native admission. Saved collector16582 completed its comparison with
`failed_directional_consistency`: both rows available, no exception, zero calls.
Its SlurmFAILED1:0 is the intended scientific-gate outcome, not a technical fault.

FD force−66.4391321372 versus analyticPBE0−66.4531999104kcal/mol/Å gives
residual+0.0140677732 above the frozen0.0116453200tolerance. MACE projection is
−9.02704252585. The strict gate remains failed; the small sampled inconsistency
does not explain the much larger model/reference directional disagreement.
No complete force qualification, physical-state or preference claim follows.

Read [FINAL_REPORT.md](FINAL_REPORT.md). Nativearchivee86624f0/collectorarchive
7c67a8c2 preserve exact inputs/results/SCF/runtime/accounting and independent
reviewsd87d880f/3c8bc24f. This phase totals518881allocatedCPU-s,zeroGPU,two new
nativeenergycalls. Sourcecd4396c0 and originalgate remain unchanged.

This finite test is CLOSED with its failure retained. No automaticretry,
step/gridchange, centralrerun, furtheroptimizer or othermetal call. The separate published-input reference scout is governed by
[its newer active checkpoint](../lbt_omol_target_reference_v1/CURRENT.md), including
the observed SCNL treatment; older source-only planning is historical.
Delayed callbacks need receipt acknowledgement only, not resubmission.
Minus15847 completion delivery is recorded in DELIVERY_COMPLETION_ACK_15847.json.
Root owns integration and launches.
