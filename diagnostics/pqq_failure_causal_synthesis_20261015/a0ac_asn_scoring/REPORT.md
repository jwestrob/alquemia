# Asn orientation changes one XoxF1 vacuum pool

**The fixed amide flip changes the ambiguous Ca-conditioned sample4 contrast by +1.460981253638387 kcal/mol toward La.** Both metals had access to the same four geometries; La selects the new flip and Ca retains its previous adaptive geometry. The other nine sources retain their selections. All20 new cells completed, with60 existing cells reused.

This passes the frozen, sign-blind1kcal/mol cutoff for the already-declared matched solvent continuation. It is a vacuum local-energy result, not a changed calibrated classification, affinity or validated biological state.

## What the comparison establishes

For Ca-conditioned4, the flip changes energy relative to origin by +4.685514510164903kcal/mol for Ca and −27.140066810943047 for La. Most of that La stabilization was already captured by its old adaptive candidate; the gain over the full existing pool is only1.460981253638387. Mathematical and operational selections agree.

La-conditioned4 also favors the flip over origin for La (−27.277089306655522kcal/mol), but its old adaptive La candidate remains better. That correctly classified source remains a counterexample to treating the distant Asn oxygen as a sufficient failure diagnosis.

All eight other fixed flips raise both origin energies and lose to the existing pool. Some introduce represented short nonbonded H–H contacts; rejecting those unrelaxed candidates does not reject every relaxed alternate rotamer basin. The changed target's nearest H record is an unchanged covalent bond, not the severe nonbonded contact seen in some peers. No angles, contacts or admission rules were adjusted after inspection.

The saved five-source Ca-conditioned vacuum range decreases17.61986454203725→16.158883288386278kcal/mol. The La-conditioned range stays16.27608960442012. These are correlated structural samples and raw vacuum contrasts. The solvent-inclusive change remains unavailable.

## Actual work and next step

Slurm4617:363elapsedseconds,7allocatedCPUs,2541allocatedCPU-seconds; schedulerTotalCPU30:36.911,batchMaxRSS8300488KiB. The warm worker records20forwards,zero failures,1821.757749174processCPU-seconds,7517480KiB peakRSS andzeroGPU. Separate geometry preparation cost6schedulerseconds on1CPU. All60 archived baseline checks passed before model calls.

[RESULT.json](RESULT.json) preserves exact values, source/implementation pins and the full ten-source denominator. The actual collection is `workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_scoring_v1/execution/FINAL_COLLECTION.json`. Next: execute the declared same-state frozen-charge ALPB continuation for all20 admitted flips, reusing all60 original solvent endpoints. No production, calibration, optimization or additional MACE change is included.
