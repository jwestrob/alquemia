# 4MAE: second real environmental response, not a class decision

Six/six endpoints completed: four embedded A/B plus two isolated baselines.
The same chemistry-only selection rule chooses Thr154 HG1, rotated+10degrees.
The electronic double difference(LaB-LaA)-(CaB-CaA) is+0.254914852kcal/mol.
The modeled motion thus relatively favors Ca, despite4MAE being a consumed
La-associated PQQ control. There is no expectation that the sign of an arbitrary
local perturbation equals a protein's preference. Do not interpret this as a
classifier failure/success or use it to recalibrate anything.

This extends the finite native response harness to a second actual source.
The conditional dry80-atom core and8774pointcharges preserve archived15P603 and
water exclusions; originalOH1.18042A remains. Coordinates, microstate and charges
were fixed. No equilibrium sample, complete experimental site or hydration
thermodynamics is represented. Exact outputs retain native gradients and isolated
endpoints for compatible future comparisons.

Numerical qualification remains unestablished: the1H4I rigid gate failed and
4MAE directional/rigid/refinement checks were not run. The intended ML comparison
is unavailable without exact model weights/interface. No cross-protein absolute
energy subtraction, biological classification, or hybrid relaxation is reported.

Job1219487:936s on344allocatedCPUs,321984CPU-s,zeroGPU;collector1219488:2s/1CPU.
Four real preparation checks including six-cell dry-run passed. Molecular workers
cumulative:48outputs,1186112allocatedCPU-s,zeroGPU, including numericaldiagnostics.
These are development costs, not routine inference throughput or utilization.

Next: review the separately identified52-atom1F6S nonPQQ representation and build
an explicit additive interaction ledger from the now-available exact-source
protein mechanics and PQQcrossLJ. Keep numerical/ML/fullhybrid limitations visible.
