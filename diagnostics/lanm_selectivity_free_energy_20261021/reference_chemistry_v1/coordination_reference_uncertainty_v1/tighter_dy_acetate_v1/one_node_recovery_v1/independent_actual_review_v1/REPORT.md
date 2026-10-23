# Tighter Dy-acetate one-node recovery: actual review

**Passed saved admission and the unchanged 0.25 kcal/mol diagnostic.** Ninety-eight independent checks bind all three raw roles and canonical/tight comparison tokens. Two original endpoints18516/18517 are reused exactly; only the missing ligand-with-ghost-ion endpoint was calculated as18693. No native calls or retries occurred in this review.

| CP interaction (kcal/mol) | Canonical | Original tight | Tighter |
|---|---:|---:|---:|
| CCSD | −513.932184 | −513.744599 | −513.892831 |
| Triples | −3.426198 | −3.170086 | −3.356843 |
| Total | −517.358382 | −516.914685 | −517.249673 |
| Total minus canonical | 0 | +0.443698 | **+0.108710** |

Tighter-minus-tight is −0.334988 kcal/mol: −0.148232 from CCSD and −0.186757 from triples. The remaining tighter deviations are +0.039354 and +0.069355 kcal/mol respectively. Independent Decimal sums agree with the collector within3.1e−11 kcal/mol. RHF energies are printed-identical to canonical for all three roles.

## Native and scientific scope

The replacement input is byte-identical to canceled18518: RHF DLPNO-CCSD(T1), TightPNO with the declared tighter cutoffs, def2-TZVP/explicit Dy basis, fixed Aux-C, NoFrozenCore, effective singlet, MaxCore3000 and eight MPI ranks. The runtime used one node, eight ranks and32000 MiB. The full-parent189-AO inventory and718 Aux-C functions match across roles. Complex/ion/ligand contain40/8/32 electrons, all correlated; ECP55 is assigned only when Dy is real. The ghost-Dy ligand role has no ECP.

CC terminal residuals are8e−9/3e−9/8e−9 against1e−8. Iterative-triples terminal energy changes are5.32e−9/0/0 against1e−8. T1 diagnostics are0.020818158/0.001482672/0.014995481. All actual triples settings and stage-specific local thresholds are retained; the stronger/weak-triples thresholds must not be confused with the changed pair/PNO/MKN values.

All three SCFs stopped by the native Energy Check. Their printed density and DIIS tests fail; the ligand's printed orbital-rotation test also fails. These warnings remain in the receipt, alongside the effective-singlet versus physical Dy-sextet/frozen-f limitation. The result establishes this finite truncation diagnostic relative to canonical under the same provisional Hamiltonian; it does not qualify physical accuracy, the full18 controls or cluster release.

## Preservation and costs

Canceled18518 and the original unavailable2/3 collector remain intact. Its59512 allocated CPU-seconds are retained separately. Replacement18693 completed in2169 scheduler seconds ×8CPUs =17352 CPU-seconds; collector18694 completed in3seconds ×1CPU. Selected three-role native elapsed×rank accounting is51739.044 CPU-seconds, including reused work, and is distinct from scheduler accounting. No canceled-attempt energy is substituted into the CP cycle.
