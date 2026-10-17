# Saved Mex DFT origin assembly

Question: can the two already admitted Mex EF3 PBE0 origins be assembled with
their exact existing MM ledger and fixed-FF-charge global ALPB into traceable
whole-scaffold energies and analytic gradients?

This is a saved-data calculation on one Slurm CPU, with zero molecular calls,
new geometries or optimization. Both La and Dy remain in the declared denominator.
Reuse job12686's pinned `MEX_REFERENCE_RESPONSE.inputs/gradients` and its exact
reference/static parents. The assembly is

`PBE0_embedded(Q209; shifted1840PC) + MM_ledger + fixedFF_ALPB(R2049)`.

MACE vacuum and the GFN2 embedding difference are absent from this new energy
ledger. No additional QM–MM Coulomb, C4, full-protein FF energy or cap population
projection is added. Preserve all 2,049 physical atoms, three Nd spectators,
waters, charge states, cap derivatives and the existing external-energy convention.

Replay the saved native DFT energy/gradient and full point-charge derivative;
verify MM component closure (existing 1e-8 kcal/mol), exact replay of the solvent
component sum, physical ordering and gradient reconstruction (existing
1e-10 absolute / 1e-12 relative tolerances). Keep the original pin checks and
native SCF convergence semantics. Record every printed terminal residual,
ConvCheckMode/ConvForced, warning lines and a separate all-residuals flag; this
flag does not replace the archived admission rule.

Write each metal's native/component energies, full physical component gradients,
total gradient and negative-gradient forces, source/state/receipt pins and replay
residuals. Report translation sums as observations without a new pass threshold.
Unavailable states stay unavailable; no fallback, rerun or subset comparison.

Physical-state and DFT total-force consistency remain unqualified. The previous
MACE/GFN finite-difference pass does not qualify this replacement surface. Do not
interpret absolute La/Dy energy offsets as affinity, create a classifier, change
defaults or launch a DFT campaign. This assembly provides auditable saved origins
for a separately declared future test.
