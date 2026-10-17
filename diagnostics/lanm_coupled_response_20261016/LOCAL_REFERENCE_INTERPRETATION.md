# Interpreting the local electronic reference

This outline applies to the frozen [four-cell plan](LBT_LOCAL_REFERENCE_PLAN.md).
Jacob authorizes DFT as the task default if the evidence establishes better
accuracy. The present test diagnoses a specific local electronic surface;
agreement between approximate methods alone does not establish physical accuracy.

## Compare each metal before their difference

For each method, report S(La) and S(Dy), where S = E(7CCN) − E(7CCO), then
D = S(Dy) − S(La). Form same-metal native energy differences before conversion.
The exact matched cores have 48 atoms and charge −1. Do not mix these quantities
with the whole-peptide, solvent-corrected exchanges or unlike-metal absolute
energies.

The admitted saved comparison from job 12825 is:

| Method | S(La) | S(Dy) | D, kcal/mol |
|---|---:|---:|---:|
| MACE, job 12471 | +77.640037 | +3.950156 | −73.689881 |
| Standalone GFN2, jobs 12488 / 12825 | −14.563185 | −22.508206 | −7.945021 |

Thus the disagreement concerns both metals' source ordering, particularly La,
and is not merely an unlike-element reference offset. PBE0 must be assessed
against all three quantities. Report signed deviations from PBE0 for both
methods; the declared 1 kcal/mol resolution is an interpretation scale, not an
error bar or a fitted acceptance threshold.

## Decisions once all four PBE0 cells are admitted

| PBE0 outcome | Scientific implication and next decision |
|---|---|
| Supports the large MACE-like source penalty and metal-dependent contrast | The local penalty has independent support. Investigate the inherited source geometry/chemistry and constrained response before replacing local electronics. More torsion steps cannot repair fixed bond geometry or differing refinement assumptions. |
| Supports the GFN-like ordering and much smaller contrast, with material MACE deviations for the individual metals | Supports replacing MACE for the next contained local-electronics calculation with PBE0, while retaining the chemical domain and other energy terms. Re-evaluate a complete paired pool before using its selections. This supports a task-specific DFT choice; it does not independently validate the physical Dy state. |
| Matches D but disagrees substantially in one or both individual S values | Differential agreement may reflect cancellation. Own-metal response/geometry selection remains unvalidated; do not choose a method using D alone. |
| Supports neither model, or conclusions depend on a flagged/missing reference cell | Retain method uncertainty. Identify the specific electronic-state, numerical or domain issue before a further finite reference test; no automatic functional sweep, larger optimizer or favorable-subset selection. |

A flatter surface or desired preference sign is never the selection criterion.
Mechanics, solvent and population questions remain separate even if the local
electronic method changes. The old whole/core comparison motivates this local
test but does not prove that a capped 48-atom result transfers quantitatively
to the whole peptide or embedded Mex protein.

## Exact output use and limitations

The GFN result provides `matrix` and `source_strains`; MACE12471 provides the
four `matched_core` rows with `energy_eV`. The new reference's `rows` are keyed
by `7CCO_La`, `7CCN_La`, `7CCO_Dy`, `7CCN_Dy`; use only admitted `available`
`energy_hartree` values and their exact source/task pins. A partial La pair can
inform S(La), but cannot establish D or a four-cell accuracy conclusion.

Retain each reference's `SCF_evidence`, including native convergence mode,
printed residuals, warnings and `all_printed_residuals_pass`. A false residual
flag remains visible under the original admission rule; do not retroactively
drop the cell or call it fully force-qualified. Analytic-gradient availability
does not prove electronic stability or suitability for a new optimizer.

PBE0 uses restricted La ECP46 / Dy ECP55 frozen-f states; GFN uses an effective
singlet, while MACE receives physical multiplicities La 1 / Dy 6. Agreement
between the two frozen-f treatments can share a model limitation. This test
does not establish physical sextet/spin-orbit accuracy, affinity, population,
or a source-robust natural-protein classifier. All comparisons and deviations
must be computed by the saved-result Slurm collector, not on the login node.
