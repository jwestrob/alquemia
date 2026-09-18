# Separate charge and kernel changes in the verified GGR self-solvation shift

Native15call validation1201172 confirms the analytical source-only component
within6.4e-13 kcal/mol. Its 2FW0-minus2FVY Ca-minus-La shift is+18.129021kcal.
Next, determine whether that difference is mainly due to projected charges,
inter-source distances, or the native effective Born radii. This is a saved-data
algebraic diagnostic, not a physical perturbation or a changed score.

Use only the complete GGR2FW0 and GGR2FVY source inventories from the verified
primary outputs. Require identical source-atom IDs and formal charges, reorder
by those IDs, and fail rather than invent a mapping. For each metal, evaluate
the verified source-source kernel at all eight combinations of charges,
pairwise-distance matrix and effective-radius vector from the two structures.
These mixed combinations are explicit algebraic counterfactuals, not computed
quantum states or new measurements.

Average each factor's marginal contribution across all six factor orderings
(the symmetric Shapley decomposition). Report all raw combinations and all
three endpoint/paired contributions, including signs. Require their sum to
reproduce the archived18.129021 difference within1e-7kcal/mol and both unmixed
endpoints to reproduce the native source-only values within1e-7. Report charge
sums, largest per-atom charge changes and effective-radius changes alongside.
No biological threshold, fitted coefficient, altered radius, deleted residue,
new solvent/DFT/MACE call or predictor change. Attribution is specific to this
mathematical expression and selected source representation, not unique physics.
