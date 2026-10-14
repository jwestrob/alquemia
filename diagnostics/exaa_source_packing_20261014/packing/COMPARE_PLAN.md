# Fixed-pool packing comparison

Read only the declared two-source preparation, the existing inconclusive-review
scores and the new12-cell packing collection. Join by exact original case_key
and record_id, and verify each old MACE request/XYZ/native-solvent input against
the preparation's original matrix. New cells must match the prepared new
XYZ/state/charge and actual MACE/native receipts. No energy calculators run.

For each context separately compute each metal's candidate work relative to its
own origin. Report MACE, total frozen-charge solvent and its four native
components, total work, same-candidate Ca-minus-La differential work, and both
mathematical/operational selected works. Mathematical selection uses ordered
origin/adaptive_Ca/adaptive_La minima; operational selection retains origin
unless improvement is strictly greater than0.1kcal/mol. Compare expanded-minus-
small **works**, never absolute total energies across context sizes. Preserve
all candidates and both sources. A missing cell prevents complete-pool selection;
any individually available origin-relative work remains explicitly partial.

Reproduce the actual old saved per-candidate works and both selection-policy
expectations within the existing1e-8kcal/mol algebra tolerance before reporting
the new comparison. Also check new complete results against the matrix collector.
These are real saved-output consistency checks, not invented molecular fixtures.
Report unavailable cells with reasons. No calibration, class call, affinity,
new optimization or production change. All numerical comparison/testing is
Slurm-only, one CPU. Root owns submission and any resulting scientific decision.
