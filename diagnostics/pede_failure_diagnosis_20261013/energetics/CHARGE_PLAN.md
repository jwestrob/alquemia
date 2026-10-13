# Saved-origin charge and solvent addendum

Question: is the source-dependent solvent contrast accompanied by different
metal/cofactor/donor charge distributions or visible SCC convergence behavior?
Read all 36 actual origin points (nine sources × two metals × two media), both
saved backends. No calculator calls, new initialization, optimization or threshold
changes. This is a descriptive follow-up to actual Slurm analysis 3341.

Preserve each metal's origin ALPB-minus-vacuum contribution, and decompose the
already-identified two-Ca sample-1 minus sample-3 solvent contrast by metal. Other
sources remain visible; these extrema do not select a predictive population.

Read raw charges/scientific JSON and require exact agreement with admitted point
records. Use the full preparation's explicit QM indices, including the three
fixed PQQ hydrogens absent from physical motion coordinates. Every atom must map,
with element and origin-coordinate agreement (1e-7 Å serialization tolerance).
Report complete PQQ, each mapped residue, metal and synthetic caps separately.
Do not treat capped fragment charges as whole-protein residue charges.

Require actual successful receipts and saved convergence; parse only SCC-table
row shapes (excluding orbital population lines), retaining final residuals and
raw rows. Fractional occupations between 0.001 and 1.999 e are descriptive at the
existing 300 K electronic temperature. They do not assign oxidation states or
prove an electronic branch switch. Backend agreement does not independently
validate the electronic model or identify the physical solution.

`CHARGES.py --config <actual CONFIG.json>` runs only in an actual one-CPU Slurm
allocation. Outputs are `RESULT.json`, `ORIGIN_CHARGES.tsv`, `ATOM_CHARGES.tsv`,
`ORIGIN_SOLVENT.tsv` and a measured execution receipt. Missing/changed files fail
visibly; no fabricated charge file or zero substitution. Root owns submission.
