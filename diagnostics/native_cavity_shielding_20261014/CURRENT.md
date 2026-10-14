# Native solvent-exclusion feasibility; no molecular calculation submitted

The completed midpoint test did not establish better discrimination. Read
../exaa_metal_midpoint_20261014/REPORT.md before designing another search.

Next physical question: does omitted protein exclusion materially change the
fixed core-charge solvent response? An isolated component diagnostic could use
native ALPB on the exact capped core plus deduplicated source atoms with zero
charges, subtracting the original core electrostatic solvation term. Zero charges
on the added atoms explicitly turn off permanent environmental electrostatics;
they are not a proposed protein charge assignment.

The existing native API has a material limitation: charge sites and cavity atoms
share one list. Synthetic cap hydrogens therefore remain descreening and ALPB
shape centers. This is a capped-union exclusion sensitivity, not a source-only
physical cavity or a complete protein-solvent score. No radius edits, charge
redistribution, whole-protein SASA addition, new Coulomb term or calibration are
licensed by this feasibility conclusion. The old failed global-PB protocol stays
closed.

The completed source audit is feasibility/FEASIBILITY.md. Root owns any later
preparation, implementation and Slurm submission. Its proposed12-update scout
uses four consumed origin endpoints, cap controls and numerical checks, with zero
new SCC/MACE/DFT. Freeze actual transformations, input pins and measured memory
placement before execution; the proposal is unrun.
No new model, preparation, executable, molecular result or production change exists
in this namespace yet. All prior jobs are terminal; do not restart them.
