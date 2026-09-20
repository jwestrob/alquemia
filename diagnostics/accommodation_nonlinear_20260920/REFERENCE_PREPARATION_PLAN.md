# Common-rule reference mapping preparation

2026-09-20. Root assigned this geometry-only preparation under the existing
[nonlinear plan](PLAN.md), before any reference accommodation energies.

Include all 25 canonical PQQ cases and the consumed 1H4I, 4MAE and 1KB0
crystal controls from the frozen compact-solvation inventory. Reuse their exact
archived original-core and complete-context coordinates, charges, multiplicities,
source graph and cap maps. No hydrogen preparation, atom additions or deletions,
folding, energy evaluation, optimization or submission is permitted here.

Use the same existing torsion-control mapping code and Kinematics for every case.
The active modes are terminal anchor-Glu chi3, plus terminal extra-Asp chi2 only
when the original requested homolog role is actually ASP. Non-ASP homologs and
every other coordinate stay fixed. Do not select modes from labels or scores.

Verify paired coordinates and state, actual retained role identity, complete
terminal groups, source/cap Jacobian checks and q=0 replay against both original
and context XYZ. Report bitwise identity and numerical differences separately;
use the existing mapping tolerance of 1e-10 Angstrom without modifying inputs.
The existing mapping checks use 1e-7 Angstrom/unit derivative and 1e-9 Angstrom
bond tolerances. Record unsupported cases explicitly in the full 28 denominator.

Emit the existing torsion `design.json` case schema (source, preparation,
role-to-mode IDs, per-metal mapping and origin pins) for compatibility with
second_shell's separate nonlinear runner. Preserve source model/software/ORCA
pins and all first-plan settings by reference. This file does not authorize
reference energies; root owns the integration review and experiment decision.
