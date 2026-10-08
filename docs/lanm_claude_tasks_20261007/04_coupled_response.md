# 04_coupled_response: Qualify the response needed by the chosen experiment

Read COMMON_RULES.md in full and pass it to every subagent.

## Start condition

mapping/ledger review may parallel01/02; new molecular tests require03 decision and an explicit finite matrix

## Inputs

Reuse metal_environment_force_assembly.py, physical cap/source mappings, force_assembly_20260928 reports, coupled_scaffold_20260928/QUARTER_STEP_REPORT.md and Glu91 native gradients. Existing successful donor and classical N83 checks are not missing work.

## Work

Define the complete energy to be moved and its analytic total forces. Address one actual environment/boundary direction and any omitted solvent/scaffold term needed for the selected experiment. Do not infer full hybrid validity from one donor derivative. Distinguish implementation consistency from metal cross-LJ and solvent model accuracy. Reuse existing passes unless a changed component invalidates them. Predeclare tolerances and displacements from numerical accuracy and physical scale, not desired labels. Prepare new code in owned modules; shared executor changes require integrator review.

## Completion criteria

Deliver interaction ledger, mapped coordinates, exact finite Slurm test manifest, real test residuals when run, supported/unsupported capability table and a go/no-go recommendation for the chosen movement. No invented springs, unrestricted Hessians, eigenvalue clipping, full-protein vacuum subtraction or claimed entropy from an unvalidated matrix.

## Ownership

Own only this task’s worktree branch and diagnostics/claude_lanm_20261007/04_coupled_response/; raw products under workspaces/claude_lanm_20261007/04_coupled_response/. Request a scoped follow-up before editing shared executors/scoring. Root owns integration; this Claude session owns and may submit its assigned scientific jobs.
