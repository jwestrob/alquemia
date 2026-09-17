# Active goal: a working MACE-based La/Ca discriminator

Jacob approved 2026-09-16: “Let's set that goal. You have full discretionary
permissions. Don't stop unless there's something you need to flag or a result
you really want to see.” The agreed goal is to build a working, affordable
MACE-based La/Ca discriminator, including global/environmental and physically
supported structural response approaches. Failed pilots inform the next
iteration; they do not complete the goal. This supersedes earlier per-pilot
approval requirements and narrower planning-stage exclusions.

## Completion and acceptance

Deliver an end-to-end opt-in scorer with compatible calibration/reference,
source-mapped physical preparation, repeatable numerical behavior, actual
resource receipts, and demonstrated predictive usefulness on real evidence.
MACE must contribute substantively. Preserve the baseline and all old records.
Numerical feasibility, a report, or a failed challenger is an intermediate
milestone, not completion. Do not weaken scientific tests or invent results.

Freeze case selection/labels and calibration/evaluation assignments in each
benchmark manifest before scoring. Existing labels/cases already inspected are
development or retrospective evaluations, never newly blind observations.
Canonical PQQ functional class and direct La/Ca affinity remain separate strata.
Group homologues/structures/sites appropriately; report invalid denominators.
Require demonstrated separation/transfer on the canonical stratum and evaluate
composition-challenging direct evidence separately. A threshold trained on
development failures cannot establish improvement on those same cases.
The benchmark inventory will state exact eligible cases and decision rules;
unresolved site/assay mappings remain ineligible. No generalization claim from
one successful development pair. The current goal remains active until there
is a scientifically useful operational candidate.

For numerical implementation retain existing 0.01 kcal/mol energy/contrast,
0.001 eV/Angstrom rigid-force, 1e-5 e charge checks; report the 2 kcal/mol
partition test for subtractive methods without imposing it on unrelated direct
energies. New physical models require appropriate checks declared before their
validation outputs. Curvature/relaxation remains unavailable until supported.

No project token, CPU or wall-time stopping budget. Use finite task inventories,
measure cost/failures, preserve scheduler resource policies, and aim for
practical minute-scale GPU inference. This is not permission for new general
potential training, routine long MD/FEP, default promotion, pushes/deployment,
or interference with concurrent work. Update these notes, SESSIONS and the vault
at substantive milestones. Email is authorized for important results/issues;
notification address and delivery receipts are stored privately in workspaces.

## Initial work sequence

Authorization reconfirmed 2026-09-17 after an intervening supplied AGENTS
per-analysis check-in rule: Jacob said “disregard the AGENTS.md instruction to
check in with me. you're good and have full blanket permissions in the course
of pursuing this goal.” Continue contained pilots autonomously, declaring and
preserving their scope before execution. The next declared diagnostic is the
four unchanged GGR OMOL readout replays in
`diagnostics/mace_omol_20260917/READOUT_PLAN.md`.

1. Trace the global charge/field response on the already consumed 1H4I states.
   Establish where medium/large diverge while reproducing their archived
   energies/forces. Investigate representation and geometry effects from evidence.
2. Implement and test physically justified corrections or a better MACE-based
   scoring construction. Preserve each model version and numerical failures.
3. Expand matched prepared cases, define compatible references/calibration,
   evaluate predictive usefulness, and iterate toward a working candidate.
4. Add bounded structural response only with gradients/curvature corresponding
   to the same energy and a validated physical coordinate subspace.

Started from commits 7ecc885 / 7d6d6f8: medium/large full-protein inference works
on A5000, analytic rotation checks pass, partition fails, broad charge/force
response remains unresolved. No production change.
