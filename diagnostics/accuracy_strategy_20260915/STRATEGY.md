# Accuracy-first development direction

2026-09-15. **Strategy review and evidence curation; no new calculation plan approved or executed.**

## User direction and division of work

Jacob: “We need to think about how to improve performance rather than quantify
existing performance because we can accomplish that in parallel.” He explicitly
requested a subagent for evidence curation, then added: “Make sure the subagent
makes a vault note when it's done.”

The evidence agent owns the challenge-panel curation and its vault note. The
main agent owns the model-improvement direction. Two other agents reviewed
existing environmental and mechanical work without running calculations.
This division does not authorize a new model, preparation, fold, or score run.
The earlier removal of compute/time stopping budgets remains in effect.

## What should change in our approach

The next experiment should implement a candidate improvement and directly ask
whether it improves discrimination across the selected development cases.
Physical checks are necessary support for that experiment; producing another
set of validated derivatives without a usable correction is not its objective.

The previously proposed four additional connected-GGR endpoints on 2FW0/2FVY
are **shelved**. They would extend a representation-sensitivity observation,
but would not supply a general improved scorer. This prospective decision
supersedes the next-step suggestion in the completed GGR report; it does not
change that report's results.

## Preferred accuracy hypothesis

**Recover the intact protein's interactions with the binding site while
keeping the expensive quantum region small.** The desired physical model
contains the whole protein; the whole protein need not receive a dense quantum
calculation. Permanent electrostatics, the protein/solvent boundary, and the
local electronic distribution must belong to one coherent energy model.

Reasons to prioritize this over another geometry-only campaign:

- Alpha-lactalbumin has the same qualitative ordering failure from both its
  Ca-derived and La-derived crystal structures. Selecting the native metal's
  geometry alone has already been tried there.
- Its carbonyl fragments omit the Asp84 side-chain charge and Lys79 side chain.
  These are concrete examples of interactions absent from the representation.
  Their individual effects or net direction have **not** been established.
- Connected versus extended GGR changes the contrast by about 7.34 kcal/mol.
  The cavity and chemical context change together, so this supports sensitivity
  to representation, not a causal assignment or proof of improved accuracy.

This is a **ranked hypothesis**, not a demonstrated explanation of either
failure. Hydration, deformation costs, and conformational populations remain
possible missing contributions. A purely electronic descriptor may not capture
all solution-affinity directions even after environmental improvement.

## Engineering conclusions from existing attempts

| Attempt | What it establishes | Consequence |
|---|---|---|
| Point-charge embedding | Moving Asp across the partition caused 63–68 kcal/mol jumps despite charge closure | Do not rerun that construction as the new model |
| Frozen-charge APBS transfer | Refinement left roughly 10.03 kcal/mol partition sensitivity; direct and isolated-reference contributions changed strongly and oppositely | A finer grid alone is not a defensible remedy |
| Native whole-protein GFN2 | No converged energies; dense startup/memory made the tested route infeasible | This rejects that computational route, not the usefulness of protein context |
| Native ORCA CPCM/B | Documented solvent coupling solves the large QM2 system, then supplies its field to small calculations | Whole-protein QM2 retains the global electronic cost; outer-MM participation in the cavity was not established |

No reviewed implementation currently provides a ready, affordable replacement.
Before proposing scientific runs, the model needs an explicit solvent and
boundary energy expression and an identified maintained implementation.
This is a concrete unresolved engineering issue, not permission to launch
another exploratory embedding calculation. A generic ONIOM equation or a
plausible input keyword does not resolve it.

If such an implementation is established, use the consumed alpha/GGR cases
for development and retain PQQ controls to detect loss of useful behavior.
Specify the actual cases, model, parameters, and runs in a subsequent proposal.
Improved relative ordering is the first performance question; a changed model
cannot inherit the released PQQ bands. The agent's challenge-panel curation
identifies potential independent and mutant-response tests for later use.

## Secondary accuracy hypothesis: metal accommodation

The strongest historical force-field negative concerned **La/Dy** denticity
and weak-site retention, not broad La/Ca discrimination. The early static
12-6-4/BVS report also predates the pair-specific C4 audit. Neither supports a
general claim that metal-specific adaptation cannot help La/Ca.

A documented TIP3P Ch-BE La/Ca parameter arm remains a possible geometry
generator; it is not a validated affinity scorer or curvature model. The
historical Dy outer-water study already found a hydration effect but returned
NO_CALL, so “try adding water” is not a new general solution. Actual mechanical
corrections remain unavailable: the new analytic gradients alone do not supply
a defensible stiffness model.

Accommodation is a secondary proposal, not an automatically scheduled run.
No arbitrary spring constants, new water states, or fitted label corrections
are introduced here.

## Sources and durable outputs

- [Completed GGR study](../ggr_mechanism_plan_20260915/REPORT.md)
- [Alpha failure audit](../benchmark_set_20260915/ALACTA_FAILURE_AUDIT.md)
- [Embedding result](../pqq_balanced_embedding_20260914/RESULT.md)
- [APBS completed pilot](../affordable_challenger_20260915/COMPLETED_PILOT.md)
- [Global feasibility archive](../global_representation_20260915/ARCHIVE.md)
- [Canonical 12-6-4 result](../../HANS_LANM_AMBER_1264_CAPABILITY_RESULT_2026-08-04.md)
- [Historical Ch-BE proposal](../../docs/canonical_amber_1264_relaxation_pilot.md)
- [ORCA multiscale solvation manual](https://www.faccts.de/docs/orca/6.1/manual/contents/multiscalesimulations/qmmm-molecules.html#solvation)
- [Challenge-panel curation](CHALLENGE_PANEL_CURATION.md)

Vault curation note:
`/home/jwestrob/jwestrob/obsidian-vault/agent-captures/2026-09-15_laca-accuracy-challenge-panel-curation.md`.

Baseline/default, historical experiments, thresholds, and the original benchmark
ledger are unchanged. No scientific executable was run for this strategy review.
