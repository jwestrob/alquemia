# September22 delivery: retain the supported fast scorer, close unhelpful searches

**The released MACE+GFN2 context scorer remains the strongest practical method
in the completed structural challenge. The new accommodation variants have not
earned routine use.** This is a completed implementation/comparison delivery,
not completion of the broader accuracy-improvement goal.

## What the recovered work establishes

All earlier preserved DFT and fold-proposal work is finished. On the207 sources
covered by both methods, released context scoring gives203correct/2wrong/
2inconclusive versus DFT198/3/6. Both retain all94 available three-fold summaries
correctly. These are repeated, consumed structures from25 protein groups, not
207 independent biochemical labels. [Recovered evidence](REPORT.md).

The previously completed native DFT displacements also confirm a real mechanism:
MACE captures metal-dependent donor-strain relief in two PLM structures. That
supports structural sensitivity; it does not establish those proteins' metal
preferences or validate the newly optimized candidates.

## New work and decisions

| Component | Actual result | Decision |
|---|---|---|
| Shared geometry pool | Original30 and full225 executed; no added discriminatory benefit over prior proposals; separately calibrated transfer weaker than released static | Working reusable primitive; omit routine extra cells |
| Force-selected four-angle search |29/30 sources score; on24 matched canonical members static24correct versus adaptive22correct/1wrong/1inconclusive | No promotion; missing member leaves new calibration unavailable |
| Joint metal/four-angle search |7/8 candidates,3/4 paired pools; native metal forces relieved, available decisions unchanged | No broader campaign from this pilot |
| Scaffold compliance | Real protein parents parameterize; capped-local subtraction/exterior cofactor coupling unresolved | Score unavailable; no invented stiffness or correction |
| Nikasha naming | Thin entrypoint and current docs added | Legacy invocations/results/defaults retained |
| PLM integration | Joinable176-protein/200-gene export, actual expression measures, methods outline and editable figures | Ready for manuscript assembly; existing PLM DFT scores remain identified as DFT |

Details: [full pool](../nikasha_shared_pool_20260922/FULL225_RESULT.md),
[angular comparison](../adaptive_accommodation_20260922/CANONICAL_POOL_COMPARISON.md),
[joint response](../adaptive_metal_20260922/FINAL_RESULT.md),
[scaffold feasibility](../scaffold_environment_20260922/REPORT.md).
Both optimizer failures are invalid trial steps, not invalid biological inputs.
Changing the optimizer would create another declared development version, not
retroactively turn these failures into successes. No proton inventory, water
state, cofactor state, geometry label or score threshold was changed to rescue
an inconvenient prediction. No new chemical ambiguity required titration work.

All launched molecular jobs are terminal; no duplicate collectors remain. The
new shared-pool/adaptive/joint work used263,224 allocated core-seconds and607
GPU-seconds in total, including the explicit numerical diagnostic and failures.
This excludes reused historical calculations and local preparation/reporting;
it is development allocation, not scanner latency. No new DFT was launched.
Final checks:27 pool/comparison/joint tests and9 PLM export tests pass; the released
standard plan validates. [Verification and runner details](VERIFICATION.md).
The preserved separate candidate-DFT campaign remains dry-run only. The native
GFN2 iteration-ceiling repair reproduces controls and keeps its original failed
attempt visible; it does not change the Hamiltonian or convergence tolerances.

## Practical delivery

Use [Nikasha's current commands](../../docs/NIKASHA.md), with the released fast
route for supported explicit PQQ sources and `--mode dft-reference` for preserved
DFT access. An actual completed plan can be checked without launching work:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/nikasha standard validate \
  --plan workspaces/pqq_fast_release_20260920/standard_1H4I_v1/plan.json
```

For a new real input request, preparation/execution/collection commands and an
actual source example are in the [release guide](../pqq_fast_release_20260920/COMMANDS.md).
The optional three-source aggregate remains available. Research collections
retain original/accommodated contrasts, selected geometries, per-metal work,
components and explicit unresolved reasons. They are not production substitutes.

Manuscript material: [methods outline](../nikasha_manuscript_20260922/METHODS_OUTLINE.md),
[exact technical supplement](../nikasha_manuscript_20260922/TECHNICAL_SUPPLEMENT.md),
[structural-comparison figure](../nikasha_manuscript_20260922/STRUCTURAL_COMPARISON_CAPTION.md),
[actual DFT/MACE response figure](../nikasha_manuscript_20260922/FIXED_DISPLACEMENT_CAPTION.md),
[joinable PLM tables](../nikasha_plm_export_20260922/REPORT.md).
No full PLM cohort rescore, remote push, package publication or production
promotion was performed. Paper-facing AI results must distinguish the MACE
reference/response evidence from the existing DFT-scored environmental cohort.

**Recommendation:** retain the released fast protocol for the PLM work. The next
research question should challenge the scoring physics on existing known-class
structural failures; simply adding more freedom to these successful local
searches has not demonstrated a classifier benefit. Numerical search safety and
scaffold subtraction remain specific engineering problems, not a reason to
delay use of the supported scorer or to claim the broader goal achieved.
