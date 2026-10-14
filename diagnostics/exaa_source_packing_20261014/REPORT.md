# ExaA: crystal scores Ca; added packing does not explain the large fold response

The two authorized comparisons are complete. **1FLG chain A gives the expected
Ca call under the unchanged accommodated frozen-charge calibration.** The second
crystal chain remains unavailable at an existing preparation cutoff. The result
supports source sensitivity of the current model, not a conclusion that AF3 is
wrong. The separate finite packing test weakens the proposed missing-contact
explanation for the ambiguous two-Ca ExaA fold.

## What the crystal comparison establishes

Chain A retains the exact Q9 local chemical-state signature, mapped deposited
heavy coordinates and current dry, conditional oxidized-PQQ3− preparation. Its
six common-pool endpoint cells pass all numerical admission gates. Both bounded
MACE searches converge with valid covalent geometry and no active-bound flag.
The solvent-containing selection chooses the **La-generated geometry for both
metals**; the shared-pool calculation therefore matters to the selected state.

| Source/context | Ca selected work | La selected work | Ca−La accommodation contribution |
|---|---:|---:|---:|
| 1FLG A, current context | −5.169804 | −7.099801 | +1.929997 |
| Two-Ca sample0 control, current context | −7.624675 | −12.934776 | +5.310100 |
| Two-Ca sample3 ambiguous, current context | −9.698819 | −26.390560 | +16.691741 |
| Same sample3 saved motions, expanded context | −7.247497 | −23.564450 | +16.316953 |

All entries are kcal/mol relative to each metal's own origin within that same
context. These are accommodation contributions, not affinities or complete
discrimination scores. The enlarged context receives no old reference bands.

Crystal A's raw contrast changes from −405468.1574529917 to
−405466.2274557295 kcal/mol: a small **La-like** change on this score scale.
Its separately calibrated static decision is inconclusive and its accommodated
decision is Ca. These methods have distinct previously frozen canonical bands;
the call transition is **not** evidence that motion shifted its energy toward Ca.
Neither calibration was changed for this experiment.

The full comparison retains all ten current predicted ExaA sources and both
crystal chains. Among the predictions, six accommodated calls are Ca and four
are inconclusive; none are La. Five original single-Ca preparation failures are
shown separately as historical recovery provenance. Folds and chains are
correlated observations of one protein, not twelve independent biological tests.
See the executable [comparison](comparison/PLAN.md) and actual result at
`workspaces/exaa_source_packing_20261014/crystal_saved_comparison_v1/execution/RESULT.json`.

## What the packing comparison rules against

Adding the structurally selected complete surrounding fragments changes the
ambiguous sample's differential accommodation contribution by only
−0.374789 kcal/mol. Both metals experience similar added penalties, and neither
selected geometry changes. This does not support omitted packing as the dominant
explanation for this saved motion's large La response. It does not test a fresh
expanded-context optimization or all possible scaffold responses. Full components
and the matched control are in [the packing report](packing/FINAL_REPORT.md).

## Limits that affect the interpretation

- Chain B's Asp-partner distance is3.511Å versus the unchanged3.5Å admission
  cutoff. This is a coverage limitation, not proof of a defective crystal.
- The crystal's experimental PQQ redox state remains unqualified. Its resolved
  mature sequence maps exactly, but full constructs differ. The dry model omits
  the deposited waters and global assembly response. A crystal input does not
  make these chemical approximations experimentally verified.
- The sole reused terminal OXT addition is remote from the scored region; all
  deposited heavy coordinates remain unchanged. Original preparation failures
  are preserved. See [the preparation report](crystal/REPORT.md).
- This is consumed development evidence. It does not establish that the large
  predicted-fold relaxation is physically wrong, nor independently validate the
  learned metal-response surface.

## Decision and next useful test

Retain the current classifier/calibration and close this narrow packing-dominance
hypothesis. The most informative next question is whether an independent
electronic-structure calculation supports the large ExaA donor-motion response
at the **same saved coordinates and state**. Select a small target/control
comparison and its Hamiltonian before execution; earlier extra-Asp checks do not
validate these Glu/coupled motions. A reference disagreement would motivate a
response-model change; agreement would shift attention toward source state and
preparation. A blanket refold or another arbitrary context expansion is not
justified by these results. No new quantum campaign is submitted in this round.

## Actual work and recovery

All owned jobs are terminal. Packing4091 used681 scheduler seconds on7CPUs
(14 MACE calls including2 controls,48 SCC calls,24 ALPB evaluations). Crystal
native4095 used428 seconds on7CPUs and33 actual MACE forwards. Solvent4096 reused
those MACE results, using128 seconds on7CPUs,24 SCC calls and12 ALPB evaluations.
These counts include independent component checks. Preparation/tests/collectors
have separate receipts; see `FINAL_SCHEDULER_RECEIPTS.txt` in the workspace.
No new DFT, fold, GPU work, threshold fitting or production promotion occurred.

All calculations used the qualified published-parameter standalone GFN2 runtime
and native frozen-charge ALPB replay; the erroneous ORCA default La assignment
was not used. Code, exact configurations, input pins and failures are preserved.
Delayed completion events require acknowledgement only, not a rerun.
