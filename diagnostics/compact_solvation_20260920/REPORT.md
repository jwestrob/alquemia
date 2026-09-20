# Matched solvent transfer improves compact MACE discrimination

**Completed 2026-09-20. Pursue this challenger; preserve the released default.**
The fixed solvent correction improves core-MACE alpha/GGR ordering from2/6 to6/6,
retains all25 canonical PQQ cases and all3 consumed crystal transfers, and widens
the complete-context PQQ gap from2.315 to5.076 model-kcal/mol. No correction
coefficient, geometry, chemistry or label was fitted to obtain this result.

This is a preliminary improvement to MACE on already-used cases. It does not
establish greater accuracy than the promoted water-prepared DFT baseline, which
already gets these alpha/GGR directions right. Six structural comparisons are
**one biological direction**, not six independent biological successes. Calibration
accuracy and the three retrospective transfers do not establish broad accuracy.

## Actual result

| Representation and model | Canonical PQQ | Crystal transfers | PQQ class gap | Alpha/GGR ordering | Weakest alpha/GGR margin |
|---|---:|---:|---:|---:|---:|
| Core native vacuum MACE | 25/25 | 3/3 | 79.031 | 2/6 | -4.406 |
| Core MACE + solvent transfer | 25/25 | 3/3 | 29.949 | **6/6** | **3.209** |
| Complete-context native vacuum MACE | 25/25 | 3/3 | 2.315 | 6/6 | 13.371 |
| Complete-context MACE + solvent transfer | 25/25 | 3/3 | **5.076** | **6/6** | **14.258** |

Gaps/margins are model-kcal/mol. Larger positive class gaps and pair margins
indicate the desired ordering; they are not accuracy probabilities. Every model
in this table orders all154 cross-class PQQ calibration pairs correctly. The
core solvent correction reduces the original large PQQ gap while fixing the
four failed alpha/GGR structural pairs. Context retains stronger alpha/GGR
margins and its initially narrow PQQ gap more than doubles. This is a useful
tradeoff, not a claim that every individual margin improved.

Crystal calls remain 1H4I=Ca,4MAE=La,1KB0=Ca, with both representations. New bands
use the frozen class-extrema procedure on the25 calibration cases only. Transfers
are excluded. Exact raw scores, class spreads, groups, bands and all pair margins
are in [RESULT.json](RESULT.json); full receipts remain pinned under workspaces/.

### Does the MACE component contribute?

The observed combination cannot be replaced with standalone GFN2 on this set:

| Representation | Electronic model | Ordered PQQ cross-class pairs | Alpha/GGR directions |
|---|---|---:|---:|
| Core | GFN2 vacuum | 150/154 | 4/6 |
| Core | GFN2 ALPB water | 135/154 | 4/6 |
| Context | GFN2 vacuum | 148/154 | 2/6 |
| Context | GFN2 ALPB water | 146/154 | 2/6 |

Their class gaps are negative (overlap), so no separating band is reported.
MACE plus the GFN2 **solvent difference** preserves the strong MACE PQQ ordering
while improving the difficult core comparison. These are descriptive component
checks on the same consumed set, not another fitted classifier or independent test.

## What changed

The archived132 native OMOL endpoints cover33 cases in both core and context
representations. Coordinates, endpoint-specific water-H preparations, charges,
singlet states, water inventories and model checkpoint were reused exactly.
No new geometry, DFT or MACE calculation was run. The new low-level calculations
evaluate vacuum and ALPB-water GFN2 on each matching endpoint:

```
E_composite,M = E_MACE,vac,M + E_GFN2,ALPB,M - E_GFN2,vac,M
R_composite = E_composite,Ca - E_composite,La
            = R_MACE + delta_solv,Ca - delta_solv,La
```

MACE eV and GFN2 hartree are each converted once to their recorded common scale.
The GFN2 vacuum molecular energy is subtracted, so it is not added to MACE twice.
The correction includes GFN2's change of density/solvation energy between media;
it is not a reaction-field energy evaluated on the MACE density. MACE produces
no endpoint charge distribution here. This is an approximate composite descriptor,
not self-consistent MACE/continuum physics or a binding free energy. ALPB uses the
supplied cluster boundary, not the shape of the entire protein.

ORCA6.1.1 Native-GFN2-xTB, installed default parameters, ALPB(Water),300K electronic
smearing and the native xTB mixer are fixed. Actual native Ca/La parameter exports,
electron counts, charge/spin and atom order were checked; no DFT ECP convention
was silently imposed. The [verified native-method manual](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb)
documents this implementation and its mixer/solvent controls.

Protocol: `native_OMOL_plus_native_GFN2_ALPB_transfer_v1`, with explicit core/context
representation keys. There is **no compatible aquo reference**; S remains null.
The very large raw MACE atomic-energy offset carries no affinity interpretation.
New calibration bands are separate from the released DFT score and decision bands.
No production/default/reference/watcher was changed.

## Numerical and physical checks

- All264 primary endpoints converged and passed actual state, energy and charge
  accounting. All66 case/representation scores are available; none substituted a
  baseline result. Maximum absolute atomic Mulliken charge was1.044852e; all
  charge sums met the declared5e-4e tolerance.
- The eight ordinary-SCF/TightSCF numerical checks all failed to converge after
 124 cycles. They have no usable energies. Status is
  `numerical_crosscheck_unavailable`: the predeclared0.10kcal/mol transfer and
 0.20kcal/mol score-agreement tolerances are **untested**, not passed or relaxed.
- Identity/unit/sign algebra uses actual archived GFN2 values. It checks exact
  zero transfer for identical inputs and the Ca-minus-La sign. This is algebraic
  validation, not a new solvent-boundary convergence experiment.
- Native output provides totals and SCF components, but no separate ALPB energy,
  cavity-radius export or electronic entropy component. Those entries stay null.
  Runtime inputs and installed executable fix the actual solvent selection.
- Nonfatal PMIX diagnostics appeared despite normal, converged calculations and
  valid receipts. Shared permissions/MPI infrastructure were not modified.
- No whole-protein partition, boundary, grid or rigid-transformation validation
  is claimed for this cluster-ALPB experiment. It does not rescue the older
  failed global environmental implementations.

**Judgments:** Predictively useful on the tested development set: yes.
Computationally practical at this scale: yes. Independently numerically qualified
and broadly validated: not yet. Internal convergence alone does not settle the
composite model's systematic errors.

## Measured cost

| Job | Actual work | Wall seconds | Allocated CPU-seconds | Result |
|---|---|---:|---:|---|
| 1203145 | Eight primary pilot endpoints | 24 | 1,536 | 8/8 complete |
| 1203147 | Eight alternate-SCF checks | 60 | 3,840 | 8/8 failed convergence |
| 1203165 | Remaining256 primary endpoints + collection | 500 | 32,000 | 256/256 complete |
| **Total new work** | **272 low-level attempts** | **584 summed job-seconds** | **37,376** | **Zero GPU-seconds** |

Each job allocated64CPUs/128GiB and used eight concurrent eight-rank tasks. Full
primary compute took454.873s; collection/startup account for the remainder of the
500s allocation. This is development-panel cost, including both representations
and the failed numerical checks, not the cost of one production site.

The median sum of four eight-rank low-level task durations per site was43.665s
for core and70.177s for context (about349 and561 task-rank-seconds respectively).
Those are serial sums; parallel site latency was not separately measured.
Recorded archived MACE endpoint-pair wall times had medians4.175s(core) and2.260s
(context), on mixed A5000/H200 versus H200 hardware. These are reused inference
timings, not new GPU allocations or a controlled comparison of representations.
Preparation/folding and allocation overhead prevent claiming a matched end-to-end
speedup from these numbers. ALPB adds a modest CPU stage; it is not free inference.

[COSTS.json](COSTS.json) retains scheduler records including failed attempts.
[PRIMARY_ENDPOINT_COSTS.json](PRIMARY_ENDPOINT_COSTS.json) and
[REUSED_MACE_COSTS.json](REUSED_MACE_COSTS.json) retain actual task timings.
Observed batch MaxRSS was4,328,500KiB(full) and4,609,668KiB(numerical); pilot0 is
unavailable accounting, not zero memory usage. These are not guaranteed aggregate
MPI-process memory peaks. Requested128GiB is recorded separately.

## Delivered and next step

Focused preparation/validation/collection code, source comparison code, fixed
manifests,264 real primary energies,8 visible numerical failures, all66 paired
scores, comparison/calibration tables and executable commands are saved. Fourteen
real-fixture parser, source, cache/state and algebra tests pass; no mock scientific
successes or skipped executable tests. Scientific integration is the actual272
attempts above, distinct from those14 software tests.

**Recommendation: pursue this challenger.** It has earned further work through
observed utility. Resolve numerical reproducibility with a contained same-state
check before promoting it, and measure the final scanner path on matched hardware.
Further generalization needs additional trustworthy biological directions; simply
rescoring these consumed cases cannot supply that evidence. Do not retune the
coefficient, state or threshold to these outcomes.

Read/replay the saved experiment using [COMMANDS.md](COMMANDS.md). The next
read-only command from the repository is:

```bash
cat diagnostics/compact_solvation_20260920/REPORT.md
```

The independent CC job1202429 remains separate and is handled by its existing
completion monitor. It was neither restarted nor expanded for this result.
