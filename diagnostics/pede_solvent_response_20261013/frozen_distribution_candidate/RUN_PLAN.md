# Finite candidate execution — 14 October 2026 UTC

This fixes the execution population before new candidate scores. The physical
model and interpretation remain in [PLAN.md](PLAN.md). Production is unchanged.

## Inputs and gates

- Inventory: `workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_inventory_v1/INVENTORY.json`,
  SHA256 `0cc05bebae04cb12c708ad67173b934884080d3c00185a17ac808db3bf0c1af1`.
- Canonical calibration: exact designated 25 members (14 Ca, 11 La), 150 distinct
  source/geometry/metal cells; 300 native component calls, with a water-charge
  reproduction check before each vacuum-charge evaluation.
- Operational selection keeps the origin unless a candidate lowers the energy
  by strictly more than 0.1 kcal/mol. Both metals see the same complete pool.
  Mathematical minima and static origins remain separate reported policies.
- Calibration uses the existing rule: Ca upper bound is the largest canonical
  Ca score; La lower bound is the smallest canonical La score. All 25 members
  and a gap greater than 0.02 kcal/mol are required. Operational calibration
  must pass before transfer. Its exact result is frozen before transfer starts.
- Transfer: 1,392 additional native inputs, of which four exact PedE origin
  replays are reused; 2,776 new component calls. These cover the expanded 135,
  original 125 Ca-source denominator, separately identified recoveries, all
  15 two-Ca folds, and the existing 100 grouped triples. These populations
  overlap; neither structures nor triples are independent biological replicates.
- Native water-component residual bounds and the maximum informative bound are
  unchanged from the qualified replay. Apply the same precision ceiling to
  vacuum-charge components before admitting them. Unavailable cells remain
  unavailable; no old score substitutes for a failed new score.

No new MACE, SCC, DFT, optimization or folding calculations. Native replay uses
the pinned successful xTB 6.7.1 static libraries and a compiled variable-size
wrapper. All compilation, tests, component evaluation and score collection run
on Slurm. Native workers use at most 12 allocated CPUs, one thread per child;
the scorer uses one CPU. No memory/time override or priority boost.

## Regression already executed

Slurm job 4006 passed saved-result regression over all 278 inventory records,
25 canonical members and 100 declared triples. It reproduced archived baseline
selection, contrast algebra, calibration and decisions; missing-input tests use
explicitly corrupted copies of real fixtures. This was one allocation-second
on one CPU, with zero new molecular calls and zero candidate scores.

Exact receipt: `workspaces/pede_solvent_response_20261013/frozen_distribution_score_tests_v1/TESTS.json`.
Tested scorer SHA256: `9589c60af4f977c012eefe55de06eba044a50582a9c994c602da44b53e054375`.

## Interpretation

The first question is reference fidelity. If calibration fails, stop the
dependent transfer and retain the failed result. If it passes, compare fixed-band
transfer on common coverage, every unavailable input, selected geometries and
accommodation works. Preserve the original primary and recovery populations.
Do not change bands after inspecting transfer, fit a response multiplier, or
interpret increased abstention as improved accuracy. A successful training gap
alone is not an improvement. No production promotion follows from this test.

Execution workspaces retain manifests, exact commands, scheduler receipts,
completion-wake commands and result hashes. Root owns submissions and scoring;
the adapter and inventory agents own only their scoped implementation/review.
