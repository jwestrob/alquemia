# Three-source executor is reviewable; predictive qualification did not pass

The new thin adapter removes the historical-canonical/A0A3 gate from the existing
warm execution route. It accepts explicitly prepared three-La-source groups,
retains source/protein identity and the frozen threefold reference, and uses the
same four-mode selector, optimizer, common candidate pool and native scalar
runner. No molecular calculation was launched.

The intended two PLM triples pass preflight: six sources, twelve paired native
origin tasks and twenty-four original-coordinate solvent inputs. The fixed
upper scope is twelve searches, twelve cross-MACE calls and seventy-two native
GFN2 cells, with no DFT or protonation. Their biological labels stay unknown.
Exact prior source preparation is reused; no old energy or force satisfies this
prospective run. Protein IDs remain compatible with the existing 176-protein
export; no biological joins, transcript analysis or historical DFT fields changed.

The broader threefold transfer reported by root has 92 correct/2 inconclusive/
6 unavailable versus 94 correct/6 unavailable for the released/tenfold comparison.
This representation has therefore not retained full reference fidelity. The two
new A8R3S4 inconclusives and all missing cases stay visible. This delivery is an
experimental capability, **not a qualified scorer or promotion recommendation**.
The molecular PLM usability run remains held; see [status](STATUS.md).

## Actual verification

Seven real-artifact tests pass. They cover explicit unlabeled requests, exact
paired origins and one-rank/MaxIter500 input recipes, unchanged shared defaults,
the real MMOL missing-member failure, actual incomplete collection, archived
energy/pool/median algebra and corrupted-copy rejection. No fabricated molecular
successes. Immutable snapshot dry-run, collection and report operations were
also exercised. The incomplete result correctly retains all six sources and
both group medians as unavailable, without installing a stale authoritative
stage collection that could interfere with later execution.

The real unsupported MMOL fixture emits zero origin tasks for its full three-
source group, preserving the original sample-4 overlap defect. A failed required
matrix cell or missing source makes the accommodated result/strict group median
unavailable; the conditional unaccommodated union contrast is retained separately.
It is not relabeled as a released local-context score or substituted as success.

The output includes both mathematical and operational R, all component energies,
selected identities, endpoint work, candidate boundary/termination information,
source evidence, frozen protocol/reference pins, ordered source scores and strict
median-three/range. No probability, affinity free energy or biological accuracy
is implied by these unlabeled inputs.

## Implementation and limits

`scripts/pqq_three_source_execution.py` configures private copies of the already
used union candidate engine. The shared executors, production defaults and frozen
preparation API remain unchanged. The new adapter owns source handoff, scope
checks and group reporting; it does not implement another optimizer or solvent
workflow. Existing physical/state/cache checks run before each stage. Current
resources remain one H200/32 CPUs/200000 MiB and32 one-rank scalar workers.

All actual molecular, protonation and Slurm counts are zero. Code/preflight/tests
used local CPU only; this development overhead was not separately metered. Prior
measured integration provides an expected minutes-scale run estimate in the
handoff, not a result or a time limit. Engine timing would include only a prepared
metadata handoff here; its historical `source_preparation_included` field does
not mean fresh protonation, and the result explicitly records that qualification.

An existing execution-start marker or receipt makes a repeated execute call stop
for collection; automatic scientific retry/rescue is unsupported. The wrapper
collects/report errors while retaining the original execution exit status.
No scientific launch commands have been run. [Commands](COMMANDS.md),
[plan](PLAN.md), and [actual pins](ARTIFACTS.json) make this checkpoint reviewable.
