# Mex fixed eight-cell embedded response adapter

Source-only implementation of the frozen inventory and proposal. No job or new
reference result is implied. `INVENTORY.json`, `SPEC.md`, and
`ADAPTER_PROPOSAL.md` remain unchanged. This executable declaration requires the
**full** original eight-point derivative qualification before any response call;
the earlier partial scout pair alone does not qualify this matrix.

## Fixed inputs and ledger

Use exactly the eight inventory cells, in the original ordering:

| Geometry | La cell | Dy cell |
| --- | --- | --- |
| Old own-La restricted endpoint | `La_terminal__La` | `La_terminal__Dy` |
| Old own-Dy restricted endpoint | `Dy_terminal__La` | `Dy_terminal__Dy` |
| Broad own-La terminal | `La_cartesian_terminal__La` | `La_cartesian_terminal__Dy` |
| Broad own-Dy terminal | `Dy_cartesian_terminal__La` | `Dy_cartesian_terminal__Dy` |

The old endpoints are the actual broad-search starts. Both broad searches ended
at 80 accepted steps without convergence. All 2,049 physical coordinates, 167
represented waters and the three Nd spectators remain fixed to these archived
states. No new geometry, displacement, optimization, solvent, MM, GFN or MACE
evaluation is part of this adapter. Eight new DFT cells is a ceiling, not a claim
of execution. The full denominator survives preparation, partial collection and
failures.

The reference is `PBE0_embedded + same_MM_ledger + same_fixed_FF_ALPB`.
It replaces both MACE vacuum and the GFN embedding difference. The four existing
MM groups are counted once, with their native exclusions; fixed-charge ALPB uses
the full unshifted physical charge vector. The shifted 1,840-charge electronic
field remains distinct. No extra Coulomb, induction or continuum term is added.
This reference shares MM stiffness with the approximate model and cannot
independently validate that stiffness.

## Existing implementation reuse

`QUALIFY.py::inputs` replays the original native PBE0 source/state/runtime and
assets. The actual common7 pool is replayed using its saved baseline, search
admission, request import, native component and selection functions. Its scoring
entry point is never called. All component-call counters remain zero and no
MACE model initializes.

Current POINT, trial, source, full physical array, core/field/boundary and cap
identities are checked. The existing mapper validates cap coordinates and
Jacobians under its unchanged tolerances; archived current cap Jacobians remain
the derivative authority. MM receipts, native group energy closure, source-ID
ordering and saved component gradients are replayed. ALPB raw stdout, fixed
charges, full coordinates, binary/command, native gradient and energy closure are
replayed. Both metals share exactly the same physical coordinates at each alias.

Each prepared task copies the original same-metal Mex `EnGrad` deck byte for
byte, including `DoEQ false`, the custom basis/ECP/auxiliary data and PModel.
Core and field files are copied from the current inventory cell. New tasks carry
current physical/map provenance plus an explicit original-task reference. The
effective singlet has 832 explicit electrons; La/Dy physical multiplicities 1/6
and ECP cores 46/55 remain separate. No SOC or stable-state claim is introduced.

## Stages and execution gates

All stages require Slurm and the existing pinned Python executable. Root owns
submission. `CONFIG.template.json` contains actual existing pins and null future
evidence. Null values do not stand in for passed gates.

1. **prepare:** zero molecular calls. Replay all saved inputs, freeze both full
   response chords, and write eight tasks, eight shared-component receipts and
   one full manifest. A successful preparation is not execution permission.
2. **preflight:** zero molecular calls. Require the actual preparation and replay
   current physical/shared components again. Report missing future gates. If a
   full derivative result is supplied, replay it through the frozen collector.
3. **endpoint:** one named native analytic energy/gradient task only, with 12
   actual MPI ranks, one worker and allocation-derived memory. Require actual
   zero-call preflight, saved replay of all eight derivative tasks/four checks,
   and a pinned root measured-cost decision. The first task is always
   `La_terminal__La`. Any of the remaining seven additionally requires that
   actual admitted response scout and root's separate measured-cost decision.
   No loop, automatic retry or changed SCF protocol exists.
4. **collect:** saved data only. Preserve all eight slots, replay available raw
   endpoints and full physical gradients, and report unavailable cells and actual
   per-allocation costs. The complete four-geometry paired pool requires all
   eight cells; partial source-path comparisons retain their missing members.

Commands, with root-supplied actual config paths and fresh output directories:

```text
PYTHON RESPONSE_REFERENCE.py --stage prepare --config CONFIG.json --output PREPARATION
PYTHON RESPONSE_REFERENCE.py --stage preflight --config BOUND.json --output PREFLIGHT
PYTHON RESPONSE_REFERENCE.py --stage endpoint --task La_terminal__La --config BOUND.json --output SCOUT
PYTHON RESPONSE_REFERENCE.py --stage collect --config ACTUAL_RESULTS.json --output COLLECTION
```

The root cost-decision records are deliberately small JSON evidence records.
Both require `protocol`, `status: approved`, `scope`, `allowed_tasks`, exact
`qualification_result`, `preparation_result`, `preflight_result`,
`qualification_molecular_wall_seconds` for every derivative task, nonempty
`measured_cost_review`, and pinned `cost_evidence` files. The first scope is
`first_cost_scout`, with only `La_terminal__La` allowed. The second scope is
`remaining_seven`, with the seven other frozen IDs, plus exact
`response_scout_result`, `first_scout_decision`, and
`response_scout_molecular_wall_seconds`. Root supplies these only after reviewing
actual measurements. There is no guessed time estimate or synthetic approval.

## Native gradient admission

The existing allocation-aware manifest runner and native energy/SCF admission
are unchanged. `QUALIFY.admit_sp` is reused for state, receipt, memory, normal
termination, PBE0/D4/field and convergence semantics. Analytic task type and a
nonempty engrad path are mandatory. The execution receipt must certify the
expected engrad and match its exact task/manifest/artifact.

The existing generic frozen-f `parse` admits actual SCF/ECP/D4 analytic-gradient
evidence and QM order/coordinates/energy. Its hardcoded q−1 state/input helpers
are not used. Native QM coordinates retain the existing 1e−6 Å check and native
gradient energy the existing 1e−8 Eh check. The separate raw pcgrad must contain
exactly 1,840 finite ordered derivatives and is pinned explicitly, since the
runner does not include pcgrad in its standard artifact list.

Both native derivatives are +dE/dR in Hartree/bohr. The existing cap scatter maps
real QM atoms, both cap anchors and point charges exactly once to physical
kcal/mol/Å. Add saved MM and ALPB gradients once; force is minus that total.
Preserve raw gradients, decimal print precision, SCF observations and residual
flags. No tighter SCF interpretation, density override or force-stability claim
is added.

## Frozen reporting targets

These are discrepancy-reporting targets, **not native admission gates**:

- **1 kcal/mol** for each metal's deformation work and for the differential.
- **1 kcal/mol/Å** for the difference between reference and approximate gradient
  projections along each of the two normalized full physical response chords.

The chords are old own-La → broad own-La and old own-Dy → broad own-Dy,
normalized by the Cartesian L2 length across all 2,049 atoms. Both are reported
at every available endpoint, without outcome-based atom or direction selection.
No alignment, source window or donor-only projection is introduced. Every atom's
reference force, approximate force and difference is retained, together with
whole-system maximum/RMS force and projection components. Small projections
cannot replace full-atom inspection.

Report all four old-to-broad paths for each metal, the own-metal-path
differential and source-local shared four-geometry menus for both methods.
Use native per-metal energy differences before conversion; PBE0 subtraction uses
the exact printed energy decimal strings. Show unchanged MM/ALPB work separately
from the replaced electronic work. Never compare absolute totals across methods.
Menus report mathematical minima, exact ties and a descriptive 1 kcal/mol band;
they do not produce affinity/classifier labels.

Numerical derivative qualification remains the frozen
`0.005 + 0.0001 × max(|analytic|, |finite difference|)` kcal/mol/Å rule at its
original fixtures. Passing those directional checks is neither a global force
qualification nor an error bound on the broad-motion chords. SCF uncertainty,
finite print precision, frozen-f limitations, conditional EF3 occupancy and
nonconverged search endpoints remain explicit. PBE0 is an independent electronic
comparison under this declared model, not an oracle.

## Verification before launch

Source-only AST parsing and independent static review precede root's numerical
preflight. The actual Slurm preparation/preflight must establish saved common7,
current-map, charge/state, MM/ALPB and native input replay. Missing full derivative
qualification, preparation/preflight or either applicable measured-cost decision
must fail before any requested native call. Actual native endpoint admission and
saved collection are separate later stages. No numerical work has been run by
the source author.
