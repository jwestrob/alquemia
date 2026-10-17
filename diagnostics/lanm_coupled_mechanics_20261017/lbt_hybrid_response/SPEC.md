# Full physical response under one explicitly selected LBT hybrid region

## Status and question

Source-only, unbound draft. This specification does not authorize either region's
execution. Root decides a finite pilot only after interpreting the actual complete
104-cell fixed-pool comparison and qualified derivatives. No submission occurs
from this wrapper. It cannot start an additional region automatically.

Question: with the hybrid Hamiltonian fixed, how much additional lowering and
source consistency can the existing full-coordinate optimizer obtain from each
source-local paired-pool minimum? This is a mechanics test, not a search for a
desired La/Dy sign or a binding free-energy estimate.

## Inputs and unchanged model

One config declares exactly one region, `baseline` or `expanded`, and all four
ordered states: 7CCO/La, 7CCO/Dy, 7CCN/La, 7CCN/Dy. The 7CCN coordinates remain
Lu-conditioned; relabelling the substituted electronic state does not change
their structural provenance. There are 256 physical atoms including the target,
no represented waters or spectator ions. Both alternatives use this same wrapper
and admission rules; they must be separate, explicitly bound experiments.

The parent must be the actual complete 26-geometry/104-cell two-region hybrid
pool. Frozen pool `prepare()` replays all prerequisite evidence, including the
nine-cell baseline-La scout, eight origins and remaining 24 derivative cells.
The wrapper then restores all 104 actual hybrid POINTs with the existing
qualifier's `restore()` and exactly reproduces the parent selections and works.
Every compatible full-system ALPB cache entry is restored. No molecular calls
are allowed during preflight. Source paths may be copied byte-for-byte; pins
compare verified SHA256 values where path aliases are legitimate.

The evaluator stays EVALUATOR_V2 `686eba84fd3dfcf28724d2d4e7957528c987445d7137b18c88a07dfec055f661`.
Its full corrected La/Dy parameter profile, native state/electron bookkeeping,
172/210-atom local regions, cap/point-charge mapping, four MM groups and fixed
256-atom ALPB charges are unchanged. Energy is local MACE vacuum + embedded-minus-
isolated GFN2 + declared MM ledger + global fixed-charge ALPB. All four component
gradients are transformed to the same 256 physical coordinates by V2. No whole-
peptide MACE term or duplicate Coulomb term is added. Every new point retains
current PHYSICAL/TASK/core/field provenance. V2's narrow, diagnosed cap-Jacobian
replay allowance is unchanged; there is no coordinate, energy or force tolerance
relaxation here.

The existing two-direction derivative qualification is required for both metals
and regions. It is not proof of accuracy for every physical degree of freedom or
of the Dy electronic-state approximation. This wrapper adds no directions,
force qualification, parameter fitting or DFT calls.

## Starts, motion and physical admission

For each metal, the start is its selected minimum in its source-local full parent
pool, under the declared region. Exact ties retain the parent candidate order;
all ties are recorded. No cross-source minimum or cross-region energy comparison
is allowed. Original-source coordinates remain the immutable physical reference.

The exact shared `CARTESIAN_OPTIMIZER.py` SHA
`d93821da485f8d078e42bf7031018931632858ea720af617f26da2dab77617ee`
performs the search. All 256 atoms move; only instantaneous whole-system rigid
translation and rotation are projected out. There are no source displacement
windows, rigid bonds, selected torsions, new directions or per-atom clipping.
Full-gradient convergence remains max atom force <=0.10 and atom-vector RMS
<=0.03 kcal/mol/angstrom. Initial trust is 0.025 A, maximum atom step 0.10 A.

The shared PhysicalGuard uses the **complete original peptide covalent graph**,
not the QM/MM-truncated ledger. Native equilibrium bond bounds are 0.65–1.35,
expanded only enough to retain any pre-existing source defect. Chirality,
cis/trans peptide basins and worsening heavy/H contacts retain their existing
checks (1.0/0.55 A, excluding graph distances 1 and 2). Source motion is recorded
but does not reject a point. All ten original site oxygen IDs are checked; target
departure beyond 4.5 A from protein N/O or from that declared site is labelled.
Current preparations contain no waters/spectators, and the wrapper fails if that
inventory changes rather than silently inventing new state rules.

State-changing trials that pass covalent/contact guards are evaluated and kept.
The unchanged optimizer stops after accepting such a trial; an energetically
rejected trial stays in the record. State-departed terminals cannot enter the
bound-state paired pool. Budget-limited but physically admitted terminals remain
available with their full force residuals and nonconvergence explicit.

## Finite stages and restarts

`preflight`: saved evidence and four start/guard checks only, zero new calls.

`scout`: only 7CCO/La in the declared region, <=12 accepted steps/40 distinct
composite evaluations. This first-state choice is fixed before reference outcomes.

`run --state SOURCE_METAL`: one explicit state, <=80 accepted steps/200 distinct
composite evaluations total. All runs require a complete, physically admitted
scout from this same bound region/model/start policy. 7CCO/La resumes its exact
scout checkpoint, so scout steps/calls count toward 80/200. Other states start
from their own selected parent points. A converged scout can complete the run
without another call. No loop launches all four states or another region.

Before a search, RUN_BINDING freezes the region/state, parent pool, source and
helper hashes, evaluator identity, settings and exact initial POINT. Checkpoints
are written by the unchanged core. Saved result admission restores every completed
request and all evaluated trials, including rejected points and prior exact cache
hits; failed native requests/checkpoints cannot be retried automatically. Point
identity includes the complete model/environment through V2.

An optional `resume_result` supports an **admitted complete run's** final
checkpoint under the same 80/200 ceiling; it never resets a budget. Arbitrary
interrupted runs without a final complete result are intentionally not accepted
by this draft: their saved checkpoints/raw receipts remain available for a
separately declared failure/preemption recovery. This prevents a stale checkpoint
from silently bypassing a later failed component. No execution failure is treated
as convergence or permission to continue.

## Final paired pool and interpretation

`pool` requires four ordered actual run result slots for the one declared region.
Unavailable/failing slots remain explicit. No cross calls occur unless all four
terminal states are admitted. All original 26 candidate geometries/52 same-region
cells remain. Add up to four terminal geometries, with exact source-local physical
coordinate-byte deduplication and alias retention; do not add all rejected trials
to the candidate menu. Restore all their completed trial caches for reuse.

At most four missing opposite-metal composite cells are evaluated. Their native
components remain bounded by <=4 MACE, <=8 GFN2, <=4 MM configurations/16 group
queries and <=4 new ALPB calls; exact cached components can reduce those counts.
The complete paired denominator is declared before calls. A backend failure stops
new calls and preserves the attempted cell plus all unavailable siblings. No
favorable subset or incompatible state is silently selected.

Source-local menus use the existing frozen `menu()` and native-component
`difference()` functions. Report additional per-metal work from that same region's
parent minimum, its four components, and the resulting Dy-minus-La contrast
change. Preserve all four final force residuals, convergence flags and selection
ties. Do not compare absolute energies across regions, mix sources, claim affinity,
or label a smaller source gap as independent accuracy validation.

## Costs and execution interface

Every allocation records actual new component/model calls, parent/native CPU
seconds, elapsed time, allocation metadata and raw request/point pins. Parent,
scout and current-allocation costs are separate. Optimizer counts are cumulative
within a state's finite budget. Scheduler receipts remain root-owned; application
wall time is not substituted for scheduler cost.

After root binds immutable JSON pins, the CLI is:

```text
python RESPONSE.py --stage preflight --config CONFIG.json --output NEW_PREFLIGHT
python RESPONSE.py --stage scout --config CONFIG.json --output NEW_SCOUT
python RESPONSE.py --stage run --state 7CCO_La --config CONFIG_WITH_SCOUT.json --output NEW_RUN
python RESPONSE.py --stage pool --config CONFIG_WITH_FOUR_RUNS.json --output NEW_POOL
```

All stages require Slurm; molecular stages require an allocated GPU. Root chooses
the actual launcher/runtime from the qualified V2 environment. No launch script,
automatic scheduler, new engine or production-default mutation is supplied.

## Review and planned verification

Only source inspection and AST parsing occur on login. Independent review covers
the actual parent/POINT schemas, full graph/state map, exact binding and cumulative
cache/cost semantics. Root's first saved-data Slurm preflight must replay all 104
parent points/selections and admit four starts with zero molecular calls before
any separately chosen finite scout. A passing script or lower energy does not
establish molecular accuracy or stationary mechanics.
