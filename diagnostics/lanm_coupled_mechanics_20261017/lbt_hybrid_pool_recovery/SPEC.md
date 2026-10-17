# One isolated-component recovery of the complete hybrid pool

## Question and fixed inputs

Complete the original 104-cell, 26-geometry, two-region LBT comparison after its
single isolated GFN2 SCC failure. This is a numerical recovery of the declared
Hamiltonian. It is not a new region selection, structural search, parameter fit,
or test selected by a favorable metal-preference sign.

Job 12992 retained 66 complete cells (eight origins and 58 new composites), one
attempted unavailable cell, and 37 unattempted cells. The failed cell is
`cartesian_03 / 7CCN / baseline / La`, original `point_059`. Its embedded GFN,
MACE and MM results already exist. All 52 global fixed-charge ALPB results exist
at the exact physical coordinates in the admitted Module2 pool.

The four-call diagnostic 13009 admitted both requested damping values at the
failed isolated state and both predeclared controls. Its original fixed gates
are energy agreement at most 1e-7 hartree and every native-QM/physical gradient
component at most 0.005 kcal/mol/angstrom. Charge vectors remain descriptive.
The primary is the predeclared `failed_005` result, never the lower-energy result.
Native adaptive behavior remains unchanged; requested damping is not claimed
to be a constant coefficient throughout every native SCC iteration.

## Saved-only preflight

`RECOVER.py --stage preflight` imports the exact frozen pool preparation, including
the original scout, eight origins, remaining derivative qualifications and
Module2 source/state/solvent admissions. It replays the ordered 104-cell ledger,
all 66 accepted raw points and the original unavailable comparison.

It then reruns native admission from the four existing diagnostic outputs,
checks their exact staged XYZ/parameter/control bytes and commands, and replays
all three agreement checks. The failed and immediately preceding fixtures must
be the last two actual 12992 requests. No native calculation is launched.

One new immutable POINT is assembled using the archived failed-point map and
exact V2 component arithmetic: local MACE vacuum, embedded-minus-isolated GFN2,
the four-group MM ledger, and exact cached full-system fixed-charge ALPB. The
only replaced component is the primary accepted isolated result from 13009.
The old failed POINT and isolated output remain unavailable and unchanged.
Frozen V2 restore and the qualifier's ALPB replay then validate the assembled
POINT, full derivative closure, state, coordinate identity and cache identity.

`SOLVER_EXCEPTION_CERTIFICATE.json` records the exact failed record, four-call
qualification, three unchanged agreement gates, primary and secondary native
receipts, all reused components, Hamiltonian cache identity, and assembly source.
Both raw POINT and point-level provenance retain its pin. The certificate says
same Hamiltonian and **different solver protocol**; it does not establish a
stable electronic ground state or change the global solver default.

Admission requires exactly 67 complete cells plus the original ordered 37
unattempted cells, all 52 solvent states admitted, zero molecular-call counters
and no MACE model initialization. An unavailable raw input, mismatched fixture,
failed agreement, altered certificate or changed derivative fails closed.

## Finite continuation

Root binds the successful preflight RESULT in a new immutable run configuration.
`--stage run` replays the certificate and expected assembled POINT exactly, then
uses the **archived preflight POINT**, preserving its provenance on every replay.
Only the original 37 unattempted slots may call unchanged HybridEvaluator V2.
The ceiling is 37 MACE calls, 74 GFN calls, 37 MM configurations / 148 group
queries, and zero ALPB calls. Stop on the first new backend failure; retain the
attempted raw request and remaining unavailable denominator. No automatic retry,
global damping change, increased iterations, altered temperature/accuracy,
changed state, or different physical/finite-difference tolerance is included.

Completion requires all 104 cells and all four original source-local paired
menus. The exact original summary functions compute works and partition effects.
No unlike-Hamiltonian absolute total is compared and no cross-source minimum is
introduced. No optimizer or DFT call is part of this recovery.

## Saved-reader contract and costs

The module exports `PROTOCOL`, `prepare`, `summarize`, `menu` and
`cell_from_point` with the frozen pool signatures. The configuration retains
`qualification_implementation`, `qualification_result`, `module2_pool_result`
and `plan`. `prepare` certifies the archived replacement before returning, so
the existing RESPONSE/comparison can subsequently use unchanged qualifier
restore on the actual pinned replacement and all ordinary cells. No evaluator
code or cache engine is replaced.

The result keeps the original 104-cell schema, plus the certificate, archived
replacement, preserved failure, remaining-slot ledger and separate parent costs.
This allocation's `actual_cost`, requests and points include only its attempted
37-cell continuation. The saved assembly costs zero molecular calls. Prior
12992 and 13009 costs remain separately pinned and are not counted as new work.
All execution, including saved numerical replay, uses Slurm. Source checks on
the login host are AST, metadata and hashes only. Root owns submission.

## Commands for root staging

Run with the pinned hybrid Python/runtime and the existing qualified H200
launcher environment. Preflight needs no GPU/model call but imports the same
qualified dependency environment. Both output directories must be new.

```text
python RECOVER.py --stage preflight --config CONFIG_PREFLIGHT.json --output preflight
python RECOVER.py --stage run --config CONFIG_RUN.json --output execution
```

The run configuration differs only by the actual `recovery_preflight_result`
pin. Verified identical source/config snapshots are accepted by SHA; scientific
input bindings and all stored evidence remain immutable.
