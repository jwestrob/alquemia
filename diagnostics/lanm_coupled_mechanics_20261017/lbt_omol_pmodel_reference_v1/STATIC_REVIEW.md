# Independent static review

21 October 2026 UTC. **PASS for the declared one-CPU, zero-call preparation.** No source blocker remains in the separately versioned PModel variant. Actual preparation and a root-bound execution decision are still required before the sole24-MPI endpoint. Review used source/diff, JSON and file hashes only; no scientific imports, numerical replay, model loading or job execution.

## Frozen bindings

| Artifact | SHA256 |
| --- | --- |
| TARGET.py | `27ac4e80b5fb53c049adc6f67822f1629cd3bd7cf372edf1d4708cedcc829856` |
| SPEC.md | `183a7f10841b1e251393f953052b2a2b7da47776806536b2709d6f1484d9e296` |
| CONFIG_PREPARE.json | `a1007083e3d3139848b18bd3157d2cf0746f52e429c42f44dc34872fcfc159ad` |

All14 direct configuration file pins match their files. The failed16602 result `88185237…`, implementation `aa62db4f…`, configuration `1c9f22df…` and input `d482826a…` are explicitly fixed. Preparation binds no future result or execution decision. Launcher/READY execution is separate from this source verdict; READY was not present at this observation.

## Reviewed invariants

- The new recipe must equal the exact failed template after replacing its single `PAtom` occurrence with `PModel`. No SCNL keyword, solver fallback, changed grid, tolerance or electronic-state change is inserted. The original XYZ comment and coordinate serialization are retained, permitting byte-identical XYZ. Existing public assets and source helpers remain pinned.
- Failed-parent identity, status, job, one attempted call and source/config/input joins are checked. Parent physical/positions/ID authority must equal the re-admitted step725 central cell. Whole256 geometry/order, q−2, UKS singlet,1030 electrons, La46 and beta rotation `{514,515,20,1,1}` remain fixed. Only the initial-guess metadata changes. The saved lineage helpers perform no new molecular calls.
- Adjacent external basis copying and verification before and after execution are unchanged. Metadata now identifies the observed external-file La ECP46 and actual SCNL/libXC treatment; it does not repeat the disproven post-SCF assumption. Native execution must still supply the method evidence.
- The existing separate24-rank allocation guard, allocation-derived MaxCore, one worker, fresh-directory policy, preparation-signature match and root-release gate remain intact. Only `prepare` and a single `run` are implemented. There is no loop, automatic retry, alternate guess, additional molecular endpoint or displaced stage.
- Native PModel atomic-density subcalculations are explicitly allowed within initialization and included in allocation/native cost. Surviving task-local `*_atom*` files are pinned after runner return, including a failed attempt. The adapter correctly does not infer internal solver-call counts from file counts. An atomic density's multiplicity must not be interpreted as the molecular singlet being changed.
- Generic native receipt/normal-termination/SCF and finite energy/full256 gradient gates remain unchanged. The forbidden old-guess check now rejects polarized atomic densities; MODEL POTENTIAL and SCNL are captured as observations. The strongest success remains `native_core_admitted_pending_method_review`, so boolean text matches cannot silently release an unreviewed method. Coordinate and energy identity tolerances, force sign, full-vector comparison and frozen-direction projection are unchanged.

## Scientific limits

The PAtom failure establishes unavailable atomic-guess support on that execution path. A PModel initialization is a coherent bounded response supported by prior La/ECP initialization, but diffuse UKS/SCNL success is untested. Changing initial orbitals can select another electronic solution. This variant is therefore not exact historical ORCA6.0.0 reproduction, a stable-state determination or a force-accuracy qualification.

Keep16602 immutable. After an actual central result, inspect PModel execution, beta rotation, SCNL, orbital/auxiliary/ECP assignments, residuals, density/gap/S² and cost before any further decision. No force-model attribution, ALPB validation, preference claim or SP release follows from preparation or scheduler success. A future SP comparison must explicitly match the observed nonlocal Hamiltonian; it cannot blindly remove EnGrad from this deck.
