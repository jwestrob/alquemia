# Independent recipe and source review

21 October 2026 UTC. **Static review passes for the declared zero-call preparation.** No source blocker was found. Actual preparation admission and a separate root release remain necessary before the sole central native call. This review used source text and metadata/hash inspection only; it imported no scientific modules and executed no molecular or numerical replay.

## Frozen evidence

| Artifact | SHA256 |
| --- | --- |
| TARGET.py | `aa62db4f37de3c50765ad81115dd75cb4aef2a99d92b4465d9bbe59df9a4a74f` |
| CONFIG_PREPARE.json | `f8a105695715b69cd0e8fe06ff10b5a032c55d0a97ace77745bc3b6c75352e74` |
| SPEC.md | `0f11bb6a187391a5284fd8962fd9738672401b1c7cc5c370c1418effb1bc0ae4` |
| TARGET_RECIPE_EVIDENCE.md | `cb20795b46d2d413b1d494cd6462d68ca0754638942a8e4993267ef39d8b30bb` |
| READY.json | `4ba6e0d37341467aa07ddcf58a8dfd48a1a6ae37e862e33096144f157d74d54c` |

All eight direct configuration file pins match their files. READY binds the same source/configuration and the inspected one-CPU preparation launcher `b41e24c9d46adcbeac6bf01c0e372a6b9a6e5131c89c5fcf9d51b3153fa92c80`; its argument list supplies the script after the launcher's explicit interpreter, avoiding a duplicate interpreter. Actual execution remains outside this static verdict.

## Scope and interfaces

- Saved input admission retains the exact step725 whole256 La geometry, atom order, charge −2, physical singlet and 1030 explicit electrons. The old direction/reference helpers supply saved evidence only; their historical 12-rank declarations do not govern the new execution. No old chemistry is repeated by preparation.
- The new deck extracts the pinned published functional, basis, simple keywords and blocks without importing the upstream calculator. It retains the explicit UKS beta rotation `{514,515,20,1,1}`, Tight/300, DIIS/NOSOSCF, thresholds, DEFGRID3 and output controls. Explicit PAtom/NoAutostart protect a fresh unseeded start; no PModel, D4, SCNL, solvent or point charges are introduced. The public wrapper's NONBO/NONPA convention is preserved.
- The complete external orbital basis is copied beside the input and XYZ. Its hash and adjacency are checked before execution and again in native admission. The inherited runner executes in that directory and preserves the input except for allocation-derived runtime directives. The generic runner does not independently bind this auxiliary file; TARGET supplies that missing binding. GTOName does not establish the executed built-in ECP/def2-J assignments by itself.
- The new allocation guard requires 24 granted CPUs and 24 local MPI slots, one worker and allocation-derived memory with the existing reserve. It calls the generic runner with `nprocs=24`; it does not call the old 12-rank allocation wrapper. Preparation has no allocation-dependent native call. Fresh directories, matching prepared task signatures and an actual preparation/root-decision binding constrain execution to one central EnGrad.
- Native core admission checks the generic receipt and exact artifacts, normal termination/SCF convergence, version, charge, UHF formalism, electron count, full basis dimension, La ECP46, MaxIter/DIIS/SOSCF, a finite final energy and complete finite256 gradient with the established coordinate/energy identity tolerances. The old RHF/PBE0 scientific parser is not reused. Gradient-to-force conversion and comparison use the full physical ordering and the previously frozen direction.

## Required interpretation after the scout

The strongest programmed success is deliberately `native_core_admitted_pending_method_review`. Functional/rotation/guess pattern flags are observations, not complete method admission. Actual orbital and auxiliary assignments, ECP details, post-SCF VV10 energy/gradient evidence, XC/COSX/nonlocal grids, SCF residuals, density, frontier gap and S² still need saved-output review. Missing or unfamiliar text cannot be filled from an input echo or repaired by rerunning chemistry. The generated force comparison is explicitly descriptive and pending that review.

This is a useful consumed-geometry test of comparator-method mismatch, with ORCA6.1.1 versus historical6.0.0 explicitly unresolved. UKS rotation does not establish a stable physical singlet. A smaller force discrepancy supports a contribution from the changed comparator; persistence does not uniquely identify a learned-model defect. Absolute energies across methods are not deformation work. No central result alone releases displaced calls, qualifies derivatives, validates ALPB or establishes metal preference. The failed PBE0 directional gate remains unchanged.

Primary recipe and interpretation sources are pinned in TARGET_RECIPE_EVIDENCE.md: [FAIRchem calculator](https://github.com/facebookresearch/fairchem/blob/d4e10d0f4e87593922fda1d1ffb7dd28e504cadd/src/fairchem/data/omol/orca/calc.py), [singlepoint wrapper](https://github.com/facebookresearch/fairchem/blob/d4e10d0f4e87593922fda1d1ffb7dd28e504cadd/src/fairchem/data/omol/orca/recipes.py), and the [ORCA6.1 nonlocal-correlation documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/dispersioncorrections.html#basic-usage).
