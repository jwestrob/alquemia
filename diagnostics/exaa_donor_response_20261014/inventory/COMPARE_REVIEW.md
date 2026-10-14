# Static comparison review

Reviewed `reference/COMPARE.py` against the actual pinned inventory, `reference/PREPARE.py` and `reference/REFERENCE.py`, including the reused native collector. No scientific execution, imports or tests were performed.

Inspected comparison SHA256: `a29ac2e2cff1acae27bb15626593757eabe8b12723f84181d17ff3674eaa8711`.

The joins match the implemented schemas: inventory rows identify the exact original sources, their metal/candidate matrices expose `MACE_eV`, and native collections provide manifest pins, task IDs, `complete_diagnostic` status, energy in hartree, output/execution pins and numerical checks. Both full eight-cell scope and partial scout collections are representable.

Individual work is correctly signed as energy at `adaptive_La` minus origin energy. MACE eV and DFT hartree are converted once to kcal/mol. The differential is Ca work minus La work; method discrepancy is DFT minus MACE. This is the common physical La-proposal path, distinct from the operational score's independently selected Ca/La proposals. No frozen-ALPB term or calibration band is mixed into the vacuum comparison.

An individual DFT response requires both of its accepted endpoints; a differential requires both individual responses. Uncollected and failed cells therefore cannot become an accepted partial contrast. Conflicting duplicate cell records and mismatched manifests fail explicitly. The full cell dictionaries retain electronic warning lines, populations and the statement that density/state qualification remains unestablished. Solver convergence does not establish the same unique electronic state along a path.

Two small follow-on gaps were sent to root before execution:

1. The frozen plan asks for target-minus-control differences, but the inspected version emits only the two per-source records. Add that derived comparison with independent per-metal and complete-differential availability gates.
2. A requested collection file that does not exist currently raises before producing a result. If collection follows failed/preflight-only jobs with `afterany`, record missing paths and retain the declared unavailable cells. A malformed or incompatible existing collection should remain an error.

No other concrete schema or algebra blocker was found. This is static review, not a claim that real reference execution or scientific comparison has passed.

Root follow-through: both identified gaps were addressed before scout collection was queued. Missing collection paths are explicitly recorded, and target-minus-control differential work is emitted only when its required native cells are available. The new join awaits actual scout outputs; no executed comparison is implied.
