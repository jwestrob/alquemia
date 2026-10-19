# Independent static review: preflight and one whole-GFN scout

19 October 2026. **PASS for these two implemented stages; no static blocker found.** Actual Slurm preflight must pass before the single scout. This review performed source, metadata and file-hash inspection only: no scientific imports, numerical execution, model loads or submissions.

## Frozen artifacts

- BASE_DELTA.py: `f760cd799f326fe193d2a2809e4434000803e04438a3358de0210d829d32a5d1`.
- SPEC.md: `269793a3ca56f6caca69f4cf491af25b0f7d819f82c1ff0b104fbbbf833295e3`.
- Preflight CONFIG.json: `8dc6cdf05bd84874fc3a6a2db016f544cc4d9e59583ed2b5be47fe33a59da572`.
- Preflight RUN.sbatch: `899397ca89842184aa41d88e3e923f6a4bb641ee293b3a616276c58a8bd30714`.
- READY.json: `d751d6bd5a1a58341f1634dc61fffc6956a9406eb30d1fd24dfa2617048f95be`.

All exact path/SHA records in the final config were independently hash-verified. The root declaration and inventory remain those reviewed in [SPEC_REVIEW.md](../whole_base_delta_design_v1/SPEC_REVIEW.md); the earlier design REVIEW.md remains immutable.

## Scope and execution

The eight-cell denominator and ten-call candidate ceiling are preserved. This executable exposes only `preflight` and `run-scout`; its only molecular-call site is one whole-256, vacuum 7CCO before-La GFN evaluation. Development, missing regional calls, assembly and challenge execution are absent and require separate reviewed source. The launcher requests one CPU, no GPU, no requeue, and scheduler default time/memory, with native thread counts fixed at one.

The Slurm guard precedes scientific imports. Inspected helper bodies confirm that preflight and `RECOVER.certify_scc` perform saved admission, map and certificate replay without invoking the recovery run path, HybridEvaluator, a MACE constructor or a native molecular call. All eight physical/state/map/solvent declarations are prepared before any scout. A scout requires the actual passed preflight, identical scientific binding and prepared signatures, and a fresh output directory. The unbound scout template cannot satisfy that requirement.

## Provenance and derivatives

Seven saved regional pairs are joined to their exact POINT, task, XYZ, physical source and model/state records. I independently checked all seven native/MACE XYZ hashes and physical-source joins in the archived metadata; they match. The 7CCN after-La exception replays all four original SCC controls and its exact primary certificate, retaining different solver protocol and the unavailable original failure. It cannot enable damping changes for the new scout.

The original after-Dy inventory's `physical_geometry_parent` alias and list-shaped supplemental collection are handled explicitly. All eight ALPB receipts use the unchanged whole-reference raw admission against source-specific coordinates, symbols and fixed charges. No solvent evaluation occurs.

Archived maps remain authoritative. The existing V2 comparison permits only its already-declared bounded derived-Jacobian roundoff; coordinates and identities remain exact. Raw MACE forces become gradients with `-f*EV`; isolated GFN gradients use `(g*HA)/BOHR`. Both are separately scattered through the retained and omitted cap anchors, with GFN negated before adding the two full-256 components. Field gradients are zero. This is the declared new component assembly, not a claim that an old hybrid-gradient bitwise replay passed. No MM or embedded-GFN increment enters it.

## Native admission and limits

The scout stages the full qualified parameter text, points XTBPATH to that fresh directory and uses the pinned xTB executable/runtime. The command preserves charge −2, effective UHF 0, 300 K, accuracy 0.0001, 500 iterations, no restart, and native population/analytic-gradient outputs. Original admission checks normal termination, SCC/settings, populations, scalar agreement, full gradient layout, atomic order and coordinates; finite full-gradient checks remain explicit. Physical MACE La1/Dy6 is separately recorded from effective GFN state.

`scout_passed` requires native admission and measured native wall time at most 600 seconds. Excess cost or failure blocks continuation and retains the attempted slot; remaining slots stay not requested. This is a continuation limit, not a scheduler timeout or a retry instruction.

The mandatory base-only comparator and both-metal 1 kcal/mol development rule are bound in the declaration/config, but their collector is not implemented here. No candidate work, force accuracy, finite-difference qualification, physical-state validity or challenge readiness follows from this static pass or a future scout pass.
