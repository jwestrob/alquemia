# Four-cell LanM runner review

Reviewed the research-only `metal_environment_lanm_reference.py` before preparation completion; no molecular execution.

Fixed concrete missing invariants:

- A/B permanent-charge values and ordering must match exactly; geometry changes alone are allowed. Previously each array was only checked individually for shape/finiteness.
- The execution path explicitly gives `METAL_ENV_WORKERS` the manifest concurrency, so allocation-derived per-rank memory does not silently use the renderer's default4 when actual concurrency differs.
- Source identity must match the pinned preparation config; execution runner/renderer must match the pinned implementation entries; resource layout must contain positive integer ranks and1–4workers.

Existing useful behavior retained: each La/Dy pair has exact common coordinates and identical nonmetal membership; A/B physical composition and charge are fixed; input/state/environment/software pins enter the cache; partial native jobs remain visible and are not silently overwritten; a missing required endpoint leaves the response unavailable. The comparison is `Dy(B)-Dy(A) - [La(B)-La(A)]` with one Hartree-to-kcal conversion, not an affinity or inherited PQQ classification. Point-charge gradient files are preserved for later physical projection.

The runner still uses the existing unqualified default-grid/SCF reference recipe, deliberately separate from claims of qualified mechanical relaxation. This review introduces no new electronic method or numerical-policy choice. Full boundary forces and solvent remain outside this four-cell electronic-component experiment.

Syntax compilation passed. Once the real Hans8DQ2 preparation became available, seven actual-fixture tests passed against frozen INPUTS SHA c5886704ee541752ec2adfe7ada0b5f8b2c0df6014d5ed8676b74447d1141c3f. The intact expanded core has La806/Dy833 explicit electrons. Explicit corrupted-copy tests reject singletDy, changedcharge/parity, false electronic metadata, mismatched paired coordinates, changed permanent-charge inventory and stale hashes. Endpoint all-electron count, oxidation hypothesis and2S metadata are now checked when supplied. No molecular calculation or new execution manifest was made by this review; root owns the actual four-cell dry-run/submission. No dummy molecular configuration or successful endpoint was substituted.

Run: `python -m unittest discover -s tests -p test_metal_environment_lanm_reference.py -v`.
