# Independent static review

19 October 2026. **PASS for the declared one-CPU saved-only analysis.** No
material blocker found. Review used source inspection, saved JSON identities
and hashes only; no scientific imports, numerical replay, tests or submissions.
Actual execution remains untested.

## Findings

- All 16 direct CONFIG pins verify. The original search result binds the exact
  optimizer and `scaffold_accommodation` helper. Actual search/checkpoint,
  collector/endpoint, physical IDs, state, common solvent and force-report
  metadata agree with the adapter's literal schema checks. Collector 15173's
  four rows also equal their pinned admitted singleton rows, and all four
  force tables retain the full physical IDs and declared availability.
- The adapter uses the terminal model gradient, all eight stored history pairs
  and saved trust 0.1 Å. It reproduces the original order: rigid projection of
  gradient, L-BFGS construction, rigid projection of direction, the original
  projected-slope fallback/history clear, maximum-atom zero check, trust scaling
  and raw-model-gradient descent check. Reference information does not choose
  or alter the direction.
- Model and reference solvent-conversion arithmetic remain separate and match
  their original producers. Exact raw-component and saved-force-table replay
  precedes the audit. Full gradient signs are correct: force error is
  `F_model - F_reference = g_reference - g_model`. Both local derivatives use
  their full gradients along the same direction. The dimensionless-alpha and
  normalized-distance units are stated correctly.
- Full atom/residue summaries and deterministic error rankings are descriptive.
  They introduce neither a new force-accuracy threshold nor a residue-energy
  decomposition. Native SCF flags, state limitations and missing Dy model
  forces remain explicit.
- Scientific imports occur after the Slurm guard. Imported helper paths expose
  no engine call in this execution path. There is no search invocation, model
  load, new coordinate trial or molecular evaluation. Fresh outputs and failure
  records are retained. READY and the inspected launcher request one CPU, no
  GPU, one thread per numerical library, default scheduler limits and no requeue.

## Interpretation boundary

The original optimizer reaches its 80-step limit before constructing another
direction. This is therefore the direction its unchanged formula **would**
construct if continuation were separately permitted, not an actual proposed or
executed next step. The adapter and SPEC preserve that distinction. A negative
reference local derivative would not guarantee finite-step lowering, establish
force accuracy or authorize continuation. A runtime replay failure must remain
visible; this static pass does not authorize tolerance changes or retries.

## Frozen pins

| Artifact | SHA256 |
| --- | --- |
| ANALYZE.py | `be75686251e2ae467a176b44c9b114c5f36ea23bc74585ab7c356c9ba3658f16` |
| SPEC.md | `cf34e0a7a54db1bb9f070aa4d6f09c998087cf0f840f8d4154e2542726e3e901` |
| CONFIG.json | `b49891153ab69d4a55874276647d8ebc759395126f7cc82e671726629139f8c1` |
| READY.json | `4800f42a469797103329022dd0505e58742ecbefd2986dbe433ce1bf06658ed7` |
| Collector 15173 | `fff95c2f235947acb24e91991915b320731992f48c2b3d6c9ca3912a6e927f38` |
| Terminal checkpoint | `4e3a612380f93e9ba7f374636ed883a659c8b0ade6b35e2f7099f5d7643e9ecf` |
