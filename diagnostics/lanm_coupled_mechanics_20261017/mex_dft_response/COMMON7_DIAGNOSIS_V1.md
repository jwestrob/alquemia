# Common7 saved POINT replay diagnosis

This is a zero-molecular-call diagnostic of preserved preparation 13063, which
passed the old terminal physical/full-force replay and failed the new adapter's
exact common7 raw-cell assertion. Frozen adapter bd29908d is unchanged.

Run once in Slurm using the original MACE environment (NumPy 2.4.6), the reviewed
Gemmi-only overlay, and the same host/thread settings as preparation 13063.
Root binds and submits. The diagnostic imports the frozen adapter's input path,
then uses the exact evaluator constructor, common5 baseline, both search
admissions, and saved cross-request imports in the same order. All those gates
remain active. Both evaluator call methods are replaced by throwing tripwires;
no model is initialized and all component counters must remain zero.

Inspect every one of the seven coordinates crossed with La/Dy. Record the exact
archived versus current cache POINT pins, native energies, typed POINT leaf
differences, physical file hashes/JSON/IDs, and assembled-gradient pins/JSON.
Also record whether replacing **only** the current POINT's `physical` child pin
with the archived child pin would make the entire POINT exactly equal. This is
a diagnostic comparison on a temporary copy: no cache or admitted row is
changed, no repaired POINT is published, and no numerical tolerance is used.

Retain all 14 rows even when the raw-cell assertion differs. Report the mismatching
slots explicitly. An observed mismatch exits nonzero with `diagnostic_mismatch`;
this status is evidence, not preparation admission. Earlier errors remain
`execution_failure`. The output always declares `replay_admitted: false`.

The known first origin-La join currently differs at `physical.path` only; its
physical SHA is identical. This observation does not predeclare the other 13
joins or authorize a broader path-normalization rule. Root reviews the actual
full diagnostic before any versioned repair.

CLI: `DIAGNOSE_COMMON7_V1.py --config CONFIG.json --output FRESH_DIRECTORY`.
No new DFT, MACE, GFN, MM, ALPB, geometry search, changed state, or source gate.
