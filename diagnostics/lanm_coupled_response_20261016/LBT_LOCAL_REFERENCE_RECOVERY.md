# Local-reference JSON replay recovery

The original La/7CCO scout 12826 remains admitted and unchanged. The three
continuations 12837, 12838 and 12839 failed before any new molecular call at
the scout's raw-replay equality check. Preserve their original results.

Saved-only Slurm diagnosis 12840 proved exactly one representation difference:
`row.executed_ecp[0]` is the native regex tuple `('La', '', '46')`, while the
persisted JSON value is `['La', '', '46']`. JSON-roundtrip equality is true,
with no canonical differences. No energy, coordinate, state or numerical
value differs. The diagnostic made zero molecular calls.

V2 normalizes the reconstructed native row to its JSON representation before
the unchanged full-row equality checks. Its config explicitly pins the
original scout implementation (SHA256 `7be8b57056d22ca4fdeabb26c991acaaefb21a4c420fb41cac9c7ad57067e784`),
the actual scout result and diagnosis. It does not accept arbitrary historical
implementations. The saved comparator binds an explicit adapter for each
result: the original source for the scout, and the new V2 source for the three
remaining cells. Each result, config and executed manifest must agree with
that declared lineage by verified content hash.

All source coordinates, physical/effective states, numerical settings,
convergence/parser gates, costs, unavailable cells and four-cell denominator
remain unchanged. No scout rerun, changed SCF settings or relaxed tolerances.
Only the original three unexecuted endpoints remain to be evaluated after
root's independent review and fresh versioned submissions.
