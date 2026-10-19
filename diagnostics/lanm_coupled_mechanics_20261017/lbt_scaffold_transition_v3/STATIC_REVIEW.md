# Independent V3 nullable-task review

Reviewed 2026-10-19T07:54:19Z by `/root/mex_health_review`.

**PASS for the narrow metadata correction; actual runtime admission remains pending.**
This review used source, saved JSON/text, diffs and file hashes only. No numerical
calculation, scientific import, model call or submission was performed.

| Reviewed file | SHA256 |
|---|---|
| ANALYZE.py | `d736a58917262f9977c6a3152b5e93feaa2183f9abb55689fef97cc657fef515` |
| CONFIG.json | `bf88807d5acd6d5e2bd20ee0859b101cc871cbb850f964768dc9a756cbd1fd22` |
| SPEC.md | `c55ada879b496434a238013b91584e60e8c5c68539a3e5d25098aec6d1858d40` |

All embedded CONFIG artifact pins passed byte-hash verification.

## Original semantics and actual failure

Actual V2 job14786 result `a66244f0…` stopped with zero MM contexts/queries when
`verify(task['ALPB'])` received null. The existing native after-Dy task and its
original signature both deliberately retain `ALPB: null`.

The pinned original CONTINUATION.py (`06e40d37…`), lines173–180, admits the later
completed whole-model after-Dy ALPB receipt (`26f188d8…`) through a temporary
copied row, explicitly preserving original cells/full manifest. Its
`prepare_task` also preserves the original task signature. The actual final
model row and reference solvent admission contain the same completed receipt;
the native vacuum task's null field does not mean that final solvent is missing.

V3 requires precisely after-Dy when the task field is null, requires the original
signature to remain null, and preserves equality of the completed model ALPB
receipt and reference solvent receipt. Every non-null task ALPB pin is still
verified. XYZ, input, physical-parent, positions and native receipt/output/engrad
pins remain mandatory. Saved metadata inspection found those other file fields
to be objects in all four native rows; no additional nullable-file mismatch was
identified in this inspected path.

## Scope and limits

The source diff otherwise only versions the protocol and binds the preserved
zero-query V2 failure. The completed raw-ALPB gates, V2 exact conversion orders,
force comparisons, MM attribution and report are unchanged. The CONFIG diff
preserves all four inputs, 187 interactions per state, 36 contexts/816 queries,
original closure thresholds and the unchanged one-CPU launcher. Both failed
attempts and prior reviews remain intact.

No outstanding static blocker was found. Complete numerical replay and all MM
closure checks still require the actual Slurm result. This correction provides
no new physical accuracy, electronic-state or force-qualification claim.
