# Native-key source-map recovery after zero-call12711

Slurm12711 failed in input admission, before compilation or any molecular call.
All four declared cells remain unattempted. The scheduler records3seconds on1CPU;
the application recorded3.2901902499143034seconds. No molecular work is duplicated.
Original source2eda8f20, configd2e5e01d and failure result
`0c336dbd6ff8c1071c93a002618d2ff055195a74f2e261a2273167f3face7568`
remain immutable in `workspaces/lanm_coupled_response_20261016/lbt_matched_model_v1`.

The exact failed lookup was `byid[a['id']]` while iterating
`PREPARATION.source_to_native.values()`. Both actual source maps use native
physical ID`A/1//H` as the key, while the deposited record stored as its value
retains ID`A/1//H1`. The physical `A/1//H` record carries that same value in
`deposited_source`. Thus the dictionary key supplies the native join; the value
is deposited provenance. No atom is missing and no chemical alias is invented.

V2 iterates `(native_id, deposited_record)` from `.items()`, requires the key in
the physical map, requires exact equality of its `deposited_source` record, and
retains the original element equality and1e−12Å coordinate replay tolerance.
This strengthens the provenance check while correcting the key/value schema.
Coordinates, native FF charges, terminal chemistry, atom order, all four MACE
records, Mex assemblies, native solvent functional/build, four-cell denominator,
common-pool rule and exchange signs remain unchanged. The original scientific
plan remains pinned. V2 additionally pins this note and the exact zero-call failure,
and requires the previous scientific config to match apart from implementation.

Static review covers the actual remaining source task/row/map schemas; there is
no adapter import, numerical replay, compilation or submission in this repair.
Root owns any separately staged V2 execution. No source or failed record is
overwritten, and no molecular retry is hidden.
