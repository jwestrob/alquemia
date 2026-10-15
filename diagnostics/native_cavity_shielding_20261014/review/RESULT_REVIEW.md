# Independent result review — native cavity scout 4562

**Close this finite exclusion test as a remedy for the ExaA source discrepancy.**
Adding the matched protein volume produces nearly the same metal differential in
both origins. This does not justify a new classifier correction or calibration.
It also does not establish that protein exclusion is irrelevant in other states.

Read-only review of actual saved outputs; no calculations, tests or submissions
were performed by this reviewer. Values below are copied from the collection.

| Saved electrostatic transfer | Eh |
|---|---:|
| Control, Ca minus La | +0.008018277707476207 |
| Target, Ca minus La | +0.007867401972925958 |
| Target minus control | −0.00015087573455024872 |
| Target minus control, shape part | −0.0000026368921155823415 |
| Target minus control, cap-volume diagnostic | +0.00000007680059671555384 |

All 12 declared native updates have actual successful execution receipts. The
four-component local replay is exact; electrostatic translation, rotation and
permutation checks pass, as do component/block/repeat and per-endpoint precision
gates. Both paired precision flags are true. **The aggregate contrast flag is
false:** its conservative transfer/cap-volume bounds are
`1.2399271007806734e-6`/`1.5402942086342367e-6 Eh`, above the frozen `1e-6 Eh`
ceiling. Matrix completion does not override this qualification limit. No
threshold relaxation or repeat for the purpose of obtaining a pass is warranted.

The small cap-volume sensitivity is not evidence that charged caps are harmless:
the physical/cap cross and cap/cap terms remain separate contributions, and the
zero-cap vector is not an electronic state. The primary transfer excludes HB,
SASA and reference terms; in particular, native HB changes substantially on
burial. It must not be called the total solvent correction. This origin-only
experiment retains vacuum populations and adds no protein permanent field or
electronic response. It does not test cavity effects on accommodation, and
unequal Ca/La fragment charges and native radii contribute to the common shift.

The actual runner records 134.5394701550249 wall seconds, 538.1578806200996
allocated CPU-seconds, 130.92 measured child CPU-seconds and 7,545,856 KiB peak
child RSS. The parent's terminal scheduler boundary is 136 seconds on 4 CPUs.
These are different timing boundaries; neither is total campaign cost or a
parallel speedup. There are zero new SCC, MACE or DFT calls.

## Evidence pins

All files are under `workspaces/native_cavity_shielding_20261014/native_v1/execution/`:

- `FINAL_COLLECTION.json`: `dba0b62ca02f9f21cda82b934fcbdd42d1d92768b96f5d30a2b09063ccc965cf`.
- `SCOUT_QUALIFICATION.json`: `8135aa5dd8360ebadabbde1dbb10b935c13b78bb7a64f6f1e80091c245f6aae9`.
- `EXECUTION.json`: `69939174c12721f145777de0d2848a2378d57989d1d29c4fb054784ad2433625`.
