# Strict225 complete

Jobs1211337–1211340 and observer1618192 are terminal. All2400 fresh scalar calls
passed; all96 exact reused cells passed. No molecular retry, new MACE, search,
DFT, protonation or preparation chemistry ran. Do not rerun these manifests.

Result:207 correct /0 wrong /1 inconclusive /17 original exclusions. Only Q4Ca2
changes from prior precision, inconclusive→Ca; this repair was already seen in
the reused strict32 development set. Both old released errors stay repaired;
all200 freshly rescored pools retain their prior decisions. C5AXV8 La3 remains
inconclusive. All94 available three-fold aggregates remain correct.

Full comparison: `workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json`
SHA6ed4c1f4acf4359010922ff3f088c5b36105caf52dc92c734e27dd7e35d0df46.
Fresh strict32 reference171bd31d… remains frozen. No fit occurred here.
`FINISH.json` records successful automatic comparison/report/tests and no
collector recovery. Six final tests pass, zero skips.

New cost168352 allocated core-seconds /18960.007 accounted CPU-seconds /0GPU.
See REPORT.md for all counts and INTERPRETATION.md for useful conclusions,
measured throughput/overhead and limitations. ARTIFACTS.json pins the source,
outputs, receipts, timing and tests. Production is unchanged; root owns the
separate strict-static ablation and motion-envelope decision.
