# Four matched strict scalar comparators: complete

All 48 fresh scalar calls and all four whole pools completed. Each retains its correct classification under the frozen strict32 fresh reference. No new calibration, MACE evaluation, geometry search or DFT calculation ran.

| Exact noncanonical source | Known class | Original loose R | Fresh strict R | Strict call |
|---|---|---:|---:|---|
| A0A3F2YLY8, Ca sample3 | La | −405455.470283 | −405455.469245 | La-supported |
| A0ACD6B9F2, La sample4 | La | −405439.809173 | −405439.810834 | La-supported |
| A8R3S4, La sample1 | Ca | −405466.922180 | −405466.928339 | Ca-supported |
| A8R3S4, La sample3 | Ca | −405467.346652 | −405467.349478 | Ca-supported |

R is the composite E_Ca−E_La in model kcal/mol. Values are rounded for display only. All original/adaptive_Ca/adaptive_La matrix cells remain available, with archived MACE energies and unchanged source coordinates/charges. Mathematical and operational calls coincide. These are consumed folds from correlated protein groups, not four new independent biological observations.

Every output confirmed fresh SAD/NoAutostart, native mixer, TolE=1e-10 Eh, effective state and unchanged parameter export. All use 300 K electronic temperature, MaxIter500 and one rank. No failed attempt or retry occurred. Fresh strict numerical qualification comes from the earlier full32 experiment; this four-source completion applies its frozen reference without recalibrating.

Job **1211311** used **51 s × 32 CPUs = 1,632 allocated core-s**, 64 GiB requested RAM, zero GPU time. Summed ORCA runtimes were 944.007 s. Two real-artifact tests passed in 0.734 s, no skips. The allocation includes execution and collection; local preparation/tests/report are separate.

Exact collection: `workspaces/strict_native_comparator_20260923/run_v1/COLLECTION.json`, SHA256 `937e468e052d0215c50b826dc0a1dde400cf03ab0ad87e33c5ed6db22437542c`. Reference: strict32 `REFERENCES.json`, branch `fresh`, SHA256 `171bd31d466ff97ef6073ce23286af8519ef23690f7d6ad23bb444cab30ecd5c`.

These completed pools supply matched scalar controls for the envelope experiment and exact reuse for the separately executed full225 transfer. The earlier tenfold ledger and production scorer remain unchanged. No further molecular work is pending in this branch.
