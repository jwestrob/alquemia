# Zero-field counterpart ready for root release

Frozen source `RUN.py` **fcf72c59**, specification **8c19e258**. Actual preparation v2 admits the 38 saved energies and constructs all 78 inputs with exact coordinate strings, inherited parameters and permanent-frame cancellation; result **545d0d77**, manifest **1f69bc6b**. The 109 saved/synthetic controls pass in 16.36 seconds wall / 15.42 CPU seconds on one thread, with zero native calls. Control receipt **544b4d97**. The initial parameter-comment parsing failure and its source/config remain preserved separately.

`READY.json` **6b5d7c35** contains exact paths, hashes and command arrays. Root should copy `CONFIG_RUN_TEMPLATE.json` into workspace `run_v1/CONFIG_RUN.json`, bind a reviewed release through `root_release`, and retain the actual preparation pin. `ROOT_RELEASE_TEMPLATE.json` gives the gate schema; root supplies the independent review and removes its template-only designation. The fingerprint excludes only `root_release` and `preparation_result`.

The launcher already executes `/usr/bin/python3`. Its positional arguments therefore begin with the absolute `RUN.py` path, then `run --config ... --output ...`. The prepared launcher uses one four-CPU task, four single-thread workers, and scheduler-default memory/time. The native ceiling is 78 analyze calls; there are no gradients, finite-field calls or internal Cartesian finite differences.

Saved collector interface:

```
admit_result(config_dict, manifest_dict, result_dict)
    -> {rows, model_cycles, candidate_result, new_native_calls: 0}
```

Every result row has `moments` containing energy/components, both native induction summaries, all per-site induced vectors, their sums, printed COM total moment, corrected laboratory-origin total moment, raw output pins and the old-energy replay difference where available. `model_cycles.parents[gid][C|P|B]` contains primary induced-vector differences, corrected-M differences and prospective precision bounds. All vectors are Debye; `READY.json` records the exact native conversion to atomic dipole units.

No molecular run or submission has occurred during this preparation. This diagnostic compares the unchanged failed-science 17676 model; no candidate promotion or revised target follows from preparation success.
