# Native stopping diagnostic: complete

All 20 calls completed. Tightening only native `TolE` to 1e-10 Eh makes all ten exact cells agree from their loose-cold and loose-pass2 seeds: maximum difference **2.55e-7 kcal/mol**, against the frozen 0.1 tolerance. No input geometry, charge, parameters, mixer, temperature or physical model changed.

The five formerly unsettled La cells and their Ca partners were all retained. P38539 La returned from the pass2 excursion by 4.11304 kcal/mol; Q4W6G0 La shifted by 0.12236 kcal/mol. All calls verified actual XTBRESTART, native mixer and effective TolE. All 20 meet printed RMS-density limits; 19 meet MAX-density limits. The remaining MAX-density residual is 1.2264e-5 versus 1e-5. Native stopping is governed by its own energy criterion; density thresholds alone do not establish an implementation violation.

This supports stricter scalar stopping on these consumed failures. It does not qualify forces, establish a unique electronic ground state, or repair a classifier/reference from ten components. Both seed branches are continuations: a prospective routine supported here costs a loose cold calculation plus a strict continuation. Neither is a fresh one-call calculation.

Job **1211229** used **17 s on 20 CPUs = 340 allocated core-s**, 48 GiB requested RAM, zero GPU time; summed ORCA runtimes were 172.72 s. Historical seed costs are separate. Three real-artifact tests passed in 2.181 s, no skips. Baseline and production are unchanged.

Full energies, per-cell shifts, traces, charges and receipts: `workspaces/strict_native_stopping_20260923/run_v1/{COLLECTION,TRACE_AUDIT_v1,COSTS}.json`; compact pins are in `ARTIFACTS.json`.

The next separately frozen experiment compares fresh strict starts with strict continuations uniformly over the existing 32-source/384-cell pool. It is not another pass on only inconvenient outputs.
