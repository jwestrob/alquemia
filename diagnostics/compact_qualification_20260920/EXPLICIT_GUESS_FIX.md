# Explicit native-xTB guess override — 2026-09-20

The MORead simple keyword is also overridden in the actual native-xTB input:
job1203284 prints "MOInp is set but Guess!=MORead. MOInp will be ignored!" and
again starts from SAD. No orbital-seeded solver agreement was tested by that
attempt. Preserve the failed output and cost.

Use the manual's explicit `%scf Guess MORead` and `MOInp "seed_source.gbw"`
instead of relying on the simple keyword. Request this after the native solver
selection in the block, with no NoAutostart simple keyword. This is another
technical input correction of the same eight cases, not a changed Hamiltonian
or target result. Exactly eight low-level calls, unchanged states, thresholds,
temperature, parameters and resources. No repeated primary panel or calibration.
Require actual orbital restart evidence and absence of the ignored-guess warning.
If native xTB still overrides this control, document that backend restriction;
the completed six-case tighter-native convergence test remains separately useful.

The user has authorized pursuing the candidate; this corrects the agreed restart
mechanism. It does not establish or claim another biological result.
