# Executed checks

2026-09-16 software checkpoint: **37 tests run, 36 passed, one explicitly skipped**, 73.180 s in unittest (74.068638 s including process startup). The skipped test requires all 25 actual surface outputs; it will be run after collection. No scientific executable was launched by this suite.

The four modules cover real archived baseline arithmetic and signs, paired coordinates/charges, boundary ownership, gas-phase ORCA parsing, native TABI parsing, actual frozen surface preparation, cache invalidation, preserved implementation compatibility, and explicit malformed-input/failure handling. Corrupt-input fixtures are labelled copies of real artifacts. They are not scientific results.

Exact command, stdout/stderr, return code and tested source hashes: `workspaces/global_electrostatic_20260916/software_tests_v1/`.

Actual calculations are separate: four ORCA vacuum endpoints and four ESP checks completed; two isolated TABI controls completed; the 23 whole-protein surface checks are still running. Their scientific gates are not inferred from software-test success.
