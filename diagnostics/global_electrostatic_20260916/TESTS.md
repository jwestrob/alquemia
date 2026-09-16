# Executed checks — final

**40 tests passed; zero failed or skipped.** All 25 actual solver tasks were
available, so the previously skipped scientific-output integration test ran.
Unittest time: 115.233 s. No scientific executable was launched by the test suite.
The four quantum endpoints, four ESP utilities and 25 TABI tasks are separately
reported scientific calculations, with their own receipts and failures.

Modules: `test_global_electrostatic_preparation.py`,
`test_global_electrostatic_accounting.py`, `test_affordable_tabi.py`,
`test_global_electrostatic_assess.py`, `test_global_tabi_controls.py`.

Tests use real archived/current coordinates, energies, native control echoes,
meshes and execution receipts. Explicitly corrupted copies test malformed-input
handling. Tests cover charge/coordinate ownership, score algebra, cache changes,
missing/duplicate/partial artifacts, native controls, actual component accounting
and source compatibility. A valid scientific rejection passes software tests;
the suite does not require the challenger to win.

Exact command, output, return code, process timing and tested source hashes:
`workspaces/global_electrostatic_20260916/software_tests_final_v1/`.

Historical checkpoint `software_tests_v1/`: 37 tests, 36 passed and one skipped
while actual solver work was pending. Three later native-control tests passed
separately and are included in the final 40-test suite. Prior receipts remain.
