# Final delivery verification

- 27 real-artifact pool/comparison/joint/recovery tests pass, zero skips:
  `FINAL_POOL_TESTS.txt`. Uses the existing scientific CPU Python and unittest.
- 9 actual PLM export tests pass, zero skips: `FINAL_PLM_TESTS.txt`. Uses the normal
  Python environment and pytest, as specified in the export commands.
- Released standard1H4I plan still validates with the Nikasha entrypoint:
  `FINAL_RELEASE_CHECK.json`. This check launches no molecular work.
- Independent joint-coordinate preparation checks (7), scaffold inventory checks
  (5) and actual scientific execution receipts are retained in their reports.

The first broad discovery command used the scientific environment's unittest
runner on a pytest module; it recorded one import failure because pytest is not
installed there (`FINAL_TESTS.txt`). No environment was modified. The subsequent
correct runners above passed. This was a test-runner mismatch, not a molecular
calculation or an omitted failing test.

New molecular execution and failures are reported separately in the shared-pool,
angular and joint reports. No new DFT, empirical labels, mock protein or fabricated
successful scientific output was used to pass these tests. All launched jobs are
terminal; production/defaults and old result files remain unchanged.
