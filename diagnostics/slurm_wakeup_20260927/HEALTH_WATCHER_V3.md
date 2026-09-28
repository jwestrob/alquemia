# SCF health v3: near-root stagnation is a review event

Version `native_scf_health_v3` adds the predeclared root-requested rule: the latest12consecutiveTRAH/NRmacro errors all exceed5times the actual printed TRAHgradient tolerance, improve by less than twofold from first to last, and have total energy span below1e−6Eh. The error threshold is parsed from `Converg. threshold (grad. norm)`; if absent, this rule is unavailable rather than guessing a value. The NRmacro rows are part of the same solver trajectory and must not be skipped.

Alert name: `near_root_stagnation_review_only`. It asks root to examine branch convergence; it never cancels, loosens acceptance, accepts unconverged energies, or imposes elapsed-time/CPU caps. Stable energy alone does not establish convergence. Version1/2policy records remain unchanged; armv3on a fresh receipt.

Eight actual-fixture tests pass. Both retained smallDy HCore/PModel logs trigger from their real molecular prefixes; HCore's window81–92 includesTRAH andNR. Completed realLa outputs do not trigger. Removing the tolerance line in an explicitly corrupted real-log copy disables the new rule. Tests retain all earlier gross-DIIS/TRAH/success/progress checks. Exact timestamps, hashes and parsed windows: HEALTH_POLICY_V3_AND_TESTS.json.

No new chemistry, jobs or queue messages were launched by this change. This policy was developed from consumed failures and is a prospective operational alert for the next capability case, not independent scientific validation.
