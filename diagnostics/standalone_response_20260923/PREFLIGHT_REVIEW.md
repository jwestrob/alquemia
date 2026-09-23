# Prelaunch integration review

Root independently reviewed the evaluator's composite signs/units, cap Jacobian,
common physical constraints, actual start-force reuse, two-by-two eight-thread
layout, limit-versus-failure admission and common cross-metal pool. Three technical
fixes requested before launch: null-safe reports, per-case matrix initialization,
and explicit SciPy wrapper/kernel/version pins. All are in immutable run_v3.
run_v2 remains preparation-only; it was never executed. Scientific settings are
unchanged from PLAN_v2.md.

Actual run_v3 manifest SHA256:
`aa8fcc7058aa89fceb467820f147973af6fdd5d7fec5207e18a2b36d13890b85`.
Six source rows/twelve searches;36real native force cells checked;56reused standalone
baseline cells and16fresh baseline calls. Maximum1000standalone/480MACE calls.

Five actual-fixture tests pass; the sixth is explicitly skipped until molecular
execution. A prior test field-name mistake was corrected using the actual stored
analytic derivative key; no molecular output or calculation was fabricated.
The Q88 source/cap projection reproduces its existing qualified standalone derivative.

SLSQP stopping on small energy changes is not stationarity evidence. Initial/final
native, solvent and total projected gradients and real evaluation counts are retained.
No elapsed/allocated CPU budget stops the project: the finite count is the declared
proposal algorithm and comparison scope. No general backend or production change.
