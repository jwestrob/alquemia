# Term collector source review

Passed; no source correction requested. This is a source and saved-control review,
not execution of the collector or a numerical qualification result.

For each of28 terms, the analytic strain slope is the after-minus-before cumulative
virial increment contracted with the frozen yy strain matrix/L0. The matching
energy is the individual term energy, so its central difference is
(E_plus−E_minus)/(2h). Reported discrepancies are FD−analytic; two-step change is
half-step minus full-step. These signs match the inherited positive dE/dq
convention. L0=1Å preserves the stated kcal/mol/Å versus dimensionless-strain units.

INITIAL virial is zero; the observer checks all nine-component telescoping sums,
zero post-term virial, unchanged component energies at FINAL, and arithmetic/print
energy-sum closure. The collector independently requires exact central projection
closure and bounded component-FD sum closure for both steps. Bounds apply to
arithmetic bookkeeping, not enlarged physical acceptance thresholds.

Exactly five states are bound: center and yy-strain ±1e−4/±5e−5Å. Each must match
its historical total energy, full gradient and total virial token digest exactly,
with separate polarization/constants/state/role checks, native source/build/input
and one-program receipt joins. The same frozen helper computes the two original
FD checks and original two-step check at0.001 kcal/mol/Å plus inherited printing
allowances. All three result dictionaries must reproduce the original failed
checks exactly. Collection completion explicitly preserves that numerical failure.
Per-term slopes are descriptive; no new per-term gate or active-response/physical
accuracy admission is introduced.

The twelve saved manufactured controls were read, not rerun. They exercise signed
term cancellation,28-term completeness, closures, missing/corrupt endpoints and
five-parent scheduler accounting. Scheduler cost includes each of the five parent
allocations once, with build and historical baseline costs separate. Actual output
bindings and root release remain required before saved collection.
