# Technical restart correction — 2026-09-20

Job1203265 completed24/24 native tight endpoints; all six transfer/pair checks
pass unchanged tolerances. The largest score change is0.025544882kcal/mol.
The eight proposed `.xtbw` ordinary-SCF restarts failed. Actual output starts
from SAD and reproduces the unseeded first-iteration energy exactly; the ordinary
solver did not consume the charge/multipole restart. This is an unsupported
restart mechanism, not evidence that the saved native solution is unstable.

Correct the same eight solver checks using the actual archived native `.gbw`
orbitals with explicit `!MORead` and `%moinp "seed_source.gbw"`. These files exist;
their source calculations, coordinates, charge, parameters and energies are
already pinned. Copy them to fresh directories, retain the source hashes, and
do not place an `.xtbw` file there. Verify actual MORead/orbital-restart output.
The new preparation stage is `orbital_restart_fix`; original32-call manifest
and all failed outputs remain immutable.

This adds exactly eight low-level attempts to repair the agreed restart test.
The Hamiltonian, tolerances, cases, electronic temperature, nuclear coordinates,
water inventory and scientific question are unchanged. No new classification,
calibration, model or default change. Same existing runner and64CPU/128GiB
allocation, no project timeout/compute budget. Actual failure cost is retained.

Explicit ORCA restart syntax was checked in the official manual:
[SCF initial guess/restart](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html#restarting-scf-calculations).
This technical correction is within Jacob's "pursue!" and the declared
same-state solver-restart qualification; root informed him before execution.
