# Native GFN2 iteration-limit diagnostic — frozen before execution

Protocol: `adaptive_pool_native_GFN2_MaxIter500_diagnostic_v1`.

Parent authorized this exact contained diagnostic on 2026-09-22 under Jacob's
standing development permission. It asks whether the single adaptive-pool SCF
failure reaches the unchanged convergence criteria with more iterations, while
two successful, exact origin controls reproduce the primary energies.

## Three fixed calculations

1. 4MAE La at the adaptive-Ca candidate, vacuum: original pool task failed after
   125 cycles. Retain that failure and all its output unchanged.
2. 1H4I La at origin, vacuum: exact successful cell linked by pilot_v2 matrix.
3. 4MAE La at origin, vacuum: exact successful cell linked by pilot_v2 matrix.

Copy each source XYZ byte-for-byte. Add only `MaxIter 500` to its existing `%scf`
block. Preserve native GFN2 parameters, 300 K electronic smearing, special xTB
mixer, normal convergence tolerances, NoAutostart, SAD fresh guess, charge and
singlet multiplicity, and all molecular coordinates. No orbital seed, geometry
search, different solver/backend, DFT or subsequent iteration-limit escalation.

## Acceptance fixed before calls

Both controls must terminate normally and reproduce their exact original vacuum
energy within 0.05 kcal/mol per cell. This engineering tolerance allocates the
existing 0.20 kcal/mol four-term score tolerance across four possible errors; it
does not measure the unknown error at the formerly failed geometry. Actual
parameter exports, state, tolerances and driver settings must agree, except for
the explicitly increased iteration limit. Failed-point convergence is separate
from control agreement. Any unavailable calculation remains unavailable. No
scientific or predictive qualification follows from iteration-limit agreement.

The primary adaptive collection is not modified. Parent decides whether a
separately named recovery overlay is justified after inspecting the diagnostic.

## Execution and deliverables

Existing pinned ORCA 6.1.1 runner; three concurrent eight-rank single points,
24 allocated CPUs, 64 GiB host memory, no GPU. Use an explicitly selected
available GPU-node CPU allocation or available standard allocation. Three finite
tasks, no artificial aggregate compute/time ceiling. Record actual job, receipts,
SCF cycle counts, endpoint energies, state/parameter checks and measured costs.

Initial trace: intended native LEANSCF/xTB mixer and SAD guess are active. The
failed 4MAE task reaches cycle125 with Delta-E=-4.96e-5 Eh, above TolE=1e-6 Eh,
after substantial earlier oscillation. Its last iterate is not an endpoint.
