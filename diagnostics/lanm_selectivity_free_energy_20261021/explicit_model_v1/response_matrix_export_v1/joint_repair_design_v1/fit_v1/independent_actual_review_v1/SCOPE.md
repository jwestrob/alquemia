# Terminal fit review

This is a saved-only review of the one declared run_v2 fit. The source
configuration deliberately has `fit_result: null`. Root may copy it to a new
runtime configuration and bind the terminal RESULT pin; its path is fixed.
An optional pinned scheduler-accounting JSON may be attached. The script
refuses `started` and any result without a terminal cost record before opening
vectors, the ledger, or cached native artifacts. No live fit artifact was read
during preparation.

The reviewer reconstructs all clipped derivative requests, the saved Jacobian,
relative-cutoff Gauss–Newton direction by SVD, trust-region scaling, recorded
BVLS gradient/KKT diagnostics, and four ordered trial requests. It verifies
strict improvement with current-point tie retention. A derivative record can
be selected only when its exact hash also belongs to the permitted trial batch.
No BVLS call, parameter optimization, residual evaluation at a new point, or
native call is made.

All complete vector scores are replayed with the frozen scorer. An independent
response normalization reconstruction checks the four context-weighted pair
blocks, four C vectors, four fixed-reference susceptibility states and all 531
weights. The original 413-vector E/F function and its 25 gates remain unchanged.
The initial model/response/score must match saved zero-call preparation; new
runtime cache identities may additionally bind the qualified binary.

All vector cell identities, exact permitted key edits, geometry/parameter bytes,
program argv/output pins, finite counts, failures and cache reuse are checked.
The selected candidate's new native cells additionally undergo the original
saved raw parser; its inherited cells must exactly match the admitted baseline.
Other candidate raw outputs receive pin/receipt checks and score replay, not a
second complete raw numerical admission. Full selected raw admission is not
fresh validation and does not add a physical response tolerance.

The reviewer did not author the fit driver or selection algorithm. It did author
the evaluator and objective. Reusing their pure parsers/scorer is disclosed;
those routines are not presented as an independent scientific oracle. The
independent source reviews, actual native qualification, exact reference targets
and source pins remain part of the chain.

Driver child CPU is retained. Per-program process-global child-CPU differences
overlap across concurrent workers and are not summed. Native elapsed times are
summed separately. Exact scheduler allocation cost is unavailable until root
binds an accounting artifact; it is never inferred from Python wall time.
The source preserves all failed cells/vectors and does not admit a failed fit as
a frozen candidate. A passed review of terminal bookkeeping is distinct from
passing the original 25 E/F gates.

Preparation: 15 focused controls use the already admitted initial point and
small selection/guard examples. No fit result, partial vector, native output
from the live run, fresh validation, or trimer regression was opened. Root owns
any eventual review execution and terminal result binding.
