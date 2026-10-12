# Guarded-dispatch code review

Reviewed the separate research dispatcher at
`workspaces/orca_parameter_audit_20261012/worktree_runtime/scripts/guarded_native_gfn2_dispatch.py`,
SHA256 `581af16b3184e883a2d46c3c0dea9396c28016e304086e41b6fab23f9c7ae7ff`.
This was a read-only code review; no molecular endpoint was executed by the
reviewer.

The dispatcher stages the exact approved full parameter file, runs shared
preflight, and requires parameter postflight plus separate process/state/SCF/
population admission. Unsupported state/profile/runtime combinations have no
default fallback. A dry run remains explicitly unscored. Task directories are
fresh, numerical children are serial, and prior results are not overwritten.

The first revision did not install the source executor's signal handlers.
SIGTERM/SIGHUP could therefore terminate the wrapper while its separately
started child process survived. The reviewed repair now installs SIGINT,
SIGTERM and SIGHUP handlers that set the owned executor's stop event, allowing
its existing process-group cleanup to run. Prior handlers are restored.
Post-execution errors also retain a failure receipt and artifact hashes through
the dispatcher's `finally` block.

No further blocking defect was found in this focused review. Real-input
dry-run and rejection evidence belongs to the runtime agent's separate
`runtime/GUARDED_DISPATCH*` records; those records must bind this final source
hash. The dispatcher's execution branch has not itself run a molecular
endpoint. The 60 earlier controls establish runtime/parameter evidence through
the original executor and must not be described as guard-dispatched runs.

This review does not make a serial, explicitly bounded research wrapper a
replacement for every historical production or MPI caller. New callers must
use the guard contract and retain their own scientific admission, resources,
job ownership and failure records.
