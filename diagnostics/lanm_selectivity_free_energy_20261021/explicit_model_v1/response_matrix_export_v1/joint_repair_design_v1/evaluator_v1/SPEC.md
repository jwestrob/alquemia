# Finite native joint evaluator

Source preparation only until the exact new native build and independent
energy/force/moment/export qualification are bound. Root owns execution.

`Evaluator(config, cache_dir, root_release=None)` exposes `initial_vector`,
`family`, `evaluate(values)`, `snapshot()`, and `stop_event`. Evaluation returns
`model`, `response`, `score`, and all `cells`. The score delegates the frozen
joint objective and unchanged original 413-residual E/F objective; it retains
the original 25 acceptance gates separately from the engineering objective.

## Exact finite scope

Use the declared 12 box variables. Per represented metal, map the four response
variables to physical a,b,c,d using the pinned domain; change only that metal's
water-O VDWPR R and epsilon, and add one POLAR-LOCAL-RESPONSE record. All other
active baseline key lines and parameter bytes, including native POLPAIR, remain
fixed. Add missing reporting flags explicitly. Each fixture keeps its original
selected17676 geometry, graph, atom order, native names and multipole frame.

The objective requires 46 unique metal fixtures: 10 original E/F cells,
34 literal extracted dimers and two W9 far trimers. Analyze EM runs once for
each changed fixture; analytic testgrad runs once for each of the 10 E/F cells.
Only the four susceptibility aliases export matrices. Analyze and testgrad use
separate directories so their outputs cannot overwrite each other. There are
56 programs per cold full vector and 28 per changed metal. No finite-difference
programs run within this evaluator. The six other trimers remain postselection
regression, outside this objective. Exact cached 46 metal-free terms and two
bare ions complete the moment cycles. No reserved outcomes are opened.

Initial response amplitudes are zero and R/epsilon equal selected17676.
The admitted 54-fixture baseline and its actual independent review supply
energies, moments and four analytic matrix susceptibilities; the original
selected17676 record supplies all 10 forces and true native producers. This
initial point needs zero new calls. The baseline moment frames and cached
metal-free values are reused from the pinned saved frame preparation. New
native zero/off equivalence is a required execution qualification, not a claim
that the appended reporting/zero records are byte-identical to old inputs.

## Native admission and cache

Invoke the pinned existing native runner and parsers. Analyze admission precedes
testgrad. Preserve strict energy/component/induction, full force, full matrix
SPD/conditioned-solve and native stopping checks. Active response rows additionally
retain exact physical coefficient, metal/O, gamma/distance, environment, spectrum,
full induced u/v and dual-channel readbacks. Full-precision conjugate dipoles use
the mean of u and v plus the unchanged permanent moment; the old .uind and native
printed total moment remain independent checks. Zero response retains the old
native output route. The four susceptibilities use the z projection of the
qualified A+H inverse, with the original native-to-atomic-unit conversion.

Cache identity includes fixture, literal geometry/parameters, rendered key,
effective values for its represented metal, source, binaries and parser pins.
Changing only the other metal reuses the same physical cell and preserves its
original producer. Disk reuse rechecks exact input bytes, argv, receipts and raw
admission. Failed or incomplete cells are retained and never retried. Every
requested vector retains all 46 slots, including stopped/unattempted slots.

One 12-CPU Slurm task runs at most 12 single-thread native processes. The hard
ceilings are 673 distinct vectors including the initial point and 24,000 native
program invocations, including failures. The count is persisted before each
call. The first numerical/native failure sets the shared stop event and blocks
new invocations; already running calls retain their results. Child process CPU
receipts use the existing process-wide counter and may overlap across workers;
scheduler allocation cost is the authoritative total runtime cost.

## Runtime binding

The scientific evaluator configuration is unchanged by the driver's explicit
root release argument. Its fingerprint excludes only `root_release` and
`preparation_result`. Runtime configuration requires exact source manifest,
completed build, analyze/testgrad executables, and actual independent native
qualification. Root release status is `released_bounded_joint_response_fit`,
with the driver's exact 24/672/24000/12 limits and additional evaluator/source/
manifest/domain/qualification/cache-directory bindings. Missing bindings fail
before any native call. The source-only template leaves these bindings null.

Controls are limited to saved initial replay, exact rendering, effective-metal
cache identity, finite box and hard gate rejection, and the readback interface.
No optimization, new candidate molecular evaluation or physical accuracy claim
is part of preparation.
