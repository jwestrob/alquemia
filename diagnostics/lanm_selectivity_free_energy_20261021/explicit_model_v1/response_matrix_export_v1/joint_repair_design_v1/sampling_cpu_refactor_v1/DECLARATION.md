# Vacuum CPU common-field and force-contraction refactor

Declared 22 October 2026 before source generation/build/native execution.

Question: can the already qualified regularized local-response functional use
common direct/mutual-field dispatch and O(Q) derivative contraction while retaining
the same vacuum energy, dual response, conservative Cartesian/frame forces and
virial? Q is the number of selected metal–water pairs. This is an implementation
refactor, not a new interaction family or periodic extension.

Exact authorities are in SOURCE_INPUTS.json. Copy the qualified original water
source manifest82886b76 into an isolated source_v1; do not change the current
fit/evaluator/native tree or the separate donor prototype. Preserve alpha,
pdamp/POLPAIR, coefficients/domain, radial/contact functions, local-frame
routines, all energy routes, both RHS, original PCG/peek and readback format.

Move h and −H*x from the vacuum kernels into common field dispatch used by the
original solver and diagnostics, exactly once. Retain original branch order and
external-field ordering. Replace the dense lrdh coordinate Jacobian by stored
per-water geometry and two-pass contraction of the same stationary dual gradient;
include all other-water derivatives, metal reaction and native symmetric virial.
Use prefix/suffix environment sums to avoid subtracting a dominant water from
the total. Floating-point summation order may differ for nonzero coefficients;
the Hamiltonian and prospective qualification thresholds do not change.

Keep enabled periodic/replica/Ewald/list/host/predictor/alternative-solver guards
and the1–28 active-site/single-metal scope. Preserve the existing dense normalized
SPD certificate and exact absent/all-zero bypass. Common dispatch alone does not
qualify or enable any previously rejected route. No GPU, MD, new geometry,
parameter fit, reference call or production change is included.

Planned next native work, only after independent source review and root release:
reuse the pinned BUILD.py/BUILD.sbatch to compile the isolated manifest, then
reuse exactly the original nine-fixture57-program/636-internal-FD qualification
inventory and unchanged numerical gates in a fresh run. All57 must run against
the rebuilt executable; none of the old executable's27 recovery-prefix rows may
be reused. Old baseline authorities remain comparisons. Preserve native failures
and costs. Build/executable and execution-review fields remain unbound now.

Source checks may compare the old dense and new contracted formula on a small
fixed synthetic geometry set and check dispatch/guard text. They invoke no native
program or optimizer and do not replace actual qualification. Root owns builds,
submissions and eventual numerical admission. One numerical thread for preparation.
