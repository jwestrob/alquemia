# Per-site field adapter reproduces native coupled response

Job1201061 completed all six declared response solves. **All nine replay and
rigid-transformation checks pass.** Supplied fields are the actual native arrays
from the preceding controls; no quantum source or biological score was added.

The derived routine differs from pinned Tinker `induce0c` only in its name,
input declaration and replacement of `dfield0d` with assignment of four supplied
arrays. The full mutual/GK solver and convergence code are unchanged. Both the
original library and its source remain intact. An exact inverse transformation
test verifies every other line. No new solver algorithm is introduced.

Induced dipoles reproduce the original native results within
1.8214596497756474e−16eÅ. Reconstructed energies reproduce native energies within
3.092281986027956e−11kcal/mol. Every solve takes16iterations at the unchanged
1e−7Debye tolerance. Frozen-source dipoles are exactly zero. Maximum rigid total
energy change6.99565134709701e−8kcal and vector error4.742733983320591e−15eÅ pass
the frozen1e−7kcal/1e−8eÅ checks. All unrounded terms remain in the result JSON.

## Actual execution and recovery

Four tests pass in7.062s: exact native source transformation, exact input field
serialization, real convergence/response parsing and retention of numerical
outcomes. Two integrations were explicitly skipped before execution. There
were no failed/retried response solves.

The first build failed because Tinker's imported integer `maxval` shadows the
Fortran intrinsic. Replacing the maximum test with the equivalent element-wise
`any` test fixes the wrapper without changing a criterion or native routine.
Both build receipts/source versions are retained. A premature preparation
correctly rejected the unavailable frontend before creating a workspace or
running science; its time is unavailable. The build CLI now propagates failure.

Allocation25wall-seconds×64CPUs=1600allocated core-seconds; actual199.544jobCPU-s,
zeroGPU. Six native kernels total4.325634061wall-s; native process wall/CPU sums
8.835634470/186.129141s. Both builds together2.323078075wall/2.244898CPU-s.
Scheduler batch peakRSS was reported0, which is unavailable sampling evidence;
the actual native child high-water memory and preparation cost are in the JSON.
No new DFT, MACE, force, geometry change or `energy()` evaluation occurred.

## Meaning for the discriminator

We now have a working, inexpensive per-site field input path with checked
native energy accounting. This permits a quantum-density driving field without
forcing it through the failed charge/dipole fits. A complete model still needs
consistent permanent-multipole coupling, a defined reaction-field approximation,
covalent-boundary treatment and partition/predictive checks. Do not present
this interface qualification as increased La/Ca accuracy. Baseline stays default.

Protocol `Tinker26_2_supplied_direct_fields_native_response_replay_v1`;
workspace `workspaces/mace_omol_20260917/native_field_input_v1`, isolated successful
build `native_field_input_software_v2`. No adapter calculation remains live.
