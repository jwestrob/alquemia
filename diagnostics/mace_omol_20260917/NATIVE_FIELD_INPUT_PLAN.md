# Per-site field input: native response replay before a new source model

Declared after the native field-accounting control passes; active-goal autonomy
applies. Six response solves, no new DFT/MACE/charge fit or biological score.

## Question and implementation

Can the native vacuum/GK mutual-response solver accept the already computed
four direct-field arrays and reproduce its original induced dipoles and energy?
This isolates the input interface needed for a later quantum-field model.

Use the three existing ion-excluded frameworks, `frozen_tight` and
`frozen_rigid` variants only. Reuse their exact coordinates, physical parameters,
source masks, GK cavity and four native d/p direct-field arrays from job1201055.
Reuse their actual no-response energies and original mutual energies/dipoles.
No fields are generated, fitted, rescaled, combined or added in this replay.

Extract `induce0c` from pinned Tinker26.2 source commit
87050685eff8840d312e2a332cc82c33f63c7c3d into an isolated derived routine. Only
rename the routine, declare an input array, and replace its `dfield0d` call with
assignment of the four supplied arrays. Preserve every convergence, masking,
preconditioner, mutual-field and GK line. Verify this exact source transformation
before compilation. The maintained library and original source stay unchanged.
Use the same native mutual-field kernels and 1e−7Debye/100iteration settings.
No new PB or induced-dipole solver is invented; this is an input adapter.

The frontend reads an explicit pinned field file, verifies atom ordering and
Born radii, then calls the derived native routine exactly once. It performs no
`energy()` call, force evaluation or atom movement. Export all four induced
arrays, convergence evidence, physical parameters and timing. Reconstruct the
vacuum polarization and solvent total from the verified native contraction and
the compatible reused static component; retain every term separately.

## Frozen checks and scope

- Maximum induced-dipole component difference from original ≤1e−9eÅ.
- Vacuum polarization and solvent total replay differences ≤1e−7kcal/mol.
- Frozen-site induced arrays remain exactly zero; all values finite.
- Same source geometry, permanent/cavity parameters and atom mapping; the
  serialized coordinates, not a separately recomputed rotation, are authoritative.
- Every native solve reaches the unchanged tolerance; nonconvergence stays a
  failure, with no increase in iteration count or fabricated successful output.
- Rigid transformation uses the previously measured/control geometry. Require
  total-energy rigid difference ≤1e−7kcal/mol and rotated induced-vector error
  ≤1e−8eÅ. This is a deterministic replay, not a scientific affinity tolerance.

Exactly six new response solves on64CPU standard/memory with the previous node
exclusions and64GB request. Existing native controls suggest seconds of kernel
work; measure all build, preparation, execution and failed-attempt costs. No
project CPU/time budget. Preserve partial failures and immutable replays.

A pass qualifies this narrow response-input interface only. It supplies no
quantum source reaction field, covalent-boundary rule, metal parameterization,
complete hybrid correction or prediction. Define those in a separate recorded
model before any new metal score; no automatic backend/source substitution.
