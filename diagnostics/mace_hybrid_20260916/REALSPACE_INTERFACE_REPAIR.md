# Isolated realspace interface repair

The first standard-GPU attempt, job1200302, failed before producing an energy.
47-atom La core; 23 allocated seconds,16 CPUs,one RTX A5000,64474 MiB host RAM.
The unmodified MACE0.3.16 forward constructs a reciprocal grid unconditionally.
Its nonperiodic neighborhood code assigns a large virtual cell based on absolute
coordinates; constructing the unused grid exhausted GPU memory. This failure
does not measure the actual realspace core calculation's memory requirement.

Further installed-source inspection and checkpoint loading, without inference,
identified interface incompatibilities with graph-longrange0.4.4:

- The calculator accepts but does not apply `pbc_handling`.
- The checkpoint predates backend dispatch attributes (`pbc_handling`, etc.).
- MACE passes obsolete `force_pbc_evaluator` and `pbc` arguments.
- MACE sends feature densities as `[N,1,4]`; this backend expects `[N,4]`.
- The checkpoint lacks the new `features_dim` metadata on four self-interaction
  blocks. It is exactly the saved `features_irreps.dim`; restored without changing
  any saved coefficient, buffer or learned parameter.

The versioned adapter `scripts/mace_realspace_compat.py` uses the backend's
public realspace dispatch setters, bridges those signatures and the singleton
channel, and skips the unused reciprocal grid. It rejects any periodic input
or request. Both maintained realspace paths ignore reciprocal tensors/volume:
`GTOElectrostaticFeatures._precompute_geometry_realspace` stores only positions
and batch; `_forward_dynamic_realspace` delegates to the original realspace
feature kernel; `GTOElectrostaticEnergy._forward_realspace` delegates to the
original realspace energy kernel. No kernel, weights, precision, coordinates,
charge, spin, solvent, cutoff or scientific task is changed.

The original isolated environment and all failed products remain immutable.
`pilot_v2` repaired the interface/grid but failed on missing `features_dim`;
job1200306 used11 allocated seconds. `pilot_v3` includes the metadata restoration
and all four cores completed in job1200308. It reuses every exact original input
byte via `repair-interface`, with fresh implementation/cache hashes and the same
scientific protocol and12 calls. Independently preparing v2 had changed only
binary roundoff in the rotated full-system files; none of those tasks ran.
The final v3 reuses the original files directly, avoiding that issue. This is
technical recovery under AGREEMENT.md, not a new analytical variant. Held then
cancelled our own still-pending H200 job1200207 and stopped only its task-owned
continuation2628149/child2628154 to replace its broken implementation. No
inference was performed by1200207. Other agents' jobs/watchers were untouched.

Official documentation recommends these package versions and realspace mode:
[MACE-POLAR usage](https://mace-docs.readthedocs.io/en/latest/guide/polar_mace.html).
Actual installed sources are pinned by the software inventory; the incompatibility
was established from those sources and the actual checkpoint, not assumed away
because the versions match the documentation. No package was upgraded.

Validation: all checkpoint state-dictionary parameters/buffers are bitwise
unchanged. On all four real computed core density/coordinate fixtures, adapter
feature values, Coulomb energies and coordinate gradients equal the original
realspace kernels exactly on CPU. Periodic input/forced-periodic requests are
rejected. These are interface/kernel checks; no original end-to-end result can
be compared because that original software combination did not execute.
