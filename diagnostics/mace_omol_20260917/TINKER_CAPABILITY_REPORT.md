# Tinker capability result — 2026-09-17

**A pinned CPU backend and real framework adapter now work. No scientific
energy or force has yet run.** This is preparation for a coupled polarization
test, not an improved classification result. Baseline/default unchanged.

## Executed checks

Built Tinker 26.2 (CPU commit
`87050685eff8840d312e2a332cc82c33f63c7c3d`) in the isolated workspace. GNU
Fortran11.4, `-O2 -g -fno-fast-math`, OpenMP; no working environment modified.
Build1201010 found an upstream CMake omission: `mdsave.f` uses `uatom`, but the
library source list omitted existing `uatom.f`. A deferred `target_sources`
include adds that unchanged upstream file; recovery1201011 succeeds. Source
checkout is pristine. Failed logs/manifest/cache and both receipts survive.

Exported and natively initialized exactly three real protein/water frameworks:

| Case | Atoms | Frozen source atoms in framework | Native keyword reenabled | Adapter mask restored |
|---|---:|---:|---:|---|
| GGR1GLG, extended-core support |4697|57|57|Yes|
| Alpha1F6S |1931|51|51|Yes|
| Alpha6IP9 |1897|54|54|Yes|

The selected metal remains explicitly unparameterized in each physical ledger;
these are not full-system scoring calculations. The masks come from already
recorded QM/cap projection support. They do not invent a new geometry or charge
distribution. Source coordinates, connectivity, monopole charges and multipole
axis mappings match exactly. Polarizability/damping differences are at most
2.22e−16 in Å³/Å^(1/2); the largest local-multipole input difference after unit
conversion is 1.071e−10 in OpenMM MD units. All differences are retained; no
energy-accuracy claim follows from this format comparison.

The native `kpolar` parser reads POLARIZABLE, then later sets `douind=true`
again for nonzero-polarizability sites. The actual outputs confirm this on
every requested frozen framework source. Our parameter adapter applies the
recorded mask after `mechanic` and verifies that polarity and damping remain
unchanged. **Whether the subsequent coupled SCF obeys this mask remains
untested.** The parameter-only executable calls no energy/force routines.

Four tests pass in2.501s, using executed outputs and corrupted copies: real
round-trip/source masks, missing completion, reenabled frozen flag and zeroed
damping. Tests make no new native calculations. CLI dry-run passes.

## Backend and primary-parameter findings

CPU source implements POLPAIR in its pair-damping matrix and uses `douind`
through both vacuum and GK induced-dipole solvers. This supports the next
numerical check, not automatic qualification of our proposed hybrid.

The GPU repository now redirects from `tinker9` to `tinker-gpu`. Inspected
commit `44bcd7c898e68f6827c112df79b8361e443730d2` pins a different Fortran
submodule (`33dcaf084a3e85730559d4e3a3188d9dd87a0f10`). Its actual energy
dispatcher has no GK/solvation calculation; the SOLV enum alone is insufficient.
Its polarization data copy transfers polarity/pdamp, not `douind`. No GPU build
was attempted. Do not claim that the CPU mask fix supplies GPU source freezing
or coupled GK. CPU remains a possible affordable route; solver cost is unknown.

Both primary supplemental files were retrieved from their Figshare records
and pinned (TXT SHA `eaa7da82c749c482d2c0b5110aa3cf1f1098e59600b016ed8e155471a32e06cd`,
PDF SHA `5f112c7aa4976e1c365d77711678b762ab672c368e2cf6e6b0b6ae8e118f2100`).
The [TXT supplement](https://doi.org/10.1021/acs.jpcb.2c07237.s001) contains
actual La09/La22 input blocks and real complex geometries. Its La09 amide
overrides differ from the article's rounded Table1, so their records must
remain distinct. The [PDF supplement](https://doi.org/10.1021/acs.jpcb.2c07237.s002)
also reports a roughly49kcal/mol monodentate-acetate error. The published model
therefore does not establish uniformly accurate donor chemistry. No supplement
parameters have yet been applied to a scored protein. Ca type358 is present
in the pinned AMOEBA2018 file; cross-version La/Ca compatibility is unverified.

Sources: [CPU code](https://github.com/TinkerTools/tinker/tree/87050685eff8840d312e2a332cc82c33f63c7c3d),
[GPU code](https://github.com/TinkerTools/tinker-gpu/tree/44bcd7c898e68f6827c112df79b8361e443730d2).

## Cost, artifacts and next action

Builds used **2560 allocated core-seconds**, **285.312 actual CPU-seconds**,
**zero GPU time**, including the failed build. Export preparation took4.314806s
wall/4.063931s CPU; the three native parameter reads took3.220401s wall and
2.049645s CPU. Peak native child RSS378356KiB; export-process peak483332KiB.
These are engineering costs, not measured production scoring costs.

Workspace roots under `workspaces/mace_omol_20260917/`:
`tinker_sources_v1`, `tinker_software_v1`, `tinker_capability_v1`.
[TINKER_CAPABILITY_RESULT.json](TINKER_CAPABILITY_RESULT.json) carries exact
hashes, receipts and differences; [commands](TINKER_CAPABILITY_COMMANDS.md)
reproduce inspection without rerunning preparation. No class/reference or
mechanical correction is available from this stage.

Next is the declared CPU framework solver check in
[TINKER_FRAMEWORK_SOLVER_PLAN.md](TINKER_FRAMEWORK_SOLVER_PLAN.md), before
spending new DFT work or treating this backend as a scorer. Retain the baseline.
