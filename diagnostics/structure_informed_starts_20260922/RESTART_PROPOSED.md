# Native charge/multipole restart continuity — proposed, not executed

The completed exact-input repeats reproduce the coordinate-sensitive Q88 La
vacuum energies. A bounded next diagnostic can test whether the two electronic
solutions persist when each geometry is initialized from the other's solution.
No job/manifest for this follow-on has been launched by this proposal.

## Actual available route

The [ORCA 6.1 native-xTB documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb)
describes a basename.xtbw restart containing atomic charges and multipole moments;
its presence takes precedence over a GBW restart. This is the native mechanism,
not an assumption that a generic MORead keyword activates the same path.

[RESTART_SOURCES.json](RESTART_SOURCES.json) pins the four actual saved solutions
from repeat 1 of job 1210115: old/new geometry × vacuum/ALPB. All have successful
native-mixer outputs, same atom order, charge −2, singlet, 438 electrons, matched native
GFN2 parameters and retained 11,864-byte endpoint.runtime.xtbw files. GBW files also
exist but are unnecessary for the proposed native test. The geometries differ
by at most 6.4623e−7 Å. Normal termination does not establish fully converged native
charges; retain their printed density diagnostics as measured.

## Proposed eight calls

For each exact destination XYZ (old, new) and each medium (vacuum, ALPB), run:

1. A same-geometry/same-medium xtbw seed control.
2. The opposite-geometry/same-medium xtbw seed.

Total: 8 native GFN2 La singlepoints, comprising 4 self controls and 4 cross starts.
No Ca, MACE, DFT, geometry search, new chemical state, changed temperature or
parameter. One 64 CPU / 128 GiB allocation, 8 × 8 MPI, common host. Keep the current
primary native mixer and MaxIter 500. Every new working directory receives only
its exact XYZ/input and the recorded
seed copied to **endpoint.runtime.xtbw**, matching the existing runtime renderer's
actual basename. Preserve an immutable separate seed copy/hash before ORCA can
rewrite the active file. No GBW, parameter-override file or other restart is copied.
Remove NoAutostart consistently from all 8 inputs to request restart; it changes
the initial-guess policy, not the Hamiltonian or physical state.

Before interpreting any energy, verify actual native-mixer execution and evidence
that the supplied charges/multipoles were read. The manual establishes the
supported mechanism; local consumption/diagnostic printout has **not yet been
tested**. If no positive restart evidence can be obtained, report
restart_not_confirmed rather than calling changed energies a successful restart.
The existing native output parser checks must not silently substitute an ordinary
SCF run. Retain full iteration traces, occupations, charges, final energies,
printed residuals and all seed-before/after/receipt hashes.

Compare all four initialized solutions per medium as a matrix. Determine whether
cross starts converge to one common energy, retain two apparently stationary branches,
or remain sensitive to residuals. Report both vacuum and ALPB and matched transfers;
do not choose the branch that yields a preferred class or overwrite the earlier
scores. An energy-based electronic-state policy would require its own consistent
model/reference and adequate convergence evidence. This diagnostic alone cannot
establish which native solution is the appropriate ground state or supply a newly
validated classifier.
