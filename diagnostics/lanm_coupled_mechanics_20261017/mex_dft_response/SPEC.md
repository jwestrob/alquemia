# Mex embedded reference: fixed eight-cell response inventory

**Source-only preparation. No DFT endpoint is submitted or authorized to start by
this inventory.** Actual derivative admission and root's measured cost review
remain prerequisites. Root owns all execution and subsequent finite decisions.

## Frozen selection

Use the four named geometries already admitted by common7 job 12981, regardless
of their energy rankings. Each is scored as both La and Dy:

| Geometry alias | Definition | Existing search status |
| --- | --- | --- |
| `La_terminal` | Old own-La restricted endpoint; exact broad-La start | 12 accepted continuation steps; finite proposal, no stationary-minimum claim |
| `Dy_terminal` | Old own-Dy restricted endpoint; exact broad-Dy start | 12 accepted continuation steps; finite proposal, no stationary-minimum claim |
| `La_cartesian_terminal` | Module2 broad own-La terminal | 80 accepted steps; not converged |
| `Dy_cartesian_terminal` | Module2 broad own-Dy terminal | 80 accepted steps; not converged |

`INVENTORY.json` binds all eight actual paired POINTs, physical coordinates,
current core/field/cap maps, MM ledgers and global solvent artifacts. It also
retains the actual search pins and previously reported full-force residuals.
There is no reselection, added candidate or replacement of an unavailable cell.
The original PBE0 source coordinates are not the restricted terminal geometries.

## Hamiltonian and physical state

The reference is

`PBE0_embedded(Q209;1840 shifted point charges) + saved MM ledger + saved fixed-FF ALPB(R2049)`.

PBE0 replaces **both** local MACE vacuum and the GFN2 embedding difference.
Reuse the existing eight MM and ALPB components after native provenance/state
replay. No new MM, ALPB, GFN or MACE call is required by this matrix. The MM ledger
is the same pinned `f37fb9d6...` policy at all four geometries. The global solvent
charge vector is the same unshifted `12b03bb0...` vector; it is distinct from the
shifted electronic embedding field. Preserve actual ledger terms and exclusions,
including its explicit ion treatment, rather than interpreting legacy descriptive
spectator text as a different force field.

Keep all 2,049 physical atoms, 167 retained waters, four synthetic caps, EF3 target
and Nd(III) spectators at EF1/EF2/EF4. Spectator identities/charges are unchanged;
their coordinates are those of each saved terminal. Whole charge is +2 and QM
charge is 0. This remains conditional EF3 exchange with three Nd spectators.
It is not pure-metal binding affinity or an equilibrated free-energy calculation.

Reuse the archived original PBE0-D4 / def2-TZVP / def2/J / RIJCOSX / DefGrid3 /
VeryTightSCF / PModel / NoAutostart deck, `DoEQ false`, custom basis/ECP/auxiliary
bytes and ORCA6.1.1. Preserve physical La multiplicity1 and Dy multiplicity6,
ECP46/55, effective singlets and 832 explicit electrons in the spin-free frozen-f
representation. No SOC or independently stable electronic state is implied.

## Exact-cache audit

The inventory audits pinned native manifests in
`workspaces/metal_environment_response_20260926` and the actual eight-state Mex
derivative preparation. It compares same-metal core symbols/coordinates and the
complete ordered external charge/coordinate records exactly, allowing only text
format differences that parse to identical values. No distance tolerance or
nearby-coordinate reuse is permitted. Core/environment mismatch rejects reuse
before any method comparison; an exact match would additionally need complete
method, basis/ECP, state, physical map, ledger, solvent and native-admission checks.

No audited PBE0 task matches both core and environment for any of the eight
selected states. Thus the current finite ceiling is eight new embedded DFT
evaluations. This conclusion is bounded to the pinned archive inventory; it does
not assert that no unregistered calculation exists elsewhere.

## Execution proposal and gates

Implement the thin adapter described in `ADAPTER_PROPOSAL.md` using the existing
native manifest runner, renderer, admission and force-scatter helpers. Freeze its
executable and numerical comparison declaration before new response outcomes.
Retain `EnGrad` so each of the eight finite calls supplies energy and analytic
QM/point-charge derivatives; do not add separate central-gradient calculations.

First proposed numerical/cost scout is `La_terminal__La`, the first cell in the
fixed ordering, independent of signs or gains. A full response matrix requires
the scout plus the remaining seven cells; a failed endpoint stays unavailable.
No automatic retry, optimizer, alternate initialization, different step size,
electronic state, convergence tolerance or additional region is included.

All new response calls remain blocked until the actual 12977/12978 La external-
charge pair is admitted, passes its frozen derivative gate, and root records
measured tractability. The remaining six derivative qualifications retain their
own prior gate and cannot be inferred from approximate-model checks. The final
scope of the two-direction/two-metal qualification must remain explicit in any
response interpretation; a partial pass is not a global force qualification.

Use 12 actual MPI ranks and allocation-derived memory under current Slurm policy,
without memory/time overrides or use of another user's allocation. Record all
attempted/failed/accepted cells, SCF observations, residual flags and actual costs.
The eight-cell denominator and existing nonconvergence survive every result.

## Comparison proposal to freeze before execution

Report per-metal work from each old restricted geometry to each broad geometry,
the own-metal-path differential, and source-compatible shared four-geometry
paired pools. Compare PBE0 electronic work with the replaced MACE+GFN field work;
show unchanged MM/ALPB contributions separately. Retain the plan's 1 kcal/mol
response resolution. Do not compare unlike absolute total energies.

Use both already-declared own-metal old-to-broad physical displacement chords as
force-projection directions. Normalize over all 2,049 physical Cartesian atoms
without picking high-force atoms or directions after DFT results. Generate and
pin the vectors in a zero-call Slurm preparation. Projection magnitudes and
errors remain descriptive unless their numerical resolution is independently
established; the existing two tiny-displacement derivative tests do not certify
the entire broad-motion chord or all gradient components. Retain raw gradients,
projection definitions, printed-gradient precision and SCF uncertainty flags.
No electronic accuracy or equilibrium claim follows from a favorable sign.
