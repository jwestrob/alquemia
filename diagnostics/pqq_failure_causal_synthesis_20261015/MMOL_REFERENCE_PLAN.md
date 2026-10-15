# Reference the opposite-sign MMOL1770 accommodation

Declared before new numerical work, 15 October 2026 UTC. Jacob's active goal is
to establish why the remaining predictions fail or abstain, with autonomous
contained experiments. This is a development reference, not a production change.

## Question and exact work

Does an independent electronic treatment reproduce the Ca-directed accommodation
of XoxF4-2/MMOL1770 Ca-conditioned sample0? Its static origin is correctly called;
the current accommodation shift is -4.186093 kcal/mol, of which -4.004941 is the
MACE-vacuum component. This contrasts with ExaA's independently supported large
La-directed response. Previous MMOL Ca1 and La4 studies are different sources.

Use the exact admitted 193-atom, 13-cap, dry oxidized-PQQ3- context, its saved
origin/adaptive_Ca/adaptive_La common geometry pool, and both endpoint metals.
Six singlepoints total. No optimization, new geometry, protonation, water,
assembly, donor or context change. Both metals see the same physical geometry
in each column. The source1 peer inventory is retained for a potential later
comparison, but is not part of this six-cell execution.

## Hamiltonian and admission

Use native ORCA6.1.1 r2SCAN-3c, vacuum, DefGrid3/TightSCF/NoAutostart, without extra
dispersion. Ca charge -3 and La -2, physical singlets. Verify actual native
def2-mTZVPP, auxiliary basis, La46-electron ECP, electron inventory, convergence,
populations and normal termination. Preserve the existing corrected native TRAH
collector, allocation-aware renderer, execution receipts and GBWs. No GFN2
parameter file belongs in these DFT calculations. Native GFN2's historical La
shell-assignment error is not an alternative model or a retry route.

The numerical admission rules are inherited from the completed ExaA reference,
including actual printed TRAH thresholds. Preserve density/state warnings and
energy-only endings; this is not proof of a unique electronic state. Do not
silently loosen gates or accept a failed energy.

## Comparison

For each candidate q, report W_M(q)=E_M(q)-E_M(origin) and the same-path contrast
D(q)=W_Ca(q)-W_La(q), separately for native DFT and saved MACE. Report the complete
vacuum pool accommodation min_q W_Ca(q)-min_q W_La(q), each selected geometry and
every cell. Also show the works at the historical composite-selected geometries
as a separate diagnostic. Do not equate those mixed selections to a common-path
response or assume vacuum row minima reproduce solvent-based choices.

Use the earlier 1 kcal/mol interpretive resolution for discrepancies, not a
fitted acceptance threshold or a DFT error bar. Sign reversal or substantially
different same-path/pool work motivates a local surface/state diagnosis.
Agreement shifts attention to source geometry, omitted environment or state;
it does not validate the experimental preference or prove a bad fold. No MACE
decision bands are applied to raw DFT energies. Missing cells remain missing.

## Staging, resources and stop rules

1. Slurm real-input/preflight and archived parser tests, all six declared cells.
2. Two independent 12-MPI origin scouts, Ca and La. Review actual native endings,
   state evidence and runtime before expanding.
3. If executable and interpretable, four saved-proposal cells using the same
   declared protocol. Root owns submission and the continuation decision.
4. Slurm collection/comparison, full six-cell denominator and actual costs.

Recent compatible ExaA192-atom endpoints used 12 CPUs and 3670-6994 wall seconds
each (8 cells consumed 522696 allocated CPU-seconds). This provides a planning
scale, not a forecast or speed guarantee for this different electronic system.
The earlier PedE slow Ca endpoint took substantially longer. These are one-time
development calculations, not the routine scanner workload.

Use the standard shared queue and current scheduler limits with allocation-derived
memory, no priority boost, no memory/time override, and no interference with
other jobs. Keep automatic terminal/health wakes. A failed scout calls for a
specific diagnosis, not an automatic iteration/allocation increase. No broad
functional sweep or new high-level optimization is included.

## Ownership

Root integrates, submits and reports. Curie owns exact saved-source inventory
and independent mapping review. Gauss owns the focused adapter based on the
existing ExaA native runner. Heisenberg independently handles saved ExaA PQQ
force projections; those are not DFT force evidence. Production and historical
calibrations/results remain unchanged.
