# Native fixed-charge solvent forces: passed on real LBT3 fixtures

Job12566 completed all26 declared native component evaluations and every frozen
gate. The returned gradient differentiates the same fixed-charge ALPB energy;
this permits a cheap solvent-force component without whole-protein SCC. It does
not validate a combined protein Hamiltonian or a metal-preference prediction.

Four archived water-charge endpoints reproduce all native components within the
declared population-rounding bounds. Eight directional checks on La/Dy metal-x
and Asn5 OD1-x agree with centered energy differences at both step sizes; largest
recorded residual is1.058e-9Eh/angstrom. Repeat is exact; joint translation changes
energy by1.11e-16Eh and the largest gradient component by1.93e-17Eh/Bohr.
The original tolerances were not changed. RESULT.json contains every residual.

The isolated wrapper reuses the pinned existing xTB6.7.1 libraries. Input charges
are held fixed; these forces exclude electronic charge-response derivatives.
Native energy includes electrostatic, hydrogen-bond, surface and reference terms.
The new energy is distinct from self-consistent GFN2(ALPB)-GFN2(vacuum).
No SCC, MACE, DFT, optimization or new molecular parameters were used.

Measured cost:5schedulerseconds on1CPU,5allocatedCPU-s; observedTotalCPU2.807s,
MaxRSS147920K. Application4.173s includes compilation and26component calls.
This is a256-atom fixture measurement, not whole-protein throughput.

Next: define a consistent full physical charge vector and combined interaction
ledger for Mex. Joint rotation and full mixed-model directional derivatives
remain untested. Do not project cap populations or add a second Coulomb term
without an explicit energy/force definition. Exact command/config/binary pins
are in the immutable workspace SUBMISSION.json and full result linked by RESULT.
