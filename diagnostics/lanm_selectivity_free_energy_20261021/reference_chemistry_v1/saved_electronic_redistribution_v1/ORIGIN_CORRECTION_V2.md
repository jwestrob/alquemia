# Pre-execution origin correction

Independent source review found that ORCA reports each fragment dipole about
its own center of mass. The original source and declaration remain unexecuted.
V2 converts every vector to coordinate origin zero before taking differences:

- total: mu(0) = mu(O) + Q O;
- electronic: mu_e(0) = mu_e(O) − N_e O;
- nuclear: mu_n(0) = mu_n(O) + (N_e + Q) O.

Verify the printed origin equals the mass-weighted native coordinates in the
explicitly labeled atomic-unit table, including zero-mass ghost centers. Native
effective nuclear charges must sum to N_e + Q. Independently reconstruct the
nuclear moment from those coordinates to check sign, ECP convention and units.

The origin/coordinate table prints six decimals and dipole vectors nine decimals.
The prospective per-vector conversion uncertainty is therefore
5e-10 + |charge|*5e-7 atomic units per component. Propagate absolute coefficients
through each C/P/B difference. Nuclear cancellation uses this bound, replacing
the unexecuted1e-6 bound that ignored origin rounding. The origin-versus-COM
check includes1e-6 atomic units for the rounded coordinates/origin plus the
propagated three-decimal mass rounding: sum(5e-4*|r_i-COM|)/(total_mass-number_of_real_atoms*5e-4).
Zero-mass ghost centers are excluded from mass-rounding uncertainty. Comparing the reconstructed
nuclear vector also includes effective nuclear charge times5e-7 for coordinate
rounding. These are derived print-precision bounds, not post-result fitting.

All population-sum and algebra tolerances remain unchanged. V2 also validates
the independent review's status and exact collector join, exact80input inventory,
and writes structured failure/cost records. No molecular call or scientific
acceptance criterion changes. Comparative dipoles are unavailable if an origin,
nuclear, or rounding check fails.
