# Evaluate the small source potential without FMM truncation

The first PCM pilot stops nonzero sources at the declared potential check,
before their forward energy solves. For example, 2FW0 Ca has maximum native
FMM-minus-direct potential errors of 1.3451e-7 au (coarse) and 1.9146e-7 au
(primary), above 1e-8. Monopole integral normalization passes at 8.9e-16.
The actual zero-source control reaches the solver and returns successfully.
Retain the first pilot, its failures and cost unchanged. This does not show
that PCM's reaction energies fail; the nonzero energies have not been computed.

The source contains only the existing small set of nonzero QM-projected charges.
Evaluate phi by the direct Coulomb sum at every ddX cavity point, using the
same source positions, charges and atomic units as before. Use `math.fsum`
for each point; audit against the existing vectorized direct implementation
within the unchanged 1e-8 au tolerance. This is a summation consistency check,
not independent evidence about the quantum density. Retain the native FMM
potential and its discrepancy as a separate diagnostic. Do not relax the gate.

The maintained ddX interface accepts externally evaluated phi; psi still comes
from its exact monopole helper for the identical source distribution. No new
PB engine, fitted charges, density, cavity, dielectric, physical model or
biological threshold. Keep the accelerated PCM operator, all numerical
parameters and all 20 planned solve roles from DDX_SOURCE_SELF_PLAN unchanged.
Use a new protocol/manifest for this numerical source-evaluation path:
`fixed_source_full_protein_ddPCM_GK_component_comparison_v1_direct_source_v2`.
Repeat the full declared inventory; count both pilots and do not reuse the
earlier zero solves across different implementation hashes. No DFT or MACE.
