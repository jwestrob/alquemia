# Exact solvent selection of adaptive candidates

The parent reviewed the prepared eight-endpoint angular pilot and launched it as
1209857 after the completed shared-pool integration1209840/1209845. The four-case
shared pool did not change any prior own-proposal score. Its primitive remains the
common selection rule for the adaptive test. This continuation is within Jacob's
September22 execution handoff and the frozen angular proposal plan.

For each of those same four sources, retain the completed origin and both earlier
proposals. Add each actual successful, physically admissible adaptive proposal;
deduplicate only under the existing1e-12-Angstrom numerical-copy rule. Both metal
rows must contain this identical candidate set. Preserve chemistry, original
state-specific charges, native OMOL and exact primary native GFN2 solvent recipe.
No minimization across sources or chemical states. If a required new proposal or
matrix cell fails, the expanded result is unavailable; earlier scores remain
separately named. Do not discard the unsuccessful metal's candidate requirement.

Reuse each adaptive proposal's own real native MACE value after checking its
state/coordinates/model. Compute the competing metal at that same geometry and
both native GFN2 media for every newly admitted row/candidate cell. Upper bounds:
eight fresh MACE singlepoints and32GFN2 singlepoints, zero extra optimizations or
DFT. Reuse all earlier matrix cells. Use the existing warm GPU and64CPU8x8 solvent
runner. CPU work may run on an available GPU-partition host without requesting a
GPU; record actual hardware and allocation receipts.

The separate protocol is
`nikasha_adaptive_angular_common_geometry_native_OMOL_GFN2_ALPB_v1`.
Retain mathematical row minima and the existing0.10-kcal/mol operational origin
tolerance separately. Report individual native/solvent/composite works, selected
geometries, paired contrast changes, physical movement, boundary and residual-load
diagnostics. These are finite proposed states, not equilibrium populations or
stationary composite minima. No change of pH, water inventory, entropy or spring.

This four-case pilot has no new absolute reference. Report old-band transfer
explicitly as transfer only, plus raw contrasts/differences. 1H4I and4MAE are known
class controls; both PLM labels remain unknown. Moving a PLM value towardLa is not
proof of accuracy. A useful mechanical response or a diagnosed boundary issue
can justify a declared next test; silent scope expansion cannot.

Implementation reuses `scripts/nikasha_pool.py prepare-expansion`, followed by its
existing validate/execute-mace/collect operations. All actual source/receipt paths
are required arguments. Source-composition differences are rejected, never filled
with zero corrections. Larger-reference calibration and fold tests are subsequent
decisions; this pilot alone does not support production promotion.
