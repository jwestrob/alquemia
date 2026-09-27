# Independent coarse partition-sensitivity diagnostic — 27 September2026

Question: does promoting the already-selected full Thr159 side chain change the
metal-dependent response to exactly the same A/B physical perturbation materially?
This is the handoff's expanded-region test, now explicitly a diagnostic rather
than a qualification claim. Fine rigid invariance remains failed and is not waived.
User requested continued work. Stop SCF/grid repair attempts that leave the same
rotation difference; do useful independent representation analysis instead.

Four endpoints only: Ca_A,Ca_B,La_A,La_B from the immutable expanded_region_v1
INPUTS.json (63QM atoms,9078 MM point charges). B moves HG1 in QM index57, so it
MUST use its own XYZ; unlike the small-core scout, A/B external files are identical.
Same real source and physical atom inventories, new cap/boundary mapping and
original N/C charge redistribution documented in preparation/expanded_region.
No optimizations, geometry rescue, new proton/water inventory or biological labels.

Use original native r2SCAN-3c NoAutostart DefGrid3 TightSCF EnGrad, DoEQfalse,
Ca-3 andLa-2 singlets, native46electron La ECP. Match the archived original
small-core scout protocol; refined/strict variations do not silently replace it.
Per-metal response is E_M(B)-E_M(A), and delta=(LaB-LaA)-(CaB-CaA).
Compare delta_expanded-delta_small; never subtract absolute totals of differently
sized QM regions as a binding energy. Small-core per-metal changes omit quantum
internal Thr energy, so their differences between partitions include that change.
Double difference cancels only genuinely metal-independent terms. Link chemistry,
redistribution and quantum response change together; do not uniquely attribute
an effect to polarization. Full hybrid mechanical/cross terms remain unavailable.

Interpretation frozen before energies: report signed changes and force projections
where actual mapping supports them, no success/failure classification. Small-core
rigid residuals~0.01kcal/mol and original0.05numerical target remain visible; these
are not a proven uncertainty bound on the expanded core. Near-zero/small effects
are unresolved without qualification. Large sensitivity is evidence against
assuming partition independence, not proof of accuracy or affinity. No automatic
LanM/other-core extension or production promotion follows these four cells.

Execute4workers x86MPI=344slots, --mem=0 exclusive verified RealMemory with usual
runtime headroom, normal priority, no PQQ dependencies and no test queue. Monday
cutoff remains. Arm collector and actual completion wake. Count all costs.
