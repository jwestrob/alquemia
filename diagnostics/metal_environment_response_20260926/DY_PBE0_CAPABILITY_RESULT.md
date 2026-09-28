# Isolated Dy PBE0: stopped without convergence

## PBE0 capability attempt stopped — 27 September 23:00 PDT

Root received the live near-root stagnation wake and stopped only worker1219868.
Twelve consecutive TRAH/NR macroiterations54–65 remained above5× the printed
1e-5 orbital-gradient tolerance, with an energy span4.83699e-7Eh. The residual
oscillates rather than improving reliably; nearly constant energy is insufficient
for a qualified analytic-force reference. No final energy/gradient is accepted.
The attempt ran about5hours before the alert; exact terminal cost follows in the
collector1219869 receipt. This was one isolated50atom explicit-f Dy sextet PBE0
calculation, not a protein relaxation or a preference comparison.

ROOT_STOP_DECISION.json in dy_pbe0_capability_v1 preserves the live evidence.
Keep collector and terminal watcher running. On completion record actual cost
and unavailable endpoint; do not restart this attempt or increase iterations.
Both native r2SCAN-3c and this PBE0 treatment encountered electronic convergence
problems in the isolated core. This rules out the protein field as a necessary
cause, but does not establish a unique cause, incorrect spin, or impossible Dy
chemistry. Numerical precision, electronic occupations and the source model remain
unresolved; the plateau alone does not distinguish them. No further method sweep
or dependent LanM expansion is justified by this result. PQQ remains untouched.

Next work uses saved outputs to define a specific electronic/numerical diagnosis
and a defensible reference strategy. The exact MACEPOL-EF weights remain unavailable;
there is no qualified field-aware ML candidate or within-series discriminator from
this campaign. Preserve the successful Ca/La environmental-response work separately.


Final scheduler accounting: worker1219868 cancelled after18101s on24CPUs; collector1219869 completed in4s on1CPU. Total434428 allocatedCPU-seconds, zeroGPU. Collector retains endpoint unavailable after cancellation; no valid quantum energy or gradient.
