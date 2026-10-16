# Matched-basis first scout

Both matched 7CCO/La searches execute and preserve physical geometry. The
seven-coordinate control lowers the already relaxed energy another 0.149229
kcal/mol; the ten-coordinate candidate lowers it 0.655434 kcal/mol. This is
one metal/source, not an improved preference result. Neither endpoint is
stationary. Both stop after finite searches that include geometric rejections
and energy-tested fallback directions.

All eight source/metal/arm start preflights and the source-derived Jacobian
checks pass before model loading. Control takes two accepted fallback steps,
candidate three; all covalent, contact, mapping and state gates are preserved.
No bounds, energy settings, proton/water inventory or source reference changed.

Actual new work: 9 MACE and 18 GFN2 calls, one warm model; job12556 completed
in128s on12CPUs/oneH200 (1536 allocated CPU-s,128 GPU-s). Allocation is not
measured utilization. See SCOUT_ACCOUNTING.tsv and SCOUT_RESULT.json for source
pins, per-arm works and residual forces.

The frozen continuation gate is execution and physical admission, not an energy
sign. It passes. Remaining six searches are job12558, oneGPU/12CPUs, completion
watcher1972597, configf5d82b06 in lbt3_basis_remaining_v1. No new paired pool
result is implied; both metals will cross-score the common compatible menus.
