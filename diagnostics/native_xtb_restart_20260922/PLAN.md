# Approved native xTB restart continuity experiment

Jacob explicitly approved the next targeted native-solver restart and collective
scaffold round; root delegated only the eight restart calls to this agent on
2026-09-22. This agreement supersedes the proposed-only execution status of
`diagnostics/structure_informed_starts_20260922/RESTART_PROPOSED.md`. Its scientific
scope and pinned `RESTART_SOURCES.json` remain unchanged.

Run exactly old/new Q88JH5 La geometry × vacuum/ALPB × self/opposite-geometry seed:
eight singlepoints. Each seed is the successful same-medium native charge/multipole
file pinned by the source record. Destination XYZ, charge−2, singlet, 438 native
electrons, default native GFN2 parameters, 300K, MaxIter500, eight MPI ranks and
solvent treatment are unchanged. Remove only NoAutostart from source input.
No GBW or parameter override is copied, and no ordinary-SCF, Tight, standalone,
MACE, DFT or geometry calculation is added.

The scoped wrapper uses the existing ORCA executor. Preparation retains a separate
immutable seed; immediately before launch it verifies a fresh directory and copies
that seed to `endpoint.runtime.xtbw`. Record all pre/post files and per-attempt
receipts. Execution uses one CPU-only gpu-partition job, 64 CPU/128 GiB and eight
concurrent eight-rank tasks. No custom time/compute cutoff or retries.

The installed `orca_guess` contains `INITIAL GUESS: XTBRESTART`. Collection requires
this positive output marker, native mixer, source-matched parameters/electronic
state, charge closure and normal converged receipt. A generic successful exit or
presence of an overwritten seed is insufficient. Otherwise report
`restart_not_confirmed` or the actual failure; retain energies separately without
claiming a qualified restart. Missing cells remain null.

Report full self/cross energy matrices in both media and matched ALPB−vacuum
transfers. Keep the prior diagnostic tolerances: 0.1 kcal/mol per endpoint and
0.2 kcal/mol for differential-transfer continuity. Retain printed density residuals
without claiming they are the active native stopping criteria. No class labels,
ground-state selection policy, new reference, production change or further run
follows automatically. Every calculation remains an inspected technical diagnostic.
