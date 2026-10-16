# Saved Mex force replay: arithmetic roundoff only

Slurm12680 completed the declared saved-only diagnostic in2s on1CPU.
All four raw energies, charge arrays, core gradients and point-charge gradients
match exactly. Recomputed cap scatters differ by at most
8.881784197001252e-16 kcal/mol/angstrom; direct scatter of the saved raw arrays
matches the parser’s fresh scatter exactly. The largest discrepancy occurs only
in floating-point derived quantities, not in molecular outputs or source maps.

The failed static12673 remains preserved. A separate v2 will retain all exact
raw checks and use a declared float64 roundoff bound for derived map arithmetic.
The accepted original mapped gradients remain the reused inputs. No molecular
settings, state, geometry, physical admission or scientific precision gate changes.
No new SCC, MACE, MM or solvent calls occurred in this diagnostic.

Actual diagnostic: `workspaces/lanm_coupled_response_20261016/mex_replay_diagnostic_v1/execution/RESULT.json`.
Runtime/package/build information and every indexed discrepancy are retained there.
NumPy2.4.6 in the pinned MACE environment was used; matching Python binary hashes
alone do not establish identical environments. No environment was modified.
