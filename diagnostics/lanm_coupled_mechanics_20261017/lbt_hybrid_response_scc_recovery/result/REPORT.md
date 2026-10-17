# Embedded solver recovery qualification

Four calls and all three predeclared agreement checks pass. This qualifies one
failed embedded GFN component numerically; it does not resume the search or
establish an electronic ground state.

The failed trial at initial Broyden damping 0.05 and 0.10 gives identical printed
energies, with maximum mapped gradient difference 1.185e-6 kcal/mol/Å. The two
successful controls reproduce their previous energies within 9.664e-13 Eh and
mapped gradients within 1.186e-5 kcal/mol/Å. Gates remain 1e-7 Eh and 0.005
kcal/mol/Å. The predeclared primary remains 0.05.

All 172 QM atoms, 84 point charges, coordinates, state, parameters, temperature,
accuracy and 500-iteration limit are unchanged. SCC iterations are 266/57/130/88;
all four native calls terminate normally. Small-gap and frozen-f limitations remain.
No MACE, MM, ALPB or DFT was run. Preflight13033 cost3 allocated CPU-s;
run13036 cost53 allocated CPU-s (application52.55s, native child CPU48.20s).

The original failure remains immutable. A separately reviewed continuation must
reuse the successful components, calculate only missing ALPB for the failed trial,
and resume the original checkpoint/history with the original80-step/200-evaluation
budget. Any further SCC failure stops that finite continuation. No global fallback.
