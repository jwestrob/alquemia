# Saved native grammar recovery

The original 18321 calculation terminated normally. Its immutable result failed
only the global TCutMKN/TCutDO comparisons and an unobserved LMP2 completion
phrase. All other physical inventory, auxiliary basis, CC residual, iterative
triples, energy closure and native termination checks passed.

This version reads saved outputs only. It preserves the original adapter,
result, input, START, release and raw streams. It changes no molecular setting,
numerical threshold, geometry, resource allocation, energy or method criterion.

The observed ORCA 6.1.1 trace distinguishes a semi-local MP2 prescreen from the
full local MP2 stage. Require the primary DLPNO SETTINGS and full stage to match
the declared tier; retain and check the prescreen's printed default scale
factors (MKN 10, DO 2, PNO 1). Require exactly the observed two-stage order.
Unknown stages or malformed numeric tokens remain unavailable.

The full LMP2 terminal row must explicitly end `=> CONVERGED`. Its finite DE
and RMAX must satisfy the printed, unchanged 1e-7 and 5e-7 tolerances. Require a
complete, sequential iteration trace and exactly one terminal marker. This is
stronger evidence than a completion phrase alone. All remaining frozen parser
checks and exact executed-input/source lineage remain required.

Scope: saved admission of the four original scouts only, zero native calls.
Single-cell recovery does not qualify the eighteen-cell approximation panel.
The old remaining-control gate cannot accept a failed original RESULT: any
subsequent execution handoff must explicitly bind this reviewed parser recovery
and preserve the original producer. Root owns that release and all submissions.
