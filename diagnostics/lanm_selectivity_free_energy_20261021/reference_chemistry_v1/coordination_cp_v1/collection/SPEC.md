# Saved coordination-counterpoise collection

This collector evaluates only the frozen ROOT_DECLARATION cycle after the native results exist. It launches no chemistry and does not change any fitting target.

Keep exactly eight fragment slots: W8/W9 × La/Dy × real ion with ghost shell / real shell with ghost ion. Require all native checks, normal successful completion, exact input/geometry/core/basis/state, finite electronic energy, absence of fragment D4, source and launch pins, and re-parsing by the frozen native adapter. Keep every missing, failed or malformed slot explicit. Re-admit all four original full-complex results using their executed native collector and extract their electronic, D4 and final energies independently from raw output. Do not compute a partial four-complex cycle.

For `ΔX=X(Dy,W9)−X(La,W9)−X(Dy,W8)+X(La,W8)`, let `G=E_ion^fullbasis+E_shell^fullbasis`. Report:

- original electronic contrast `ΔE_electronic`;
- original, unchanged dispersion contrast `ΔE_D4`;
- original total contrast `C_original=ΔE_final`;
- basis-sharing term `ΔG` and correction `−ΔG`;
- alternative electronic contrast `ΔE_electronic−ΔG`;
- alternative total contrast `C_CP=C_original−ΔG`;
- unchanged model contrast and model-minus-original / model-minus-CP discrepancies.

Use the original conversion 627.509474 kcal/mol per hartree. Preserve original D4 exactly; no ghost dispersion is added or inferred. Final/electronic/D4 closure uses the existing 1e−8 hartree parser tolerance. Own-ion cancellation across sizes and own-shell cancellation across metals must be proven by actual prepared/source identities, not by assuming the newly calculated ghost-fragment energies cancel. No individual absolute CP-corrected complex energy is claimed without its own-basis fragment terms.

An unavailable native slot or invalid cancellation makes the corrected cycle unavailable. Available original-complex and fragment evidence remains in the receipt/report. Actual native costs are read from admitted receipts and distinguished from reused costs and the saved collector's one-thread CPU/wall cost.

The result is an alternative two-fragment reference diagnostic. It does not correct intra-shell water–water BSSE, establish basis/core convergence, qualify gradients or free energies, select a model, refit parameters, or release validation/sampling. A finite saved test set checks cycle signs, D4 preservation, finiteness, missing denominators and provenance failure; no molecular calls are included.
