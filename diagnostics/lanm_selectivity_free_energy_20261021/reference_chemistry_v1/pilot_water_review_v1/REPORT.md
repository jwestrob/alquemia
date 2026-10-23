# Independent first-water result review

Both frozen-f PBE0 water endpoints pass the unchanged nine native admission
checks. Their energies and full four-atom gradients are available. This is one
shared geometry, not a water-exchange free energy or a physical-state validation.

| Quantity | La | Dy |
| --- | ---: | ---: |
| MACE minus reference force atom-vector RMS, kcal/mol/Angstrom | 1.454888 | 6.176160 |
| Reference metal radial force, kcal/mol/Angstrom | +4.167242 | -32.664979 |

Positive radial force points away from the donor. The RMS discrepancy in the
matched **Dy-minus-La force difference** is **5.068962 kcal/mol/Angstrom**.
This establishes a metal-dependent discrepancy on this fixture; neither method
is established as physically correct by the comparison.

The optional residual field in the frozen runner is empty because it looks for
YES/NO labels; these ORCA outputs instead print numerical residual/tolerance
rows. RESULT.json supplements all six rows directly from the raw outputs,
without editing native results or changing admission. Both orbital-gradient
residuals pass their printed tolerance; several density/DIIS criteria fail.
Thus `all_printed_residuals_pass` is false for both. In particular, the Dy
orbital-gradient residual is1.9601e-6 against2e-6. These diagnostics must remain
visible in later force qualification; native convergence is not sufficient.

No molecular calls were made in this review. Saved-array arithmetic used one
Python process with all numerical thread controls set to one. Exact native and
MACE result pins and complete force-difference arrays are in RESULT.json.
