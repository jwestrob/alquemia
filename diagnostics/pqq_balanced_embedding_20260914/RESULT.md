# Balanced PQQ embedding result

**Frozen primary-map answer: NO. Final answer: NO.**

| Structure | embedded score 3.3 | embedded score 3.6 | R_PC |
|---|---:|---:|---:|
| mxaf | +76.761 | +8.777 | -67.984 |
| c5_monomer | +78.311 | +15.197 | -63.114 |
| c5_dimer | +112.656 | +49.118 | -63.537 |

D_PC monomer = +4.870; D_PC dimer = +4.447 kcal/mol.
Sensitivity branch triggered: no.

Both selective-margin gates passed, but the generic MxaF gate failed
catastrophically. Relative to the bare calculation, embedding changed the
3.3-to-3.6 boundary response by mxaf -70.188, c5_monomer -67.290, c5_dimer -68.514 kcal/mol.
The near-common jump despite exact total-charge closure shows that a static
embedded-QM energy is not continuous when Asp(-) changes representation from
fixed MM charges to an explicit QM fragment. This model is rejected for
production scoring; the non-triggered sensitivity map cannot rescue its gate.

Machine-readable energies and raw-output hashes are in `result.json`; complete
charge, repair, and distance ledgers are in `preparation.json`.
