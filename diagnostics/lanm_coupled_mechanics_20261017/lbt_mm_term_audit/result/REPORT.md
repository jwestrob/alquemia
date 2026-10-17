# LBT terminal MM force attribution

**All four fixed-coordinate cases pass. Bond stretching and angle bending dominate the terminal MM restoring force.**

Exact native parameter and term-inventory checks found no copying or double-counting defect. Full-system and per-term force closure pass; this does not establish that the MM stiffness is physically accurate.

## Components along the total MM force

Phe1 CA, kcal/mol/angstrom. Signed projections add to the total MM force magnitude; they are not independent energies or causal percentages.

| Source / metal | Total MM | Bonds | Angles | Torsions | CMAP | Nonbonded |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 7CCO / La | 52.0094 | 28.2615 | 22.5231 | 0.7223 | 0.0000 | 0.5026 |
| 7CCO / Dy | 48.2092 | 24.9580 | 22.1753 | 0.6052 | 0.0000 | 0.4707 |
| 7CCN / La | 53.1341 | 30.5330 | 18.9813 | 1.2426 | 0.0000 | 2.3773 |
| 7CCN / Dy | 51.1441 | 30.6170 | 16.6994 | 1.2267 | 0.0000 | 2.6010 |

## Evidence and limits

All208 local terms are retained in ALL_TERMS.tsv. Largest closure residuals: force 1.28e-13 kcal/mol/angstrom; energy 2.84e-14 kcal/mol.

Job13122 completed256 MM queries in16 Reference contexts,9 scheduler seconds on1 CPU. This report performs saved arithmetic only. All four original starts, parameters, charges, terms and coordinates are preserved. No DFT, MACE, GFN2, ALPB, optimization or production change was made.

The whole-peptide DFT scout13051 independently tests the accuracy of this force. This attribution alone supports no parameter change, geometry repair, force-field verdict or metal-preference claim.
