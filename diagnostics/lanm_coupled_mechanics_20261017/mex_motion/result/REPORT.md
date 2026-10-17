# Mex saved motion and mechanics

Segment selection: candidate_for_feasibility_review.

| Metal | Bonded | MM LJ | MM Coulomb | QM/MM LJ | Total MM work (kcal/mol) |
| --- | ---: | ---: | ---: | ---: | ---: |
| La | -730.792190 | -389.119324 | -2014.301539 | -36.600565 | -3170.813618 |
| Dy | -724.010750 | -386.767826 | -2017.221829 | -35.646444 | -3163.646849 |

| Motion group | La aligned RMS / maximum (Å) | Dy aligned RMS / maximum (Å) |
| --- | ---: | ---: |
| peptide_outside_QM_heavy | 0.1845 / 0.8968 | 0.1926 / 0.9027 |
| peptide_QM_heavy | 0.1544 / 0.6138 | 0.1536 / 0.5245 |
| water_oxygen | 0.6521 / 1.9755 | 0.6889 / 2.0984 |
| target | 0.0982 / 0.0982 | 0.1612 / 0.1612 |
| Nd_spectators | 0.1447 / 0.1698 | 0.1539 / 0.1875 |

Candidate window: **A/77/, A/78/, A/79/**.
It accounts for 7.42% of pooled outside-QM peptide heavy-atom squared motion after whole-protein alignment.
La/Dy contributions (Å²): {'La': 1.6935925874933018, 'Dy': 1.889927539000897}.
Effective participating residues: 49.265139429021815.
Residues covering 50/80/90% of motion: {'0.5': 19, '0.8': 46, '0.9': 61}.

Global MM group lowering is not a residue energy decomposition. Water/peptide motion and full component force locations are recorded separately in MECHANICS.json.

The selected segment is one local stiffness-test candidate. Its coverage and QM connectivity do not establish that one expanded region can explain the full release; complete chemistry, charge boundaries and the local hydrogen-bond network still need review. No localization threshold, preference sign or new reference outcome selected it.

Searches remain finite and their raw force residuals remain explicit. No model evaluations or new reference launches were performed.
