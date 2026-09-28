# Completed original-H states do not show an obvious charge-accounting failure

All four completed embedded Hans8DQ2 states contain806 printed occupied electrons,195 mapped atoms, and charge sums consistent with−1 within7e−6e (six-decimal printed-charge rounding). Their printed frontier gaps are positive:0.009533–0.009541Hartree. This does not prove orbital stability or qualify the physical electronic state, but the completed outputs do not show the occupied/virtual inversion seen during some running solver iterations.

| Endpoint | Mulliken metal charge/e | Loewdin metal charge/e | Printed gap/Hartree |
|---|---:|---:|---:|
| La_A | 1.536619 | 0.794799 | 0.009535 |
| Dy_A | 1.523712 | 0.822764 | 0.009541 |
| La_B | 1.531869 | 0.794956 | 0.009533 |
| Dy_B | 1.519425 | 0.822926 | 0.009539 |

The schemes even give opposite signs for the metal-charge change under substitution. These partitioned electron populations are not oxidation states. The largest grouped La→Dy population change is0.015085e for Mulliken (Asp85 in A) or0.027970e for Loewdin (metal in B). A→B changes are smaller (largest grouped change0.005250e Mulliken,0.002569e Loewdin). No abrupt population change is apparent across this prescribed small Asp85 motion; no universal acceptable-charge criterion is imposed.

**Decision:** do not diagnose the repaired runs as charge leakage or a wrong oxidation state from their initialization warnings alone. Continue the live numerical recoveries. If completed repaired states become available, their own populations and residuals can be compared; old-H results cannot establish the cause of repaired-state convergence difficulty. Full SCF stability, field/boundary sensitivity and physical-state accuracy remain distinct questions.

The source hashes, per-atom populations, grouped source populations and printed orbital data are in POPULATION_RESULT.json. POPULATION_AUDIT.py checks normal termination, SCF convergence marker, mapping hash, ordered atom inventory/elements and charge closure before emitting data. All four real outputs passed; no mocked output or new electronic call was used. This is parser/accounting evidence, not new electronic qualification. Population analysis ran on the login host as a short file-reading operation, no Slurm/GPU allocation or molecular calculation. Physical quantities are limited by the precision printed in ORCA output.

Reproduce into a fresh output path:

```bash
python diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/POPULATION_AUDIT.py --manifest workspaces/metal_environment_response_20260926/lady_frozen_embedded_hans_v2/manifest.json --output workspaces/metal_environment_response_20260926/population_audit_repeat.json
```
