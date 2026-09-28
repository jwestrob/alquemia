# Repaired Hans-La completes; Glu91 load decreases

The same-metal orbital-seeded repaired Hans8DQ2 La_A endpoint completed normally and passed the existing collector: energy−5728.001479272300Hartree,25SCF cycles, native atomic and point-charge gradients available. Worker1220332 took5183s on112allocated CPUs; collector1220333 took4s on1CPU. Total580,500allocatedCPU-s,0GPU. This is research validation cost, not a claimed affordable production rate.

Exact physical force assembly with matching repaired classical terms passed source-ID, coordinate, atom-order, gradient-energy and sign/unit checks. The already-declared Glu91 CB–CG tangent gives:

| Contribution | Normalized gradient/(kcal mol−1 Å−1) | Angular gradient/(kcal mol−1 rad−1) |
|---|---:|---:|
| Embedded electronic | 0.104286 | 0.384603 |
| Classical | 0.262097 | 0.966606 |
| Total finite model | 0.366383 | 1.351209 |

The original-H La value was0.948861 normalized or3.465243kcal/mol/radian. Each projection uses its own actual physical hydrogen coordinates and tangent; heavy coordinates/chemical inventory are retained. Repair reduces this one-metal load. Do not carry the original-H force magnitude unchanged into the new preparation, and do not infer the metal differential while repaired Dy remains unconverged. The separately prepared Glu91 displacement matrix remains unsubmitted pending that evidence.

This is not a minimum, affinity, population, full-hybrid derivative validation or solvent-consistent relaxation. Original and repaired energies describe different hydrogen/field coordinates; their difference alone is not metal selectivity. No new quantum calls were made for force assembly or projection.

Actual collection: workspaces/metal_environment_response_20260926/lady_repaired_la_seeded_v1/FINAL_COLLECTION.json.
Assembly: workspaces/metal_environment_response_20260926/repaired_DQ2_La_A_force_v1.json.
Machine-readable projection and exact source/implementation pins: REPAIRED_LA_PROJECTION.json. PROJECT_REPAIRED_LA.py accepts explicit assembly, motion-inputs, old-projection and fresh output paths; the adjacent implementation and original force-assembly CLI reproduce the analysis. Three remaining workers1220323,1220336,1220342 remain live. Their automatic comparison1220344 remains pending.
