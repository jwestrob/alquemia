# Water force qualification: result

**Changing the quantum method or applying counterpoise does not resolve the La-water mismatch.** At this fixed geometry, AMOEBA pulls the water toward La; both PBE0 and canonical CCSD(T) push it away. The paired La/Dy response also remains too small in AMOEBA.

## Matched comparison

These are coordinate derivatives dE/dq in kcal/mol/Å, with positive q moving rigid water **away** from the ion. Positive derivatives mean inward water forces; negative derivatives mean outward forces. The prespecified h = 0.005 Å summary is eligible because both step sizes agree within the declared 0.1 resolution.

| Method | La | Dy | Dy − La |
|---|---:|---:|---:|
| AMOEBA analytic | +12.139 | +32.808 | +20.670 |
| PBE0-D4, uncorrected | -4.172 | +32.661 | +36.833 |
| PBE0-D4, counterpoise | -4.833 | +31.901 | +36.734 |
| CCSD(T), uncorrected | -9.765 | +28.417 | +38.182 |
| CCSD(T), counterpoise | -12.463 | +25.833 | +38.296 |

The paired AMOEBA residual is −16.064 against PBE0-D4 CP and −17.626 against CCSD(T) CP. The quantum methods differ by 1.562 on this paired derivative. Dy’s close agreement with uncorrected PBE0-D4 alone is therefore insufficient: its AMOEBA-minus-CC CP residual is +6.975. The La sign conflict is robust across all four displayed quantum comparisons.

D4 derivative contributions are −0.639 for La and +0.142 for Dy; they were evaluated only on real complexes. Unchanged isolated-water and ion dispersion cancels in the derivative. Electronic and D4 terms remain separate in the full result.

## Numerical and provenance review

- All 48 native cells pass their frozen admission checks. Independent extraction from the raw output reproduces every derivative and confirms final-energy/CC triples closure. The original central results were reused.
- Largest individual two-step change: 0.01599; largest paired change: 0.00533 kcal/mol/Å. PBE0-D4 uncorrected finite differences agree with saved analytic forces within 0.02080 at h = 0.010 and 0.00481 at h = 0.005.
- All eight displaced PBE0 complexes retain the central COSX grid/final-integration stage signatures. The expected EnGrad-to-SP change and native exchange corrections are preserved. This is a sampled directional consistency result, not a universal gradient guarantee.
- All canonical CC residuals are at most 8e−9. T1 values range 0.001375–0.010105; raw amplitudes, triples and SCF flags remain available. No new diagnostic cutoff was imposed.
- Native cost: **48 calls, 7,324 allocated CPU-s**, zero GPU use and zero repeated centers. All use four MPI ranks with the actual 16,000 MiB default grant. Saved collection took 13.54 s on one local thread. This independent review made zero molecular calls.

## Scope and next physical-model step

This remains a finite-basis, frozen-f effective-singlet comparison. The shared ECP/basis, omitted explicit Dy f-shell/spin-orbit physics, and absent bulk solvent prevent a claim that CC is physical truth or that this local derivative determines affinity. The result does establish that the current water-sector mismatch survives the tested method and counterpoise choices.

**Recommended bounded next step:** declare one paired La/Dy water-sector representability pilot, varying each metal’s water-contact radius and polarization damping while holding well depth, ion polarizability, charge and ligand parameters fixed. Test energies and radial response together using these existing water curves; check identifiability and the existing hydrated training controls before accepting any candidate. Both metals must be allowed to respond. Preserve the published baseline and keep reserved validation closed until a candidate is frozen. Pending donor CC comparisons provide an additional balanced exchange check; do not silently combine CP energies with uncorrected force targets. This report authorizes no fit, parameter change or new call.

## Evidence

- [Full collection](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/water_force_qualification_v1/collection_v1/RESULT.json).
- [Independent raw-output review and exact evidence pins](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/water_force_qualification_v1/ACTUAL_REVIEW.json).
- [Frozen specification](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/water_force_qualification_v1/SPEC.md).
