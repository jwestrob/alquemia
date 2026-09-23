# Installed capability and primary documentation

Consulted 2026-09-22. The installed binary is pinned in each finite manifest;
actual calculation output determines what was executed.

- [ORCA 6.1.1 native xTB documentation](https://orca-manual.mpi-muelheim.mpg.de/contents/modelchemistries/semiempirical.html): native GFN2 supports ORCA solvent models including CPCM. Native xTB's default special mixer uses charge/multipole-specific shortcuts; ordinary SCF is separately selected by `UseXTBMixer false`. The documented default electronic temperature is 300 K. The native and external xTB interfaces must not be conflated.
- [ORCA 6.1.1 CPCM documentation](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/solvationmodels.html): the default cavity uses Gaussian surface charges, scaled vdW radii, scale 1.2 and a surface density of 5 charges/Angstrom^2. Bondi/Mantina radii are used with a special hydrogen radius. CPCM describes a self-consistent electrostatic reaction field; the additional SMD/CDS terms are not part of this pilot.

Historical installed CPCM defaults were read from the original A0A3F2YLY8 Ca/La
outputs under `diagnostics/pqq_pmdh_fixed_core_calibration_20260914/prepared/01_a0a3f2yly8-pqq-la_model/`.
These establish an expectation, not proof of native-GFN2 coupling. The new pilot
must itself print its cavity and reaction-field contribution. No local manual was
bundled with this binary; the official 6.1.1 and FACCTs 6.1 manuals were checked.

Prior numerical work is in `diagnostics/compact_qualification_20260920/REPORT.md`:
24 tighter-native endpoints passed; genuine ordinary-SCF MORead converged only one
of eight endpoints. This history motivates first retaining the working native
mixer and requiring matched vacuum controls if CPCM changes that solver.
