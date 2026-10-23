# La baseline audit

## Finding

No implemented La typing, units, parameter mixing, or pair-override defect was found in the saved training calculation. The La force mismatch cannot be assigned wholly to force-field error: the current PBE0 target differs from the published calibration target. No fit or new molecular calculation was performed.

## Implementation evidence

- All 28 training parameter files reproduce the original SI AMOEBA09 parameter block byte for byte. Their key prefix reproduces the SI AMOEBA22 key (relative parameter path normalized); appended records are the unused Dy type and numerical output/polarization controls.
- Water types 36/37 map to classes 34/35, SI charges −0.51966/+0.25983 and polarizabilities 0.837/0.496. The La calculations do not use the separate Dy/GEM water model.
- La amide oxygen type 105/class 55, amide carbon type 104/class 53 and acetate type 403 are mapped correctly. Saved native output reports the intended VDWPR and POLPAIR overrides: classes 55/500 (3.3000, 0.3019), 53/500 (3.9541, 0.8295), 403/500 (3.6318, 0.2280); damping types 105/500 (0.2987) and 403/500 (0.1518).
- All four saved La dimer interactions reproduce the published SI to about 0.01 kcal/mol. La-water and bidentate-acetate Cartesian finite differences agree with native analytic gradients to approximately 1e−8 and 7e−8 kcal/mol/Å, respectively.
- Tinker forces are minus its kcal/mol/Å gradient. ORCA forces use minus its hartree/bohr gradient times 627.509474/0.529177210903. Saved training joins, rigid-body force sums, centroid torques and category weights are consistent. Shared internal water forces cancel in a rigid translation and cannot explain the observed radial force mismatch.

## Reference and representability limits

The original paper (DOI 10.1021/acs.jpcb.2c07237, Methods; local article lines 181–210) optimizes ligand structures at MP2, then fits distance-dependent **counterpoise-corrected DLPNO-CCSD(T)** interaction curves with relativistic ZORA/def2 and metal SARC bases, together with La hydration free energy. Calling the fitted energy target simply MP2 is incorrect. SI2 lines 17–32 specify the fragment/ghost counterpoise subtraction.

Our current force target is uncorrected frozen-f ECP PBE0-D4. At the exact central water fixture, La AMOEBA interaction is −91.271778, canonical CCSD(T) CP is −89.467424 and PBE0-D4 CP is −97.365454 kcal/mol. These same-basis data show target sensitivity; they do not establish an exact force reference.

The original SI independently reports monodentate acetate underbinding of 49.09 kcal/mol (QM −466.67 versus AMOEBA22 −417.58). The article discusses environment-dependent inequivalence of the two oxygens versus equivalent fixed oxygen parameters. This is a published representability limitation, separate from implementation correctness.

## Decision

Keep the La model and all failed eligibility criteria unchanged. Pause Dy-only fitting. The next declared test compares matched water radial derivatives from canonical CCSD(T), PBE0 electronic and PBE0-D4, with and without counterpoise, against saved AMOEBA forces. Full evidence paths and hashes are in EVIDENCE.json. No reserved validation outcomes were inspected.
