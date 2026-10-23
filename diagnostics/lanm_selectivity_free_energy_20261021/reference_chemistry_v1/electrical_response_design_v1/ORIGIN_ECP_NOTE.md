# Native dipole origins, effective charges and rounding

**Confirmed for the actual ORCA 6.1.1 outputs:** the unlabeled `Position of electric origin` values are in **Bohr**, although explicit origin input follows the geometry input units. This conclusion uses native arithmetic, not an assumed equivalence of input and output units.

[ORIGIN_ECP_EVIDENCE.json](ORIGIN_ECP_EVIDENCE.json) pins three raw outputs, property files, inputs and RESULT files. The property geometry explicitly declares Bohr. Its coordinates, combined with the masses in the native A.U. geometry table, give:

| Existing case | Reconstructed COM, Bohr | Printed origin |
| --- | --- | --- |
| W9 Dy full | z = 0.269042368716622 | z = 0.269042 |
| W9 Dy water 1 only | x = 4.541562705113483; z = 1.696874786506304 | x = 4.541563; z = 1.696875 |
| W9 Dy metal + water 1 | x = 0.453237969878511; z = 0.169344371819024 | x = 0.453238; z = 0.169344 |

The real Dy center has native `ZA=11`, reflecting 66 minus the 55-electron ECP core. Ghost Dy and ghost O/H centers have zero ZA and zero mass. Reconstructing `sum ZA*(r-COM)` reproduces the full-precision native nuclear dipoles within 2.8e-12 au across these examples. Using atomic number 66 for the real Dy dipole would fail this closure. The property electronic, nuclear and total vectors also close and reproduce the nine-decimal output vectors to their printing precision.

For moving the reporting origin O to zero, with O in Bohr and moments in atomic units:

```
mu_total(0)      = mu_total(O)      + Q*O
mu_nuclear(0)    = mu_nuclear(O)    + sum(ZA)*O
mu_electronic(0) = mu_electronic(O) - Ne*O
Q = sum(ZA) - Ne
```

The full W9 cluster has sum(ZA)=101, Ne=98, Q=3; metal+water has 21/18/3; water-only has 10/10/0. A neutral fragment's total moment is origin invariant, but its separate nuclear and electronic contributions are not. No frozen-core electronic density is added back to these valence/ECP moments.

Six-decimal origin printing contributes at most 0.5e-6 Bohr per component, and nine-decimal dipole printing contributes 0.5e-9 au. Thus converting the full W9 printed vectors adds per-component bounds of 1.5005e-6 au (total), 5.05005e-5 (nuclear), and 4.90005e-5 (electronic). Linear combinations must propagate every coefficient and term. Root's bounded printed-data analysis can retain these bounds and avoid interpretation near them; no expansion to a new property-file analysis is required. The separately pinned full-precision files are available if printing later becomes limiting.

The official manual establishes the default center-of-mass origin and the property-file dipole fields; it does not explicitly label this output-origin line. The Bohr/output and effective-charge assertions above are independently verified from the actual native A.U./Bohr geometry and nuclear moments. [Electrical properties](https://orca-manual.mpi-muelheim.mpg.de/contents/spectroscopyproperties/electric.html), [property-file schema](https://orca-manual.mpi-muelheim.mpg.de/contents/utilitiesvisualization/property_file.html), [ECP conventions](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/basisset.html).

This is a units/provenance check on three existing calculations, not admission or interpretation of the full common-basis comparison. Zero new molecular calls; no finite-field release.
