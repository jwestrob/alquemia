# Native interaction decomposition — declared before execution

Under the [autonomous research agreement](AGREEMENT.md), examine which physical
interaction terms are missing when Asp303 is classical. The preceding completed
diagnostic found an8.06kcal/mol monopole error and only0.974kcal/mol change from
relaxing the core density in the permanent field, leaving9.67kcal/mol without
solvent. Those are consumed development results; this is not a blind test.

Run exactly **two native ORCA6.1.1 EDA-NOCV tasks**, one La and one Ca, using the
existing54-atom qm36 coordinates. Fragment1 is the47-atom qm33 core; fragment2
is the actual seven-atom capped Asp303 sidechain already present in qm36. Source
IDs establish membership; shared heavy atoms and caps must match coordinates.
Fragment1 charges are−1/−2 forLa/Ca; fragment2 is−1; total−2/−3. All singlets.
No new cap or atom is introduced. This isolates a real existing fragment pair.

Both fragments and adduct use native r2SCAN-3c, DefGrid3 and NoAutostart. No MBIS,
solvent, external field, geometry change, basis/ECP substitution or threshold
fit. The native EDA operation performs the fragment SCFs internally: at least
six SCFs across the two tasks, plus any package atomic-reference work. Count
and retain actual generated inputs/outputs. Two16-rank workers in a32-CPU
allocation; likely minutes at this size, with actual costs recorded and no
CPU/time stopping budget.

Report the Ca-minus-La change in electrostatic, Pauli, exchange-correlation,
orbital-relaxation and dispersion terms. ORCA reports Pauli and its separate
delta-XC term; do not omit the latter. Retain native Hartree values and convert
once. Track gCP separately because this native composite correction may be
outside the EDA table. Check component closure and compare total interaction
with the independently parsed adduct-minus-fragments energies. Fragment2 is
identical in both metal cases, so its energy must cancel within0.01kcal/mol.
Native adduct/fragment1 energies must agree with their archived gas endpoints
within0.01kcal/mol before attributing a component difference.

EDA's orbital term mixes polarization and charge transfer; it is not a unique
measurement of either. This experiment will guide a physically motivated next
approximation, not add a fitted scalar to predictions. A dominant short-range
or orbital contribution would argue against solving the problem by finer
classical dielectric meshes alone. It does not validate a particular cheap
polarization model.

If native r2SCAN-3c/ECP/composite accounting is unsupported, retain an explicit
unsupported result. Do not substitute another functional just to obtain a
decomposition. The production baseline and previous failed gates stay intact.

[ORCA6.1 EDA-NOCV interface and component definitions](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/nocv.html).
The retrieved page is preserved in the workspace source_inspection directory.
