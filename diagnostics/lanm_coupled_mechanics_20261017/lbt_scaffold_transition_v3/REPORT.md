# LBT mechanics: several classical terms drive the disputed motion

19 October 2026 UTC. The exact four-state audit is complete. All original
input, parameter, energy and force closure gates pass. Job 14787 used one CPU
for 20 seconds; the two preserved zero-query startup failures add 5 allocated CPU-seconds.
There were 36 OpenMM contexts and 816 energy/full-force queries, with no new
electronic calculations, structures or optimization.

## What the audit establishes

The classical model lowers its energy along the tested deformation through
several interaction classes. At the same before/after coordinates, this lowering
is almost identical for La and Dy. The saved whole-peptide reference still gives
a cost of +12.279593/+11.722656 kcal/mol, while the complete La hybrid gives
−28.244859 kcal/mol. Full after-Dy hybrid work remains unavailable.

| Classical contribution to after−before work | La, kcal/mol | Dy, kcal/mol |
|---|---:|---:|
| Retained bonded interactions | −10.682960 | −10.682960 |
| MM Lennard-Jones | −6.180039 | −6.180039 |
| MM Coulomb | −3.189036 | −3.189036 |
| Electronic-region/MM Lennard-Jones | −3.229524 | −3.228347 |
| Total classical contribution | −23.281560 | −23.280382 |
| Bonded interactions touching either terminus | −9.532816 | −9.532816 |

The last row overlaps the first; it must not be added to the total. The terminal
union contains 187 unique interactions. Its lowering is distributed among bonds,
angles, periodic torsions and backbone correction maps. The largest individual
contribution is −1.338957 kcal/mol. These are complete interaction energies, not
atomic or residue energies, and the audit identifies no counting or parameter-copy
defect under its gates.

## Both the force offset and its change matter

For La, the hybrid's full-system force-error RMS relative to PBE0 plus the same
ALPB term is 10.141916 before and 8.850041 after, in kcal/mol/Å. The RMS of the
change in that error is 4.463637. Whole-peptide MACE has endpoint errors of
3.511987/3.516438 and an error-change RMS of 0.530091 on this fixture.

Thus the hybrid discrepancy changes during motion as well as being large at
the start. A constant force offset alone is insufficient. Two endpoints cannot
determine a Hessian or uniquely separate equilibrium-position and stiffness
errors. The nearly metal-independent classical contribution also cannot by
itself explain a La/Dy preference at these shared coordinates.

## Decision and next work

Keep the challenged hybrid optimizer closed. The evidence does not justify
tuning one terminal bond or angle to the reference. The next mechanics design
should test a consistent treatment of the scaffold and its coupling to the
electronic region, using this saved deformation as an independent acceptance
fixture. Record a finite intervention and its expected diagnostic outcome before
new evaluations; retain the full four-state denominator and original failures.

Before-Dy whole-MACE forces and after-Dy full-hybrid forces remain unavailable.
All reference limitations persist: conditional free termini, vacuum dianion,
frozen-f treatment, nonstationary endpoints, SCF qualifications and no independent
whole-reference derivative qualification. This result admits neither a metal
preference nor production promotion. The overall mechanics goal remains active.

The [generated tables](result/REPORT.md), [complete result](result/RESULT.json)
and [full force transitions](result/FORCE_TRANSITIONS.json) retain all 187 terms,
both frozen projection directions and all 256 atom vectors. Every raw query and
serialized system remains at the immutable execution paths recorded in RESULT.
