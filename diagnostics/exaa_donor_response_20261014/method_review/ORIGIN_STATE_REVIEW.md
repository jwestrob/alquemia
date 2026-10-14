# Target La origin: numerical admission with electronic-state caveats

Read-only review of completed job4109. The origin is admissible under the declared native diagnostic rules. Its electronic-state uniqueness remains unestablished. No calculation, submission, input edit or change to acceptance was performed; proposal4110 was still incomplete at this review and was left untouched.

All line references below concern
`workspaces/exaa_donor_response_20261014/reference_v1/tasks/sample3_La_origin_vacuum/endpoint.out`.

## What completed

- The final NR macro5 has orbital residual **8.804145e−7** (1879), below the printed **1e−5** TRAH/TolG criteria (1646,1690). ORCA reports SCF convergence after97cycles (1883) and normal termination (6313).
- The final scalar is **−6073.096697325370Eh** (6108). It includes the native D4 and gCP corrections printed at6099–6104; the earlier SCF total at1894 is not a substitute final energy.
- The output confirms charge−2, multiplicity1 and828 explicit electrons (1625–1627), with La Def2-ECP replacing46 core electrons (1457). The atomic-guess multiplicity2 at1699 is an atomic La initialization, not the molecular state.
- Integrated alpha/beta densities are equal at414.000085568601 each (1908–1910); printed Mulliken atomic charges sum to−2 (2556). This supports the recorded electron/charge bookkeeping, not a physical cofactor oxidation-state measurement.
- The collection marks all17 native/identity checks true, retains the GBW and complete populations, and preserves the original optional-summary failure as `previous_TRAH_check:false` beside the corrected actual-threshold pass.

## Frontier orbitals and remaining limitations

The final orbital table prints orbital413 occupied2.0000 at **−0.007810Eh**, followed by orbital414 unoccupied0.0000 at **+0.092190Eh** (2343–2344). Unlike some PedE references, this table shows neither a positive occupied HOMO nor a frontier occupied/unoccupied reversal.

That table does not erase the near-final warnings. A negative Hessian diagonal and negative HOMO–LUMO gap−0.000097 are reported after macro4 (1867–1868). A further small negative-gap warning appears **after SCF convergence**, immediately before full diagonalization (1886–1887). Earlier large-gap/Hessian warnings are transient iterations; these late warnings warrant retention. The discrepancy between the late gap warnings and final orbital table has not been independently resolved. The printed table cannot establish a robust unshifted gap, validated Aufbau occupation or a unique stable ground state. The normal header enables level shifting (1667–1669), and the runtime explicitly sets a1Eh shift after an early small-gap warning (1744). The apparently orderly final table is therefore not independent electronic-state evidence; this review does not claim a particular remaining shift or diagonalization detail uniquely explains it.

The large final density/DIIS summary values (1920–1922) remain visible. They must not replace the actual ending TRAH/NR gradient when applying the declared admission rule; their presence also means this result should not be advertised as an all-criteria density-convergence or repeatability qualification. There is no independent stability calculation, alternate occupation/spin solution or physical redox validation here.

**Decision:** retain the accepted conditional origin and inspect the independently finishing proposal under the same rules. No new blocker to the declared paired scout is demonstrated by this origin alone. There is no donor work, paired-metal response or classifier conclusion until the required endpoints are complete and compared on Slurm.

## Immutable records

- Raw output SHA256: `7bce8d0c9ec59a7d0210224e2b85b89179ce911dd8b137aee6850f71caf7173c`.
- Native execution receipt SHA256: `5d3f081f35c3b282010e5bf7ff8c4cd533d8550062ac323e959e006ec9f72cae`.
- Collection: `workspaces/exaa_donor_response_20261014/reference_v1/attempts/4109_execute_sample3_La_origin_vacuum/COLLECTION.json`.
- Collection SHA256: `17302bb43143159c5790c47218c2cae937704fae6e078d7e4347978d1444387f`.
- Manifest SHA256: `d37733923c3203f41b73d5a3c1789d396e22f8da69e624df2c2d724d0141bcf1`.

Actual native output reports1h40m19.294s (6314); root owns scheduler/stage cost accounting. No runtime projection or follow-up solver modification is made here.
