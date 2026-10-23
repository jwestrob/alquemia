# Peptide fragment preparation: explicit NME hydrogen identities

This version retains the four fragments, eight metal configurations, exact central-residue/deposited-metal coordinates, protein parameter file, charges and zero-energy scope in [PEPTIDE_FRAGMENT_SPEC.md](PEPTIDE_FRAGMENT_SPEC.md).

The original preparation is preserved at `peptide_fragments_v1/PANEL.json`, SHA256 `131f46e2c505d1a90638452091ba777c936c299f7fa60e19674cd28813b72953`. It stopped before parameter analysis or energies: native `pdbxyz` used the supplied NME amide atom named `H` as a fallback for missing methyl `H1`, producing two atoms at the same coordinates. The exact-coordinate uniqueness check rejected that structure.

For each NME cap, retain the neighboring deposited N, CA and amide H as before. **Also retain its deposited Cα hydrogen explicitly as NME H1**: use HA for the Thr4/Asn6/Ile10 neighbors and HA2 for Gly12. This is an input-name/retained-cap-atom repair. The amide H remains attached to N; H1 remains attached to the methyl carbon. Native Tinker builds only the remaining methyl hydrogens. No central residue, metal, parameter, tolerance or optimized geometry changes.

Run at most four corrected `pdbxyz` and four parameter-only `analyze GP` invocations, one thread each. Check:

- Unique atom coordinates, a reciprocal covalent graph and native type valences.
- Every central-residue atom and retained cap atom joins exactly to its deposited coordinate.
- NME amide-H and methyl-H1 map to distinct native atoms/types231 and233, attached to nitrogen/carbon respectively.
- Each newly constructed atom is a cap hydrogen with one carbon parent and native1.10 Å bond length (within1e-6 Å serialization error).
- Native donor types143/161/153/11 remain unchanged. No undefined parameters are printed.
- All native printed monopoles are present and their total agrees with the formal fragment charge within the existing printed-rounding envelope.
- Matched La/Dy copies have identical nonmetal atoms and coordinates, charges+2 for Asp/Glu or+3 for Asn/Trp, and no metal covalent bonds.

All energies, forces, minimization, dynamics and reference/model scoring remain unrun. This preparation makes the proposed transfer-qualification fixtures concrete; it does not select or validate a protein metal Hamiltonian.
