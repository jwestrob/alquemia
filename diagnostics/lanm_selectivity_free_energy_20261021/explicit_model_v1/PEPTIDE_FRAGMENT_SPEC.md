# Four native protein donor contexts: geometry-only declaration

Prepare four capped peptide fragments from the same 9CEQ model1 used for the full amidated topology:

| Fragment | Target donor | Source central residue | Native donor type |
| --- | --- | --- | --- |
| ACE–Asp–NME | Carboxylate | Asp3 | 143 |
| ACE–Glu–NME | Carboxylate | Glu11 | 161 |
| ACE–Asn–NME | Side-chain amide | Asn5 | 153 |
| ACE–Trp–NME | Backbone carbonyl | Trp9 | 11 |

The central residue retains all deposited coordinates and hydrogens. Its preceding residue's CA/C/O supply the ACE methyl/carbonyl coordinates. Its following residue's N/CA and amide H supply NME coordinates. Rename cap atoms only; native `pdbxyz` constructs the missing methyl hydrogens. There is no geometry optimization or selection by energy. Record every retained and constructed atom and check exact central-residue coordinate joins.

Use the unchanged pinned `amoebabio09.prm` for native protein multipoles, polarizabilities and covalent typing. Verify parameter assignment with `analyze GP` only. Check the sum of printed monopoles against the declared fragment charge using its explicit printed-rounding envelope. Asp/Glu fragments are −1; Asn/Trp fragments are neutral. Add La(III) or Dy(III) at the exact deposited ion position, generating **eight paired fixed configurations** with charges+2/+3. Physical La1/Dy6 labels describe potential future QM inputs; the classical model has no spin degree of freedom.

This declaration permits zero-energy preparation only: at most four `pdbxyz` and four parameter-only `analyze GP` invocations. No new reference/model energy, force, implicit solvent, relaxation or dynamics. Native protein donor parameters remain distinct from SI acetate/acetamide; no small-molecule POLPAIR, multipole or polarizability is silently transferred.

These are finite reference-qualification candidates, not validated models. A subsequent declaration must select the actual metal/protein Hamiltonian and reference recipe and preserve identical geometries. First test the native protein donor context before considering a separately defined parameter modification. Existing single-water and hydrated mixed-donor controls provide complementary evidence but cannot establish protein transfer alone.
