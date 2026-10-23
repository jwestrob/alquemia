# Protein donor fragments prepared; scoring remains unrun

The corrected V2 preparation completes four capped protein fragments and eight paired La/Dy geometries. It made four native `pdbxyz` and four parameter-only `analyze GP` invocations, serially on one CPU thread, in 5.600381 seconds. **No energy, force, optimization or dynamics calculation occurred.**

| Fragment | Atoms without metal | Fragment charge | La/Dy complex charge | Protein donor type |
| --- | ---: | ---: | ---: | ---: |
| ACE–Asp3–NME | 24 | −1 | +2 | 143 |
| ACE–Glu11–NME | 27 | −1 | +2 | 161 |
| ACE–Asn5–NME | 26 | 0 | +3 | 153 |
| ACE–Trp9–NME, backbone oxygen | 36 | 0 | +3 | 11 |

## Preserved geometry and repaired cap input

Every central-residue atom and the metal retain their exact deposited 9CEQ model-1 coordinates. The ACE cap retains the preceding CA/C/O. The NME cap retains the following N/CA/amide-H and, in V2, its deposited Cα hydrogen as explicit methyl `H1`. Five remaining cap hydrogens per fragment are constructed by native Tinker; no central atom is constructed or moved.

V1 remains preserved. Its native input supplied NME amide `H` but omitted methyl `H1`; `pdbxyz.f` uses `H` as the fallback for absent `H1`, producing distinct native atoms at the same coordinates. V1's uniqueness gate rejected that structure before energies. V2 supplies the separate deposited Cα hydrogen explicitly. The amide H and methyl H1 now have distinct native types 231/233 and the correct N/C parents. No force-field parameter or native source changed.

## Actual checks

- All retained atoms join to exact source coordinates; all coordinates are unique.
- Native covalent graphs are reciprocal and have the expected type valences.
- Each constructed atom is a cap hydrogen with one carbon parent; its native bond length is 1.10 Å within the declared 1e−6 Å serialization allowance.
- Native parameter output has no undefined parameters. Printed monopole sums are exactly −1.00000, −1.00000, 0.00000 and 0.00000 for the four fragments.
- The six named donor sites have the same printed local permanent multipoles, local-frame convention, polarizability and Thole parameters as their exact atoms in the complete protein preparation. This checks local parameter identity, not induced-dipole or environmental equivalence.
- All four La/Dy pairs have identical nonmetal coordinates, symbols, native types and bonds. Metal atoms have no covalent bonds.

The metal types 900/901 are reserved identifiers; a complete executable protein metal Hamiltonian has not been selected. Physical La1/Dy6 labels are retained for later reference preparation; this classical model has no electronic-spin degree of freedom. These are isolated capped fixtures, not hydrated binding sites or equilibrium structures. Native protein donor parameters remain different from the small-molecule SI parameters; copying its donor corrections would require separate qualification.

## Exact evidence

- [V2 panel](../../../workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/peptide_fragments_v2/PANEL.json): `42db5a66e95edf2af5a5fd1eb77c779dc18cac4bcfeaadde39d6beba9fe71b67`.
- [Saved parameter comparison](../../../workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/peptide_fragments_v2/PARAMETER_COMPARISON.json): `0a48c11834b0440d5d420b85615fba5f686473e5441a3dba3eb78463cc191f00`.
- [Preparation source](PREPARE_FRAGMENTS_V2.py): `8189614a92079fc56f3ed591cf7504b0c7c9548c914fb036d301a988adfbd8d4`; [specification](PEPTIDE_FRAGMENT_SPEC_V2.md): `921d26632fbf9c57ce6286f96c65bab7926f1662ba7273c588ac86fcb5fb87ac`.
- Preserved V1 panel: `131f46e2c505d1a90638452091ba777c936c299f7fa60e19674cd28813b72953`.
- Protein parameter file: `eda92c3683846e340c9b07991ed15824c00023d3f3b64effee08d1b53583b6ca`; 9CEQ source: `956ed5bf63bc4e4c05c223d9e34912a0f57c0d1b757a1648c7c385cc4fab967b`.

All native input/output hashes, atom joins, declared eight-cell denominator and subprocess costs are retained in the panel. Root review precedes fragment scoring.
