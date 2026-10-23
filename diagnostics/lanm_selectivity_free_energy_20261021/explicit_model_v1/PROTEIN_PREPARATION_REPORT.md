# Amidated LBT3: paired explicit-water geometry prepared

The first 9CEQ model has been converted into a native AMOEBA biopolymer topology and two matched La/Dy explicit-water coordinate boxes. **No protein Hamiltonian has been selected and no protein energies, forces, minimization or dynamics were evaluated.**

Native `pdbxyz` retained all **257 deposited peptide atoms**, including every hydrogen. The only input naming change was the terminal NH2 cap's HN1/HN2 → H1/H2 convention. Every source atom has an exact coordinate join to the native topology. The sequence is PHE–ILE–ASP–THR–ASN–ASN–ASP–GLY–TRP–ILE–GLU–GLY–ASP–GLU–LEU–LEU–ALA–NH2. Native parameter-only analysis assigns all atoms and prints monopoles summing exactly to **−4.00000 e**.

Each paired box contains **11,955 atoms**: the 257-atom peptide, one +3 metal, 3,894 intact waters, eight Na+ and seven Cl−. The cubic box is 49.324 Å wide; solute padding is at least 15.0 Å. Seven salt pairs correspond to approximately 0.09687 M NaCl, with one extra sodium neutralizing the −1 complex. Both metals have identical coordinates and atom order; only the metal name/type differs.

Water positions come from the native example's TIP3P packing, assigned AMOEBA water types. They are **unequilibrated initial coordinates**, not an AMOEBA equilibrium sample. No MES buffer is represented. The water/ion packing and a single deposited conformation do not establish experimental state populations.

## The remaining parameter-transfer gate

The SI's small-molecule atom numbers must not be copied into the protein file. The native protein water types 402/403 have the same water parameters as the SI's 36/37 after renumbering, but the donors are materially different:

| Function | La SI small-molecule parameters | Native protein parameters |
| --- | --- | --- |
| Amide oxygen | Type105/class55; polarizability0.914 Å³ | Backbone5 or11/class5:0.837; Asn153/class5:0.834 |
| Amide carbon | Type104/class53; polarizability2.065 Å³ | Backbone3 or9/class3:1.334 |
| Carboxylate oxygen | Type403/class403; polarizability1.931 Å³ | Asp143/Glu161/class31:0.837 |
| Permanent donor multipoles | Fitted small-molecule contexts | Different residue/context-specific moments |

A La `POLPAIR`/`VDWPR` mapping would therefore be a **physical transfer modification**, not just relabeling. It must be declared and qualified on protein-like donor contexts. No executable protein metal key or fabricated Dy correction was selected. The paired cluster gate and mapped-donor reference tests precede protein scoring or sampling. The current boxes use reserved metal types900/901 and remain preparation artifacts.

## Execution and preserved failure

The first preparation completed native topology conversion, then failed while a Python parser attempted to tokenize an unrelated parameter-file comment containing an unmatched quote. That result/source/output remain in `protein_preparation_v1`. Version2 reads only `atom` definition lines and **reuses the already completed topology**; it then runs native `analyze protein.xyz GP` for parameters/system information only. The native energy and force options are off.

Version2 preparation took30.90seconds on one CPU thread; its parameter-only native invocation took0.52seconds. The source/receipts explicitly record zero energy/force calls and zero optimization/dynamics steps. This is a contained parser recovery with no changed parameter, geometry or scientific acceptance threshold.

Evidence under `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/`:

- `protein_preparation_v2/RESULT.json`: `bf73c677294f962701d1eff7e2f63a77898e495b1a71ae2be053771363411fb4`.
- `protein_preparation_v2/PARAMETER_REVIEW.json` records the exact native monopole sum and output pin.
- La box: `0599921655fe2ed3eaf604f9ac1289d7d4862c1e36cac1efe0eb8ccf03bb666a`; Dy box: `d99ca71c929a0992b0af166277cd7d5a53634b485f499d0e0a0784e8dcfc9015`.
- Preserved original failure: `protein_preparation_v1/RESULT.json`, `ff567c1eede1a1310e7d878f9725d48f0c28159fc3d529cfe82471a31d82c033`.

Native topology behavior is documented by the [Tinker utility manual](https://dasher.wustl.edu/tinker/distribution/doc/sphinx/tinker/_build/html/text/analysis-programs.html); exact source, binaries and parameter file are pinned in the saved preparation result. Target provenance and terminal-state qualifications remain in the experimental target report.
