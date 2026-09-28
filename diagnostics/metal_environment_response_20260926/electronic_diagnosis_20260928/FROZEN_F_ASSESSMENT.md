# Frozen-4f Dy: a defensible separate reference candidate, not an SCF rescue

28 September 2026. **Pursue one contained energy/analytic-force qualification of the existing ECP55 model.** The archived implementation demonstrably converges, but its force consistency and LanM selectivity remain unqualified. No molecular call, production change, or reserved-outcome inspection was performed for this assessment.

The latest DY_PBE0_CAPABILITY_RESULT.md records failure of the isolated explicit-4f PBE0 sextet as well as native r2SCAN-3c. This makes it reasonable to investigate a deliberately simpler electronic representation. It does not prove that 4f electrons are chemically irrelevant or that the converged approximation is accurate.

## Exact model and state distinction

The installed legacy recipe is PBE0-D4/def2-TZVP on H/C/N/O plus **lcecp-1-TZVP on Dy**, its matching 55-electron ECP, explicit Dy AutoAux Coulomb-fitting basis, light-atom def2/J and RIJCOSX. QZVP/QZVPP versions retain the same ECP; they are existing sensitivity assets, not requests for a new sweep. The native 3c recipe, its gCP term and its reference offsets do not carry over.

The actual Dy `NewECP` says `N_core 55`, `lmax g`. Its core comprises the 46 inner electrons plus a prescribed **4f9** occupation. Neutral Dy therefore has eleven explicitly treated electrons; formal Dy(III) contributes eight. The molecule's electron count must be recomputed from the actual atoms, charge and ECP, not obtained by copying the sextet input.

Physical Dy(III) remains an open-shell 4f9 ion with S=5/2. The ECP55 restricted valence calculation uses effective multiplicity1 because those open-shell electrons are omitted from the variational space. This is a **spin-free closed-shell valence model**, not a physical singlet Dy state, not the explicit-f sextet solution, and not a spin-orbit-resolved calculation. Store both the physical state hypothesis and backend effective state. It cannot predict 4f occupation changes, multiplet splittings, magnetic anisotropy or explicit environment-induced 4f response.

The valence electrons and ligands still respond self-consistently to geometry and embedding. Freezing 4f therefore does not freeze the whole ion or reduce the calculation to a point charge. It deliberately limits which electronic response can occur; any success must be described within that approximation.

## Primary-source basis, and its scope

Lukanowski and Weigend developed occupation-specific large-core basis families and assessed geometry, frequencies and exchange energies for120 molecules against basis-limit and scalar-relativistic all-electron calculations. Their accompanying assets include the matching potentials. This supports a conventional energy-surface approximation; it does not validate protein metal selectivity. [Primary paper](https://pubs.rsc.org/en/content/articlehtml/2026/cp/d5cp04944j), DOI10.1039/D5CP04944J. Direct HTML retrieval returned403 in this audit; the publisher's indexed abstract/results and primary PubMed abstract were accessible.

A separate lanthanide chloride-cluster study used occupation-specific large-core models, optimized structures with several functionals, and compared selected PBE0 results with all-electron scalar-relativistic treatments. This is relevant precedent for geometry-sensitive ionic coordination, not evidence for LanM affinities. [Primary cluster study](https://pubs.rsc.org/en/content/articlehtml/2025/cp/d4cp04057k).

Our scientific inference: for fixed Dy(III), fixed proton/water inventory and oxygen-donor coordination, this is a plausible way to retain ionic size, valence polarization and ligand strain while removing the difficult explicit-f optimization. The importance of missing f-shell response in small La/Dy selectivity differences remains unknown. Ionic-looking coordination is a rationale to test the approximation, not a quantitative error bound.

## What the real archive establishes

`FROZEN_F_INVENTORY.json` pins the actual basis/ECP/AuxJ assets, reviewed protocol and91 large-core stage outputs in the consumed August calibration shards. This count includes stages/starts and is **not91 independent molecules or accepted results**.

-46 outputs contain final energy, normal termination, and a native stable-wavefunction statement.
-None of those91 outputs prints a Cartesian gradient; no `.engrad` exists anywhere in the inspected calibration shards. The inputs are energy calculations. This is no archived analytic-force qualification.
-The q155 B10-medoid PModel output explicitly reports ECP55, charge−1,654electrons and multiplicity1, then stable convergence. The matched ECP28 recipe has681electrons: the27-electron difference is exactly55−28, not a changed formal charge.
-Seven task identities have both completed HCore/PModel starts. Only two meet the **original**1e−7Eh energy-agreement threshold; five miss it, by differences up to4.79669e−7Eh. These differences are small on a chemical energy scale, but the historical gate stays failed. No occupied-subspace-overlap qualification was established in this audit.
-All eight explicit-4f PBE0 molecular outputs previously inspected lacked convergence, so no completed matched frozen-f/explicit-f correction exists here. Do not compare their raw totals or fill that gap with the frozen-f answer.

| Archived matched-start identity | HCore−PModel, Eh | Old energy gate |
|---|---:|---|
| B10 extreme q155 r2SCAN TZ |−2.44045e−8|pass, energy only|
| B10 extreme q86 PBE0 TZ |+1.99900e−7|fail|
| M9 extreme q86 PBE0 TZ |−2.97236e−7|fail|
| M9 medoid q86 PBE0 TZ |−3.13722e−7|fail|
| B10 medoid q86 PBE0 TZ |−4.79669e−7|fail|
| M9 medoid q155 whole-residue PBE0 TZ |−5.41877e−9|pass, energy only|
| M9 medoid q155 PBE0 TZ |−2.47082e−7|fail|

These demonstrate reusable executable model assets, not completed historical calibration, within-series discrimination, or proof that its force field is suitable for protein relaxation. No reserved SpyCI-LAMBS results were needed or opened.

## Small next test and interpretation

Use the already consumed real50-atom Hans EF3 source, retaining its exact geometry, proton/water inventory and physical Dy(III) hypothesis. Declare a new frozen-f protocol identity and the effective electron/spin convention before evaluating anything. A bounded analytic-gradient endpoint plus selected real donor-direction finite-difference checks can establish whether this executable provides a numerically consistent local surface. Keep independent-start agreement visible; do not substitute easier acceptance thresholds into the historical experiment. Actual execution belongs to root and is not launched by this report.

If that passes, a matched source-derived small displacement is more informative than a favorable absolute energy. It can test the restoring response that failed to qualify under explicit-f methods. Extension to embedded A/B then tests environment response, separately from local chemistry. Neither requires new whole-protein optimization or repeats of the old long-DIIS ladder.

For eventual La/Dy preference, choose compatible element-specific reference treatments and a balanced hydrated-ion exchange cycle, or a same-stoichiometry protein–protein exchange contrast where element offsets cancel. A Dy-only successful endpoint is not a La/Dy discriminator. Never subtract an ECP55 total from a native-ECP28 or native-3c reference offset; electronic-core constants are different. Full protein relaxation still requires a coherent solvent-containing energy and force model, and state/occupancy populations remain separate from a minimum-energy candidate.

**Decision:** investigate frozen-f as an explicitly approximate, separately qualified reference. Preserve explicit-f failures and all historical gate failures. Do not present it as a solved Dy state, adopt it in production, or infer biological affinity until the relevant comparisons actually exist.
