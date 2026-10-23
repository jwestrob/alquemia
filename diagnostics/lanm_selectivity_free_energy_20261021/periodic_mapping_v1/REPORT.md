# Exact small-molecule → periodic protein mapping

**Mapping is concrete; transfer is unqualified.** All four paired-leg boxes and the full protein parameter namespace were inspected. The record covers85 represented nonmetal types across25 classes, every peptide atom context and all atom/VDW/polarization/multipole source records. It writes no executable key, selects no metal candidate, and performs no native, energy, force or dynamics calls. Metadata preparation took0.636s on one CPU thread.

Machine record: `workspaces/lanm_selectivity_free_energy_20261021/periodic_mapping_v1/prepared_v1/MAPPING.json`, SHA256 `edaa4cc6cbeb8b82670e0314423432667102640eb4cf0df84e2eb0bc300f34e7`. Protein parameter SHA256 `eda92c3683846e340c9b07991ed15824c00023d3f3b64effee08d1b53583b6ca`. Source/config and all four XYZ/atom-map pins are retained in the record. A first metadata attempt expected the wrong parent status string and stopped before output; `preparation_attempt_v1` preserves it. The corrected gate requires the exact existing `geometry_prepared_hamiltonian_unassigned` status and unselected-metal flag.

## Complete recipients of the chemically relevant pair terms

VDWPR addresses **classes**. POLPAIR addresses **types**. Proposed metal type/class namespaces900 for La and901 for Dy are unused throughout the entire BIO parameter file, not merely absent from the current structure.

| Chemical family | SI source type / class | Actual periodic type / class | Bound occurrences | Preserved protein polarizability, Å³ |
| --- | --- | --- | ---: | --- |
| Water O | 36 /34 | 402 /73 | 3894 | 0.837 |
| Water H | 37 /35 | 403 /74 | 7788 | 0.496 |
| Asp carboxylate O | 403 /403 | 143 /31 | 6 | 0.837 |
| Glu carboxylate O | 403 /403 | 161 /31 | 4 | 0.837 |
| Gly backbone carbonyl O | 105 /55 | 5 /5 | 2 | 0.837 |
| Other backbone/terminal carbonyl O | 105 /55 | 11 /5 | 15 | 0.837 |
| Asn side-chain carbonyl O | 105 /55 | 153 /5 | 2 | 0.834 |
| Gly carbonyl C | 104 /53 | 3 /3 | 2 | 1.334 |
| Other backbone/terminal carbonyl C | 104 /53 | 9 /3 | 15 | 1.334 |
| Asn side-chain carbonyl C | 104 /53 | 152 /3 | 2 | 1.334 |

The final Ala17 carbonyl belongs to the same backbone families; the terminal NH2 is type228/class1. There is no terminal COO− in this construct. ACE types226/227 are absent from the full paired boxes; they belong to the separate capped-fragment tests and must not be silently inserted here. Class-wide overrides necessarily reach all19 represented amide oxygens/carbons and all10 carboxylate oxygens, independent of current distance or which donor happens to bind.

Protein carboxylate C types142/160 are **class30**, not amide class3. No published La carboxylate-carbon override exists; do not invent one by analogy. Thr4 hydroxyl O42/class10 (alpha0.834), Trp9 N98/class24 (1.073), peptide/Asn/terminal amide N classes1, and N-terminal N234/class41 remain native BIO types with no borrowed oxygen correction. All hydrogens and aromatic/aliphatic atoms remain represented through their original parameters.

Na405/class76 retains charge+1 and alpha0.120; Cl414/class85 retains charge−1 and alpha4.000. Counts are Na8/Cl7 bound and Na7/Cl10 aqueous. Their metal cross-interactions must be present through declared native mixing unless a separately qualified override exists. They are not water or donor types. Salt/ion-association physics remains unqualified.

## Exact proposed remapping of the published La pair records

These are source-preserving **mapping proposals**, not activated parameter lines:

| Original SI operator | Proposed periodic destination | Original values retained if that model is selected |
| --- | --- | --- |
| VDWPR55,500 | VDWPR5,900 | radius3.3000, epsilon0.3019 |
| POLPAIR105,500 | POLPAIR5,900;11,900;153,900 | damping0.2987 each |
| VDWPR53,500 | VDWPR3,900 | radius3.9541, epsilon0.8295 |
| VDWPR403,500 | VDWPR31,900 | radius3.6318, epsilon0.2280 |
| POLPAIR403,500 | POLPAIR143,900;161,900 | damping0.1518 each |

The amide-carbon VDWPR is essential to an exact transfer of that candidate; omitting it changes the model. There is no published La-specific water pair override in this SI key; its ordinary La/water mixing must remain explicit. A future qualified water correction would target class73 with metal class900/901, and type402 with metal type900/901. A future carboxylate correction targets class31 and types143/161; amide oxygen targets class5 and types5/11/153. Future Dy numerical values are deliberately null. Never copy La corrections onto Dy by default.

**Namespace collisions prohibit copying the SI key.** Small acetate403 becomes water H403 in BIO; small acetate C402 becomes water O402; small amide O105 and C104 become Trp aromatic C/H; small water O36 is a protein Ser oxygen type. Each parameter reference must be translated by chemical role and operator, including polarization group/axis references if a new atom definition is ever introduced.

## Physical differences that the mapping must preserve

Small acetate O403 alpha1.931 differs from protein O143/161 alpha0.837; SI amide O105 alpha0.914 differs from protein carbonyl0.837/0.834; SI amide C104 alpha2.065 differs from protein1.334. Their permanent multipoles and local frames also differ. The inventory preserves every native multipole record; it does not replace protein values with SI values to make the transfer appear equivalent.

Water vdW, permanent multipoles, alpha/Thole, bond and angle terms match after renumbering, but **Urey–Bradley equilibrium distance does not**: SI `35 34 35 -7.60 1.5337` versus BIO `74 73 74 -7.60 1.5537`. Both periodic legs retain BIO1.5537Å unchanged. This is a water-model transfer difference, even though the internal term cancels from a rigid translation of a fixed water. It can matter for flexible-water distributions. The negative force constant is an existing parameter, not a new correction.

## Minimal future binding design

1. Bind one explicitly selected paired candidate, this mapping record and the exact BIO file. Preserve every nonmetal topology/type/mass and all BIO bonded, multipole, polarization and global mixing records byte for byte. Do not include the SI atom/bond/multipole blocks.
2. Append only new metal definitions900/901 and the explicitly selected mapped pair records. Keep charge+3 and the previously declared common metal kinetic mass convention in both legs; intrinsic metal polarization/VDW values require the candidate pin. No polarization-group bond connects the metal to a chosen donor or water.
3. Enumerate every represented class pair with each metal, marking native mixing or one explicit override; reject duplicate/conflicting records. The present classes are1,2,3,4,5,6,7,8,9,10,11,17,18,22,23,24,25,30,31,41,42,73,74,76,85. Regenerate the existing alchemical path in this namespace, including water H, Na and Cl; small-system path success does not bind these pairs automatically.
4. Read back resolved native atom/class, pair and multipole assignments before energies; require unchanged nonmetal definitions and exactly the declared coverage. Then require mapped CPU/GPU energy/force/path checks and the declared physical water/protein transfer gates before any periodic pilot. These are future gates, not additional executions requested by this inventory.

No mapping operation resolves the water48 reference uncertainty, the failed fixed-La control, donor transfer or missing bound hydration population. The whole boxes remain geometry/chemistry preparations without a qualified periodic Hamiltonian.
