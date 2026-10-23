# Exact periodic parameter mapping inventory

Read the actual paired-leg coordinates, atom provenance and pinned AMOEBA-BIO-2009 parameter file. Inventory every represented nonmetal type/class and its unmodified atom, van der Waals, polarization and permanent-multipole definitions. Compare the small-molecule SI water/donor records explicitly. This is metadata/source preparation: zero energy, force, dynamics or fitting calls, no candidate freeze and no reserved validation outcomes.

Retain protein/water/ion bonded terms and permanent multipoles/polarizabilities exactly. Reserve metal types **and classes**900/901, requiring absence from the entire protein parameter namespace. Map VDWPR by class and POLPAIR by type. No spatial contact cutoff selects recipients: all occurrences of a chemical family receive its hypothetical pair override, including non-coordinating backbone/terminal carbonyls. All remaining represented classes retain explicitly inventoried default mixing, not invented donor overrides.

Do not copy the SI key into the protein namespace. In particular its acetate type403 collides with protein water H403, and its amide types104/105 already mean Trp aromatic atoms in this protein. Published numerical La pair values are source evidence, not a selected periodic model. Any future candidate must supply a separate source pin and transfer decision. This inventory writes no runnable periodic key or model parameters.

Water is not an exact full-Hamiltonian renumbering: SI Urey–Bradley distance1.5337Å differs from the retained BIO value1.5537Å. Both legs keep BIO water unchanged; qualification of a small-molecule water candidate must account for that transfer. Donor polarizabilities and permanent multipoles also differ. Semantic mapping does not qualify transfer.
