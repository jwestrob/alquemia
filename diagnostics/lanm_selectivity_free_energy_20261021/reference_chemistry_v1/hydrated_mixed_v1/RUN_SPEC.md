# Hydrated donor competition: finite scoring pilot

Root requested this separate pilot after the physical preparation. It contains
exactly six native central cells: monodentate acetate with seven waters,
bidentate acetate with six waters, and acetamide with seven waters, each La/Dy.
The existing two aquo8 native references are reused after successful admission;
they are not runnable cells in this pilot. No explicit-f cells are included.

Use the frozen REFERENCE.py executor, its exact ECP PBE0-D4 recipe, basis assets,
DefGrid3, VeryTightSCF, PModel, MaxIter300 and native admission parser. Four MPI
ranks per singleton, one thread per rank, standard_expand and scheduler default
memory/time. MaxCore comes from the actual memory grant as before. No fallback,
new optimization, parameter fit or automatic extension is included.

The separate MACE config contains exactly the18 prepared mixed configurations,
using root's existing pinned RUN.py, native helper, model, runtime and physical
La1/Dy6 inputs. This excludes the two already-computed aquo8 states. The root
runner reads its finite count from config; no modification is needed.

Preserve the first preparation and prepared_v2. Bind the latter's exact XYZ
bytes and topology. The native pilot has effective singlets and physical-spin
metadata; the MACE inputs retain explicit physical multiplicities. This is a
method distinction, not equivalent electronic-state treatment.

The central contrast for each donor is E(Dy mixed)-E(Dy aquo8)-E(La mixed)+E(La
aquo8). Free ligand and displaced-water energies cancel in this balanced cycle.
The original native collector is specialized to the initial56 cells and must
not collect this pilot. Root's new saved comparison must admit both provenance
sets explicitly, confirm each paired XYZ, and require all four terms.

The frozen parser's radial_metal_force is only the metal force projected onto
the stored outward ligand axis. The prepared distortion moves the entire ligand;
its generalized force is the sum of all ligand forces dotted with that axis.
Do not label the metal projection as the derivative of this ligand displacement.
Keep full force arrays for both evaluations.

This tests many-body water/donor competition at fixed stress configurations.
Neither a sign match nor native SCF termination establishes parameter accuracy,
equilibrium hydration, protein selectivity, buffer speciation or reference truth.
Scoring remains pending root's review/release; this binding makes zero calls.
