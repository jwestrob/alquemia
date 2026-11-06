# Test 1: does a distant water alter the metal-swap energy?

Frozen 6 November 2026, before new predictions. Jacob approved the proposed
one-CPU-hour diagnostic with “sure.” This authorizes Test 1 only. The earlier
screen remains recorded; its missing reference comparison did not establish
Tauri failure. No DFT, MD, new experimental target, fitting or production change.

Question: does adding the SAME distant neutral water to both endpoints change
the energy difference of a same-composition La/Dy swap? This directly tests
cancellation of the earlier single-complex spectator effect. Passing establishes
this limited consistency property, not correct selectivity or solvent sampling.

Endpoints use existing exact internal coordinates:
A = LaAcW7 + DyW8; B = DyAcW7 + LaW8.
Translate the acetate metal to(-20,0,0) A and aquo metal to(+20,0,0) A. There are
54 atoms, one La and one Dy in either endpoint, total charge+5, multiplicity6
(five unpaired electrons). A/B positions, atom order, elemental inventory and
global state match; only the two metal identities exchange. Formal fragment
charges are+2/+3. Do not impose fragment-charge constraints or change embeddings.
This charged, isolated prototype does not reproduce the later neutral periodic
ligand box; conclusions are limited to the tested prototype.

Take the first actual OHH triplet from archived LaW8. Keep its geometry; place
its oxygen at(0,d,0) A for d=100,200,400 A, using the original orientation and
a180-degree rigid rotation about y that reverses its dipole. Add the identical
triplet to A and B:57 atoms, same charge+5 and multiplicity6. No optimization.

Use unchanged pinned so3lr-tauri-l-2606 as primary, M-2606 as fixed capacity
sensitivity. Source70707f5, GLP2294f4f and existing Python3.12.0/JAX0.10.2 env.
Float64; x64 enabled before loading; highest matmul precision; head5; no PBC.
Use a5-A neighbor list, preserving native model cutoffs L5/M4.5 A. Primary LR
cutoff1000 A, damping2 A; fixed no-cutoff sensitivity with dampingNone. All
nonperiodic pairs fit inside the loader's1000-A complete pair list.

14 joint frames x2 models x2 cutoffs =56 predictions, plus four free-water
predictions for frozen-charge electrostatic estimates, and two real LaW8
replays against saved arrays to check the changed padding/intermediate outputs.
Maximum62 new ML energy/force calls. No selection of models/settings by result.
Native loader and standard make_so3lr, not Flash or metadata-limited ASE.

For each model/cutoff, calculate full-precision:
D0=E(B)-E(A), Dwater=E(B+water)-E(A+water), delta=Dwater-D0.
Save raw float64 arrays before formatting. Also record the change in the
metal-swap force difference and model partial-charge totals on the two complexes
and water. Partial charges are model diagnostics, not oxidation-state labels.

Primary criterion, declared before values: max|delta| over both water
orientations at200/400 A <=0.5 kcal/mol passes this limited check. The100-A
frames show approach with distance. The no-cutoff repeat tests whether the
conclusion depends on finite-cutoff shifts; M never substitutes for L.
Report inconclusive if cutoff sensitivity changes the primary conclusion or
ordinary remote interaction cannot be shown small. No threshold/pose tuning.

Estimate ordinary remote electrostatics using charges frozen from the bare
joint endpoints and isolated neutral water, Coulomb sum k*q_i*q_j/r_ij. This is
an explanatory estimate, not a correction to the scored energy. A magnitude
<=0.05 kcal/mol (one tenth of the proposed tolerance) supports treating the
remaining far-water effect as irrelevant coupling for this diagnostic.
Use the model's native Coulomb constant14.399645351950548 eV A/e^2. Record
energy-component changes and the old separate-cluster C_Ac as descriptive
context only; neither changes the0.5-kcal/mol criterion or validates accuracy.
Opposite water orientations additionally distinguish a decaying dipole response
from a persistent shift. No reference-energy accuracy claim follows.

Checks: exact translated parent coordinates; endpoint composition/positions;
charge/spin/head propagation; complete LR pair coverage; float64 finite outputs;
sum of predicted charges; native energy-component closure; saved LaW8 replay
within1e-6 eV and eV/A. Failures remain explicit, never replaced by zero.

Resources: Slurm standard_expand, root-owned,12 CPUs, no GPU. One allocated
CPU-hour total cap. An internal240-second timeout plus5-second kill grace leaves
headroom below3600 allocated CPU-seconds. No scheduler memory/time override.
Source/metadata preparation is local; scientific numbers are produced on Slurm.
Record actual costs and preserve attempts. No automatic Test 2 or MD afterward.
A failure here would establish non-cancellation in this charged prototype; it
would not by itself establish failure in a neutral periodic ligand box.
