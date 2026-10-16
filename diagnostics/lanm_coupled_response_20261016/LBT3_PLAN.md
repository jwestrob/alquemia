# Whole-scaffold monomer control: LBT3 La/Dy, finite first scout

Declared16October2026 before new molecular evaluations. Root owns execution.
This follows Jacob's recommendation to begin with monomers with measured
preferences. See MONOMER_TARGETS.md. The active goal remains a LanM profiler;
LBT3 is a small opposite-preference control, not a substitute second LanM.

## Question

Does including the entire real17-residue LBT3 scaffold change the La/Dy response
and its dependence on light/heavy-conditioned source geometry compared with the
old donor-fragment representation? This is a cheap whole-system electronic
diagnostic preceding a covalent-preserving response model. Both source structures
are mandatory; no favorable-source selection, score fitting or new full series.

Use consumed solution-NMR model1 from La7CCO and Lu7CCN. The old14-ion carved
DFT/aquo experiment failed and remains closed. This is a new finite representation
comparison using the existing native MACE-OMOL model, not its restart.

## Source preparation and state

The deposited model1 sources each have252 peptide atoms plus one metal, no
deposited waters or spectators. Both lack terminal OXT and two terminal amine H.
Use the existing deterministic C-terminal OXT construction, ff19SB/OpenMM pH6
missing-H preparation with seed20260918, and preserve deposited coordinates.
Its bounded internal minimization of new H is real preparation computation,
counted separately. Do not call the preparation zero-energy work. Report actual
new-H work if available; otherwise state that internal force-call count is unknown.
No missing internal residues or donor heavy atoms may be invented. Existing H
coordinates remain unchanged; audit their bond lengths rather than automatically
normalizing them. Retain topology, atom maps and native protein mechanics.

Five deprotonated sidechains D3/D7/D13/E11/E14 plus charged N/C termini give the
intact peptide charge-5 and trivalent holo charge-2 (256 expected atoms). The
zwitterionic terminal choice is our declared chemical assumption; the main
paper's MD holo charge-2 is compatible with it but does not establish the
experimental terminal chemistry. Supplementary FigS5 separately says peptide
charge-4, an unresolved primary-source discrepancy retained here. This scout is
conditional on its explicit state, not a microscopic titration or experimentally
verified terminal assignment. No explicit waters are added.

Construct the matched48-atom donor-fragment control with the archived exact
recipe: D3/N5/D7/E11/E14 sidechains plus W9 carbonyl; ligand-4, holo-1. Use the
new whole source's shared physical coordinates exactly, and the legacy mapped
caps. Do not change the old archived clusters. The omitted D13 charge is a real
representation difference, not an accidental mismatch to suppress.

Each source/representation is evaluated at IDENTICAL coordinates for La and Dy,
with target oxidation state+3, multiplicities1 and6 respectively. Audit all-electron
parity, actual native model charge/spin batch, element support and source mapping.
Multiplicity is a declared physical hypothesis, not validated f-electron physics.

## Eight calls and explicit interpretation

Two sources x two representations x two metals = eight native MACE vacuum
singlepoints, one warm calculator. Reuse the exact existing100M OMOL checkpoint,
native charge/spin embeddings, float64, omol head and no extra dispersion.
No optimization, GFN2, ORCA or new solvent calculation in this scout.

For each representation, retain the complete two-metal/two-source energy matrix.
Report the change in Dy-minus-La contrast upon including the scaffold and the
difference between the two source contrasts. These balanced differences cancel
the common elemental offsets. Raw unlike-metal energies are not affinities.
The compatible two-source pool may be summarized mathematically for each
representation, with individual-source results retained; counts are not populations.

Where exact native model/state compatibility can be established, reuse the three
consumed Mex EF1-3 MACE endpoints for the ordered exchange vector
`R_i = (E_Dy-E_La)_Mex,i - (E_Dy-E_La)_LBT3`.
Positive R is the qualitative direction expected for greater relative light-ion
preference in Mex. This combines a protein-level qualitative expectation with
conditional site calculations; there are no independent EF-site labels and no
site summation. Mex/LBT3 assays differ in conditions and observable. No numerical
experimental exchange target, Kd, free energy, class threshold or success rate.

Missing or nonfinite cells remain unavailable; no substituted metal, artificial
spin, zero solvent energy presented as a completed composite, or favorable fallback.
The scout can reveal a useful scaffold/structure effect or reject this particular
representation. Either outcome guides the next finite mechanistic experiment.

## Resources, provenance and next step

Prepare geometry on one standard Slurm CPU; evaluate on one allocated GPU with
current scheduler-granted CPU/RAM resources and one warm worker. No memory/time
override, login numerical work, new environment installation or PQQ change.
Pin all sources/preparation/model/runtime/code, record actual calls and costs,
and arm terminal wakes. Tests of parsing/maps and executed science stay distinct.

If the prepared whole-system model is executable, inspect actual response and
source dependence before a separately declared bounded bond-preserving search.
New solvent work must use a qualified Hamiltonian: the corrected-La native GFN2
registry DOES NOT qualify Dy. Never import the old incorrect native solvent
energies or infer Dy support from a successful La parameter audit.
