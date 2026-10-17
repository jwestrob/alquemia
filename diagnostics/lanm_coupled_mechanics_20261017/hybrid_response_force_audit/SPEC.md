# Saved same-coordinate force changes at four hybrid response starts

## Frozen question and denominator

Before inspecting force differences, compare the existing whole-LBT Hamiltonian
and baseline local-electronic/MM hybrid Hamiltonian at each of the four exact
selected starts admitted by response preflight 13015: 7CCO/La, 7CCO/Dy, 7CCN/La,
7CCN/Dy. Selection is inherited unchanged from the complete common pool. Retain
all four slots; do not select atoms, sources or a different geometry for a
favorable result. This is a saved-data diagnostic with no molecular calls,
optimization, geometry modification or changed response settings.

The question is which mapped coordinates and energy components account for the
change in restoring forces when the Hamiltonian changes at fixed coordinates.
It does not identify a uniquely correct force model, an equilibrium structure,
or a binding affinity.

## Required evidence and signs

Replay the actual response preflight and its complete parent/certified recovery
through the frozen response `prepare()` interface. Use exactly its selected
candidate and archived POINT for each source/metal. Admit the corresponding
whole-model cell from the same Module 2 candidate through existing raw-component
helpers, requiring a saved full 256-atom gradient. Require exact physical atom
order and coordinates and raw component-gradient closure on both sides.

If a required whole gradient is absent, retain an unavailable slot with its
actual artifact/reason. Do not synthesize a force from energy-only receipts,
reuse the other metal's force or perform a new calculation.

All vectors are physical **forces** `F = -gradient`, in kcal/mol/angstrom.
MACE native force is eV/angstrom; native ALPB gradient is Hartree/bohr. Existing
qualified conversion constants and raw admission arithmetic are retained. Report
`delta F = F_hybrid - F_whole`. Whole-model field and MM contributions are
declared zeros by that model. Hybrid components are local MACE, embedded-minus-
isolated GFN2, the declared MM ledger, and the same global fixed-charge ALPB.
Their force changes must sum to the total. The same-coordinate ALPB term must
replay from the exact compatible saved receipt rather than be recomputed.

## Predetermined descriptive outputs

Write all 256 atom rows per available state, with IDs, residue, element, both
total force vectors, their difference, norms and all component vectors. Sort a
separate top-20 view by decreasing delta-force norm, breaking ties by physical
atom index. A top-10 residue view uses decreasing delta-force L2 norm, with ties
broken by first physical atom index. Include all residues and component/category summaries; the top
view never determines the denominator or interpretation.

Labels come only from the actual physical partition and cap map. Distinguish the
target metal, retained electronic real atoms, retained and omitted cap-boundary
anchors, and other MM physical atoms. Boundary flags are explicit and may
overlap the exclusive electronic/MM membership. Synthetic caps are not physical
atoms and are not additional rows. No radial distance cutoff assigns membership.

For each full population, residue and exclusive partition category, report atom
count, vector sum, sum of atom norms, RMS atom-vector norm, maximum and fixed
0/25/50/75/95/100 percentiles for whole, hybrid and difference forces. Also report
component vector sums and norms. These describe vectors and distributions;
they are not additive residue energies or independent causal percentages.

Keep the finite full-force residuals and prior convergence flags visible.
Numerical agreement of energy derivatives is separate from physical accuracy;
physical La1/Dy6 and effective frozen-f GFN2 assumptions remain unchanged.
No new physical acceptance tolerance, bad-bond label, or force-selected motion
follows from this audit. Source-only edits and AST parsing occur on login;
all artifact replay and force arithmetic run under root's Slurm allocation.

The two frozen helpers multiply/divide the same native ALPB array in different
floating-point order. Require exact common raw ALPB receipt and constants; retain
each helper's admitted arithmetic and report any resulting last-bit component
difference and component-sum residual. Do not introduce a force tolerance or
interpret that conversion roundoff as changed solvent physics.
