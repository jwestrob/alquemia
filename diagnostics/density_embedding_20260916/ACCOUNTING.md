# Energy accounting before execution

The potential probes are actual permanent-charge coordinates, in their recorded
order. ORCA vpot returns the quantum nuclear-plus-electronic potential from the
saved vacuum density; La uses its native ECP effective core. Existing real ESP
checks verify that this convention agrees with the endpoint charge convention.
Both MBIS and density Coulomb sums use Bohr and atomic units before conversion
by 627.509474 kcal/mol per Hartree. The old TABI constant is retained separately.

External point charges contribute to the electronic Hamiltonian and nuclear
repulsion. Explicit `DoEQ false` omits external-external self interaction.
Therefore `E_embedded - E_vacuum - C_density` is the electronic response to
that fixed field. Coordinates, charge, basis, native D4 and gCP are unchanged.
A positive response above 0.05 kcal/mol triggers an accounting/convergence
review: variational relaxation should lower energy, with this small allowance
for the existing SCF/integration settings. This is not an accuracy threshold.

This diagnostic has no solvent energy or affinity classification. Per-residue
coupling differences locate charge-representation errors; they do not establish
which residues require quantum treatment or authorize changing the carve.
