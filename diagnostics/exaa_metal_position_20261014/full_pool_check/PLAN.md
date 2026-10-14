# Saved full-pool comparison: ExaA samples 0 and 3

Declared after reciprocal solvent job4132 and before this saved-data calculation.
No new molecular evaluations, geometry changes, searches or calibration.

Question: do the two transferred-ion geometries lower either metal's energy
beyond the actual historical three-candidate pool? The earlier four-versus-two
comparison omitted adaptive_Ca and cannot answer this question.

Use both exact consumed sources, their192-atom dry oxidized-PQQ state, source
maps/caps, Ca−3/La−2 singlets, and the current MACE-vacuum plus native frozen-charge
ALPB recipe. The donor-response inventory already pins all six original cells
per source, including adaptive_Ca. Job4132 pins both transferred geometries for
both metals. Revalidate exact XYZ/state/model/parameter/runtime and accepted
component-v2/native-diagonal/vacuum-precision evidence. Historical absolute-v1
point failures remain explicit; no failed component is replaced by a different
model or by a vacuum-only score.

For each source, compare these complete common pools:

1. Historical3, tie order: origin, adaptive_Ca, adaptive_La.
2. Extended5, same first three, then transplanted_origin,
   transplanted_adaptive_La.

Both metals must have every declared cell for a pool to be available. Preserve
all valid individual cells and explicit failures. Reuse the existing frozen
candidate scorer's choose function, parameterizing only its in-memory ORDER for
the separately named extended pool; restore the original order after each call.
Exact ties retain the first candidate. Mathematical selection always chooses the
minimum; operational selection uses that minimum only when its work relative to
the original origin is strictly less than−0.1kcal/mol. Do not move the origin,
inherit the old selected geometry, or substitute the two-candidate baseline.

Energy is E_OMOL,vac(eV)×23.06054783061903 plus
G_nativeALPB(q_GFN2,vac)(hartree)×627.509474. Report per-metal selected coordinates,
energies, origin-relative MACE/solvent/total works, the Ca-minus-La scalar, and
extended5−historical3 differences separately for operational and mathematical
policies. Reproduce the actual saved historical3 result as an identity/regression
check before interpreting any gain. Reuse its existing1e−8kcal/mol algebra
tolerance. No new class calls or band transfer: classification stays null.

This is development evidence from two consumed sources. Lower selected energy
does not establish better discrimination or a correct ion position. The prior
reciprocal-response finding remains separate from incremental full-pool utility.

All numerical parsing, validation and collection run in a single one-CPU Slurm
job, submitted by root. Login work is restricted to source/metadata inspection,
pinning and code preparation. Every original/transferred cell is saved evidence;
maximum new molecular calls is zero. Output is versioned and immutable.
