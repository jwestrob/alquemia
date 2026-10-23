# Frozen 26-cell local-correlation method-sensitivity design

Declared before any new molecular call. Root owns releases and submissions.
All results and failed attempts are immutable. This is consumed-data development.

## Exact finite matrix

18 TightPNO controls: the six canonical_water_cp_v2 central cells and the twelve
`cc__` donor_cp_v1 cells. Retain every original coordinate, atom order, ghost
mask, charge and basis/ECP. Acetate means the archived bidentate geometry.

8 cluster cells: aquo8_central_La/Dy and aquo9_central_La/Dy, each at TightPNO
and the tighter tier below, full real complexes only. Reuse the original four
PBE0 results. No EnGrad, force, geometry optimization, CP cluster fragment,
additional basis, state scan or repeated canonical calculation is included.

## Electronic settings

RHF DLPNO-CCSD(T1), exact archived La ECP46/Dy ECP55 and 47-function metal basis;
def2-TZVP light orbitals; effective singlets with physical La1/Dy6 retained as
metadata. NoFrozenCore: correlate all 18 water-dimer/40 donor-complex electrons,
8 ion-fragment electrons, 10 water/32 donor-fragment electrons, and 88/98 W8/W9
electrons. Ghost metal has DelECP and its complete unchanged orbital/AuxC basis.

VeryTightSCF, Guess PModel, MaxIter300, NoAutostart; conventional RHF (no new
RIJCOSX/RIJK approximation). CC MaxIter100, STol1e−8, PrintLevel3. Iterative
triples TriTolE1e−8 Eh; no automatic iterative-triples substitution. The latter
is a prospective numerical choice, not the ORCA default 1e−4. No D4 on CC.

AuxC: def2-TZVP/C for H/C/N/O; exact existing metal AutoAux expansion with only
`NewAuxJGTO` renamed `NewAuxCGTO` (La351, Dy369 spherical functions). Orbital and
ECP coefficients remain byte-identical. PrintBasis retains actual expansions.
No new auxiliary generation, contraction, deletion or linear-dependence recovery.

| Tier | TCutPairs | TCutPNO | TCutMKN | TCutDO |
|---|---:|---:|---:|---:|
| tight | 1e−5 | 1e−7 | 1e−3 | 5e−3 |
| tighter | 1e−6 | 1e−8 | 1e−4 | 5e−3 |

Both start from TightPNO; use its full iterative local-MP2 guess, default AHFB
localization without randomization, and unchanged default TNO scaling. Retain
all native local/semicore/PNO/TNO cutoffs and pair/domain counts. The three
explicitly tightened cutoffs are a declared sensitivity experiment. No
extrapolation or domain-convergence certificate is claimed.

## Native and approximation admission

Each cell: normal termination/exit0; finite unique final CC energy and exact
CCSD+triples closure within1e−8 Eh; declared RHF/state/charge/electrons/NoFrozenCore;
exact orbital counts and real/ghost nuclear/ECP assignments; complete centers
and coordinates within the inherited1e−6 Å print tolerance; complete positive
AuxC counts, including metal351/369; declared CC residual≤1e−8 and actual
iterative triples terminal |Delta-E|≤1e−8. Require active local-MP2 convergence
and localization completion; preserve printed SCF residual flags, warnings,
T1 and available amplitude diagnostics. A missing native field is unavailable,
not an inferred pass. Native parser changes, if needed after the first output,
must be versioned saved-output repairs; they authorize no molecular repeat.

All18 controls must pass native admission. Compare the RHF energy to its
canonical counterpart within1e−6 Eh. Require all six per-metal CP interaction
errors, all three Dy−La CP interaction errors, and the two donor-versus-water
paired exchange errors to be ≤0.25 kcal/mol against canonical. Retain individual
total-energy errors and all denominators. This working resolution is inherited
from CORRELATED_QUALIFICATION_PLAN; it is not a physical error bar. Across each
three-role parent, require identical actual AuxC counts/contractions per center.

Only then admit the tight W9Dy cluster scout. After root's actual
cost/memory/native review, release the other seven cluster cells. Both quartet
contrasts require all four native results. Require |C_tighter−C_tight|≤0.25
kcal/mol; report all four endpoint shifts too. Failure means local sensitivity
is unresolved; do not pick the tier closer to a model.

Compare each raw CC quartet C against the corresponding original raw PBE0
electronic C and original PBE0+D4 C. Keep the fixed D4 contribution separate;
never add D4 to CC. No local energies are mixed with old PBE0 fragments into CP/B.

## Execution and finite release

26 maximum one-attempt cells; four MPI ranks, one thread each, standard_expand,
default scheduler memory/time. MaxCore=floor(0.75×actual granted MiB/4).
Root releases `controls_scout` (Dy water complex/two CP roles and La acetamide
complex,4), `controls_remaining` (14 after actual four-scout native review),
`cluster_scout` (tight W9Dy,1), then `cluster_remaining` (7). The initial four
cover real/ghost metal, fragment electron counts and the largest small control.
The aggregate control set remains18, with no repeat.
Each subset binds exact panel/source/review/preparation and unique execution
paths. Cluster stages also bind the full passed18-control qualification; the
last stage binds the admitted actual W9Dy result and root cost/memory review.
No stage uses a pending result as evidence. Missing/failed slots remain visible.

Preparation, collection and review make zero native calls. Original targets,
physical-state qualification and force/free-energy qualification remain unchanged.
