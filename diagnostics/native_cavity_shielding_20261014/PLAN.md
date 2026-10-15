# Fixed-charge protein-exclusion scout — declared before execution

Jacob reaffirmed autonomous execution: “why wouldn't you queue it”. Root owns
submission and integration. This diagnostic follows the completed, closed
midpoint probe; it does not reopen that geometry search or modify production.

## Question and energy

Does adding the actual omitted protein volume materially change the local core's
electrostatic solvent response, and how much comes from synthetic cap artifacts?
For each metal and source, report

`T_M = Gelec(capped source union, q_core ⊕ 0) − Gelec(local core, q_core)`

and `T_Ca − T_La`, then its sample3-minus-sample0 difference. These are matched
component transfers, not affinity or calibrated class scores. No second full
solvation energy, permanent protein field, direct Coulomb or mechanical term is
added. Report native HB/SASA/reference energies but do not use them in T.

Use unchanged qualified xTB6.7.1 native ALPB factory: water/P16/gsolv,298.15K,
230 angular points,zero salt,existing elemental radii/descreening/35Å Born cutoff.
Exact archived vacuum GFN2 populations, generated at their original300K SCC setting,
remain fixed. Do not normalize charges, move atoms, substitute incorrect ORCA La
parameters, alter PQQ state or tune radii. Element-specific Ca/La radii remain
native and must be recorded. The entire matched source is included, not a tuned
new radial shell. Zero environment charges are an explicit isolation of excluded
volume, not a chemical charge/protonation assignment.

## Sources and construction

Use already consumed two-Ca ExaA samples0(control) and3(ambiguous source), origin
geometries only, Ca−3/La−2 singlets. The exact local192 atoms,14 synthetic caps,
three carver PQQ H, dry oxidized-PQQ3− state, both source-metal positions and
all available standard-protein hydrogens are preserved. Append each omitted
protonated-source atom exactly once through source identities, never proximity.
The active source Ca is replaced by the intended local metal; spectator Ca remains
as a native-radius zero-charge cavity center. The first192 charge values are
unchanged and all appended values are exactly zero.

Synthetic caps remain cavity/descreening/shape centers because the native API
cannot separate charge sites from cavity atoms. Record cap-to-source contacts;
this is a capped-union diagnostic, not a source-only physical cavity. No radius
zeroing, omitted physical atoms or cap-charge redistribution is permitted.
Preparation checks paired coordinates, atom/source/charge/state correspondence,
exact duplicates, file pins and the declared cap/PQQ H inventory on Slurm.

## Cap attribution and numerical checks

Retain native electrostatic real–real, real–cap and cap–cap blocks and the ALPB
Q² shape contribution separately; their sum must recover the native scalar.
On the same union matrix, also evaluate exact zero cap charges. Compare that
algebraic vector with a cap-removed cavity retaining all real charges/coordinates.
Both controls conserve their own Q_real. This measures cap-volume/descreening
sensitivity of the real charges; it is not a physical state or an error bound
for the charged-cap calculation.

First scout: sample0 Ca. Six cavity updates: local replay, union, cap-removed,
translated union(+11,−7,+3Å), rotated union(+90° about z: x,y,z→−y,x,z), and
reversed atom order with every charge/mask/map permuted identically. No random
perturbations. Repeat unchanged energy-vector calls within the same matrix.
Native-versus-block algebra, repeat, local replay and each rigid/permutation
Gelec residual must satisfy the predeclared1e−9Eh allowance. Preserve the existing
analytic bound for eight-decimal charges; only archived uncertain charges receive
5e−9e uncertainty, and added exact zeros receive none. The bound must remain below
1e−6Eh. Replay comparisons include their recorded charge-precision bounds.
SASA/HB rotational changes are reported separately and cannot silently pass an
electrostatic failure. No numerical tolerance depends on a desired biological sign.

Only after this scout passes, execute sample0La/sample3Ca/sample3La union and
cap-removed variants: six further updates,12total. Their old local native energies
are reused. Keep all four endpoint rows, explicit missing/failure reasons, actual
counts and cap-control limitations. No automatic parameter sweep or third source.
There are zero new SCC, MACE, DFT, optimizations or folds.

## Execution and continuation

Preparation uses one shared Slurm CPU. A thin new caller links the existing pinned
native static libraries without engine changes/installation. It emits scalar
components, core-sized blocks and per-atom metadata; never a whole-protein dense
matrix text dump or redundant full matrix copy. Internal TBorn memory remains
quadratic and is estimated from actual prepared N before calls.

The qualified build is serial; requesting extra slots does not imply parallel
native computation. Start a single native fixture at a time, with one numerical
thread and sufficient scheduler-granted memory. Initial plan is four shared CPUs
(default16GB), subject to actual atom-count allocation check; increase to12 only
if justified by the explicit memory model. No memory/time overrides, exclusive
allocation, test-partition science or priority changes. Record requested resources,
actual CPU/RSS/time and the serial-engine limitation. Reuse measured scout memory
before admitting the finite continuation. Root arms exact completion wakes.

A failed mapping/numerical gate stops the dependent native continuation and remains
visible. A diagnosed implementation defect can be repaired in a new immutable
version without repeating successful molecular cells where compatible. No unseen
biological class or new score threshold is assigned. Once the finite matrix is
collected, decide whether any component change is substantial relative to cap/shape
and numerical effects before proposing a classifier alteration.

Ownership: Curie prepares actual sources; Heisenberg implements the thin native
caller/runner; Gauss independently reviews; root freezes plan, launches and reports.
