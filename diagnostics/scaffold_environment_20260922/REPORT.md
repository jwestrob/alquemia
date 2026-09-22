# Scaffold feasibility: protein parameters available, additive score unresolved

**No-go for an additive scaffold score in the present implementation.** All three
real standard-protein parents can be parameterized, but the matched capped-local
reference needed to avoid overlapping the MACE boundary response is unavailable.
No energy, force, Hessian, optimization or new molecular job was run. No classifier,
default, H preparation or biological label changed.

This closes a bounded feasibility question, not the broader possibility of useful
protein mechanics. A hard displacement limit does not identify which missing
physics matters. Metal translation is being investigated independently.

## What is concretely available

The pinned ff19SB XML and OpenMM 8.5.1 instantiate these standard-protein systems:

| Source | Parent atoms | Protein atoms in MACE context | Exterior atoms | Real cut bonds / caps | Parent charge, e |
|---|---:|---:|---:|---:|---:|
| 1H4I | 9,060 | 114 | 8,946 | 12 | −9 |
| 4MAE | 8,826 | 163 | 8,663 | 11 | −3 |
| PLM8344 | 8,678 | 149 | 8,529 | 13 | −1 |

The crystals require the **already archived** terminal OXT completions at LYS595
and GLU577, respectively. Their coordinates and exact source/FF hashes match the
older whole-protein preparation. Those OXT atoms lie 54.95 and 42.46 Å from the
site. No new terminal choice is made. The original sources alone fail terminal
template matching; that failure is retained in each inventory. PLM8344 matches
without completion. These conditional parents have no assigned metal/PQQ forces.

Every original protein atom, including every H, keeps its source coordinates.
The old whole-protein radial H normalization is explicitly **not** reused.
Source-to-physical mapping differences are at most 1.43e−14 Å from floating-point
unit conversion; the frozen mapping tolerance is 1e−12 Å. Paired Ca/La physical
maps are identical. Every cap replaces exactly one real crossing protein bond.
The full per-atom mapping, charges, source proton/template identities, every
bonded-term support, nonbonded exceptions and parameterized System XML are saved
under `workspaces/scaffold_environment_20260922/inventory_v2/` and hash-pinned in
[RESULT.json](RESULT.json).

All three MACE contexts contain their 27-atom PQQ and selected metal, plus the
protein fragments and caps: 154, 202 and 190 total context atoms. Each source has
24 PQQ heavy atoms; its three carver-generated PQQ H remain local-only and fixed.
Source PQQ and metal are explicitly inventoried, not silently discarded as known
force-field chemistry. There are no other source heterogens or explicit waters
in these three prepared systems. Actual protein histidine templates and disulfide
states are retained, with no reassignment: HID/HIE counts 8/2, 5/4 and 2/2; each
parent contains four CYX residues. PLM biological metal preference is unknown.

## Exact parent energy and interaction partition

The available *protein-only* parent is the ordinary fixed-charge ff19SB function

`U_P = U_bond + U_angle + U_periodic_torsion + U_CMAP + U_nonbonded`.

The serialized native parameters define bond/angle harmonics, periodic proper and
improper torsions, CMAP tables, Coulomb/Lennard-Jones interactions and their exact
exclusions/1–4 exceptions. Nonbonded method is `NoCutoff`, constraints are absent,
and no solvent or metal/cofactor force is included. `CMMotionRemover` is recorded
but is not an energy term. No temperature-dependent free-energy claim is made.

Let C be the actual protein atoms included in the local MACE context, and E the
remaining parent protein atoms. Every native term is classified by its entire
atom support: wholly C, wholly E, or mixed C/E. For nonbonded terms, the inventory
retains all exceptions and the exact complementary ordinary-pair sets, rather
than writing tens of millions of redundant pair records.

| Mixed term count | 1H4I | 4MAE | PLM8344 |
|---|---:|---:|---:|
| Bonds | 12 | 11 | 13 |
| Angles | 66 | 60 | 72 |
| Periodic torsions, including improper terms | 149 | 156 | 173 |
| CMAP terms | 6 | 6 | 6 |
| Explicit nonbonded exceptions, including zero exclusions | 242 | 227 | 260 |
| Ordinary nonbonded C/E pairs | 1,019,602 | 1,411,842 | 1,270,561 |

This partition is exact for the specified **parent force field**. It does not
prove that every mixed term is absent from the local learned energy. MACE already
sees the capped fragment's covalent geometry and nonbonded interactions. Each cap
tracks the real retained/omitted bond direction through a physical chain-rule
map. Its length is fixed, so not every boundary response is duplicated; nonetheless
its angular and many-body response overlaps some proposed boundary mechanics.

Consequently, `U_P − U_terms_wholly_C` is calculable but is not established as
“only omitted scaffold energy.” Simply adding it to MACE could count boundary
stiffness twice. A matched protein subtraction would instead require a defined
low-level energy of the **same capped local protein** with the same cap map and
compatible parameters. The actual local capped topology has 10, 9 and 11 residues
that ff19SB cannot match. No cap charges, atom types, equilibrium terms or dummy
cofactor parameters were invented to complete that subtraction.

There is a second boundary to a full relaxation model: exterior protein–PQQ/metal
interactions. The standard-protein XML does not define them. Fixed exterior and
fixed metal/PQQ geometry can make some omitted fixed-charge interactions constant
within a conditional candidate-work comparison; this does not justify omitting
their forces when the exterior or metal moves, nor ignoring density response.
No claim of a complete hybrid parent follows from protein charge closure.

The compact scorer already uses `E_OMOL + E_GFN2,ALPB − E_GFN2,vac`. No CPCM/ALPB
or full-protein solvation term has been added here. A future electrostatic coupling
also needs a consistent subtraction/reference; the protein FF inventory is not
a license to add a second solvent energy or FF charges to high-level atoms.

## Why a raw parent force would be ambiguous

The source H-coordinate convention is far from ff19SB bond equilibrium. Of
4,439 / 4,331 / 4,271 H-containing bonds, **3,647 / 3,694 / 3,550** differ by more
than 0.1 Å; median absolute differences are 0.1073 / 0.1079 / 0.1069 Å. These are
coordinate/parameter comparisons, not computed forces or evidence of biological
strain. The corresponding local-context counts are 39/52, 62/74 and 50/68.
Some radial components disappear when projected onto rigid bond-preserving
torsions, so these findings do not prove that all projected mechanics are invalid.
They do rule out interpreting a large unprojected parent force as a clean missing
scaffold signal, or assuming the starting structure is a stationary basin.

## Conditional next test — proposed, not executed

If a scaffold diagnostic remains useful after the joint metal/angular result,
the bounded test is **nine protein-only parent evaluations**: q0 and the two actual
archived four-angular candidate geometries for each of these three sources.
Use the completed `proposals_v1/after_proposals_1209857.json`, not a new search.
Transfer mapped physical source coordinates only; keep exterior/proton inventory,
OXT records and all native parameters fixed. Evaluate the same FF on one CPU
worker, retaining per-term energy/forces and work relative to q0. Project through
the current physical Jacobian and retain all local/mixed/exterior components.

Purpose: determine whether the independent protein model reports reproducible
boundary/exterior work along already sampled motions. It is **not** a classifier
correction, accuracy test, minimum search, or permission to ignore the cap overlap.
The exact runtime has not been measured. No allocation or executor is prepared;
root review is required before starting this optional diagnostic. No new model is
needed merely to demonstrate the existing coupling blocker.

For any later elimination of moving exterior coordinates, retain the nonzero
starting forces and constant as well as curvature:

`K_eff = K_aa − K_ab solve(K_bb,K_ba)`

`g_eff = g_a − K_ab solve(K_bb,g_b)`

`c_eff = c − 0.5 g_bᵀ solve(K_bb,g_b)`.

The inverse is only on a physically constrained stable subspace. No rigid-mode
clamping, assumed zero exterior gradient, invented springs, whole dense Hessian,
or entropy term is justified by this inventory.

## Execution and checks

Successful inventory: 21.416 s wall, 21.464 process CPU-s; zero molecular energy
or force calls. An earlier 13.225 s preparation attempt failed while assembling
the *local template probe* because OpenMM requires atoms in each residue to be
contiguous. The technical fix orders existing cap atoms within their residues;
all failed records remain in `inventory_v1`. It changes no source chemistry.

Five actual-fixture tests pass: source/proton preservation, paired context/cap
mapping, complete parameter/term coverage, a corrupted real map missing one cap,
and an explicitly corrupted PQQ identity. The first test run used bitwise equality
across a unit conversion and failed at 3.55e−15 Å; it now uses the inventory's
already frozen 1e−12 Å criterion. Both logs are retained. No scientific integration
test was replaced by fabricated output. [Commands](COMMANDS.md) reproduce the
inventory with explicit inputs and a new output directory.

**Recommendation:** keep the scaffold score unavailable. The concrete missing
piece is a validated boundary/reference model, not more parent compute. Continue
the separately declared accommodation tests without adding an arbitrary scaffold
penalty or attributing their constraints to this one cause.
