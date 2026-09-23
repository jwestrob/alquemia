# Collective matched-target scaffold proposals — approved finite execution

Jacob approved the native restart and collective scaffold round. Root delegated
this branch and recorded the agreement in
[the round plan](../scaffold_restart_round_20260922/PLAN.md). The earlier
[proposal](PROPOSED.md) defines the scientific comparison. This plan freezes
implementation choices before any molecular energy or force call.

## Inputs and physical comparison

Use exactly the eight already consumed sources in
`diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json`. For each, independently
start from the original, actual Ca-adaptive and actual La-adaptive geometry.
Fix the corresponding actual protein donor positions throughout that search.
These are three matched targets per source, at most 24 searches, each with one
start. Failed parent/state/mapping/initial constraints remain unavailable.

Use the same role-defined donor policy for every source: anchor Glu OE1/OE2,
anchor Asn OD1, and extra-acid Asp OD1/OD2 when that role is actual Asp. Both
carboxylate oxygens are retained regardless of contact distance. Source atom
identities and actual distances are recorded; no score or distance cutoff selects
a target oxygen. Catalytic Asp and cation-partner roles remain mobile scaffold
unless they are the same explicit target atom. Root resolved this input-schema
choice before forces; each source uses identical donor IDs for all three targets.

Water-basins' independent parent inventory pins the exact source protein atoms,
ff19SB System, bonds/proton inventory, donor identities and archived target
coordinates. Reuse the two crystals' already documented distant terminal OXT
completions only; add no new chemistry or force-field parameters.

Choose mobile residues from original source coordinates: any complete protein
residue with a heavy atom within 8 Å of a canonical-core heavy atom, every
core-contributing residue, and one layer of directly bonded peptide neighbors
of that initial residue set. Actual bonds define neighbors. Fix all exterior
protein atoms and the target donors exactly. Fix metal, PQQ, generated PQQ H,
waters and all other nonprotein inventory. Context membership, source protonation,
charges, multiplicities and cap recipe remain unchanged.

## Frozen solver and guards

The proposal objective is the exact protein-only ff19SB parent energy:
bond, angle, periodic torsion, CMAP and NoCutoff Coulomb/LJ with native exceptions.
No metal/PQQ force-field coupling, solvent, artificial spring, entropy or Hessian
term is added. Starting affine forces are evaluated; the source is not assumed
stationary. This energy generates proposals and **never contributes to the
composite score**.

Use sparse projected L-BFGS with eight stored secant pairs. Eliminate fixed
coordinates, project gradients/directions into the tangent space of all original
bond-length constraints using LSMR, and restore nonlinear bond lengths with
Newton/LSMR retraction. No dense whole-protein Hessian is constructed. Trial
directions have at most 0.05 Å atom displacement before retraction. Armijo
coefficient is 1e−4, at most 24 halvings per accepted step; retain all failed trial
reasons. Maximum accepted iterations: 200. This is a finite algorithm, not an
allocation time cap or invitation to add restarts.

Numerical settings: LSMR atol/btol 1e−11, at most 3000 iterations; retraction aims
for maximum bond error 1e−8 Å, with at most 30 Newton steps. A projected vector
must satisfy the tangent equations within 1e−7 × max(1, its Euclidean norm).
Final bond error must be ≤1e−6 Å. Fixed-atom error must be ≤1e−12 Å.

Every accepted geometry must retain source carbon stereochemistry and peptide
cis/trans signs. Source four-neighbor carbon signed volumes must exceed 1e−5 Å³;
source peptide |cos omega| must exceed 0.1, otherwise preparation is unsupported.
All heavy atoms remain within 0.8 Å of the original source (numerical tolerance
1e−7 Å). No new severe nonbonded overlap is allowed: 1.0 Å for heavy pairs or
0.55 Å when H is involved, including contacts to fixed metal/PQQ. Existing short
source contacts may not get shorter below these guards by >1e−7 Å. Exclude actual
one- and two-bond parent neighbors only. No overlap guard is a substitute for
the omitted cofactor force.

Use the installed OpenMM CUDA platform in **double precision**, deterministic
forces, device 0 in the allocated GPU namespace. Record actual platform/version/
properties. No silent backend fallback. Allocate one H200, 32 CPUs and 200000 MiB
using existing scheduler conventions; searches run sequentially. Measure the
first declared search and all subsequent actual force calls; no separate
throughput or fictitious protein experiment.

## Admission, interchange and interpretation

Stationarity requires projected maximum atom gradient ≤0.1 kcal/mol/Å and RMS
Cartesian-component gradient ≤0.03 kcal/mol/Å. Preserve actual termination and
boundary status. A finite, feasible final accepted endpoint can still be scored
after a limit or failed line search if its parent energy is no higher than its
own starting value by >1e−6 kcal/mol. Unchanged origins remain explicit; downstream
coordinate deduplication may reuse them. Do not call a nonstationary proposal a
minimum. All three targets remain required for the primary source comparison.

Each candidate pins its parent inventory, source atom order, initial and final
full-parent coordinates in Å, fixed donor target, execution receipt and physical
checks. Reconstruct exact existing context coordinates from moved physical source
atoms and the unchanged retained/omitted-bond cap map; caps never move independently.
Ca/La context files differ only at the metal identity, with unchanged endpoint
charge/multiplicity. Collective positions are not falsely expressed as old
four-angle coordinates.

Root separately owns cross-scoring of admitted geometries for both metals with
the unchanged compact composite and old candidates: at most 48 new MACE and
96 GFN2 calls. This branch launches no MACE, GFN2 or DFT jobs. Report matched target
effects and q0 response separately, then pooled known-class margins and fold
spread. These consumed sources provide development evidence, not new biological
validation; no threshold refit, default change or force-field score correction.
