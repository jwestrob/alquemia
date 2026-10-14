# Fixed-pool all-heavy packing preparation

Question: does explicitly representing packing omitted by the existing polar-only
context change the metal-dependent cost of the **same saved donor motions**?
The two declared sources are ExaA Q9Z4J7 two-Ca samples3 (target) and0 (first
same-conditioning correct control). These are consumed development sources.

Keep every fragment in the original canonical-La4.3Å union. Identify real donor
heavy atoms displaced by >1e-8Å between origin and either saved proposal. This
is the existing saved-motion mapping tolerance, excluding floating-point replay
noise; it is not a physical movement threshold fitted to results. Across all
three geometries, select every full-source heavy atom within4.3Å of any such
atom. Take the union of the resulting **complete fragment identities across
both sources**, then use that same fixed membership for both metals and all
three candidates. No iterative expansion and no Leu-specific exception.

Translate contacts using the existing source-graph conventions: sidechain atom
selects its complete sidechain, C/O selects the complete following peptide, N
selects the bonded preceding peptide. A CA contact selects both adjacent peptide
groups plus its sidechain (GLY has none), permitting normal graph closure through
CA. Proline ring selection includes its bonded preceding peptide. Unsupported
nearby species/topology makes preparation unavailable. Existing water inventory
is retained without adding waters. Spectator Ca stays in the full source and
must remain >4.3Å from every expanded retained heavy atom.

Reuse `consistent_context.forced_expansion`, `SourceGraph.materialize` and the
existing PQQ/cysteine/peptide/polar-sidechain handling. A private helper adds
neutral ALA/VAL/LEU/ILE/PHE/MET sidechains, validating complete heavy membership
and per-atom hydrogen counts against the actual pinned ff19SB residue template
and source H-parent graph. No source hydrogen is generated or repositioned.
No protein force-field energy or new metal parameter is introduced.

Prepare origins from the archived protonated sources. Graft each existing real
proposal through its exact physical source map; preserve all shared real atoms,
cofactor atoms and surviving caps exactly. New caps follow retained/omitted source
atoms using the established sigma-cap construction; caps have no independent
coordinates. Newly included real atoms remain at their source positions unless
already part of the actual saved physical motion. Preserve old outputs and
write new candidate/state/map identities. No optimization or new folding.

Checks before scoring: replay each old origin within1e-12Å; reproduce saved
kinematics within1e-8Å; preserve source identities/shared coordinates exactly;
retain PQQ/protons/waters; preserve source covalent bonds within1e-8Å across
candidates; reject the existing <0.45Å overlap failure; verify paired Ca/La
coordinates/charges/spin/electron parity and constant source-pool composition;
record explicit fragment charge ledgers and matching state identities across
both sources. Invalid input leaves the corresponding source unavailable.

Output is two-source preparation and a declared12-cell uncomputed matrix.
The matching12 old cells are traceable reuses, not repeat evaluations. Root owns
MACE/solvent execution and interpretation. No calibration refit, band transfer
claim, production change or affinity inference belongs to this preparation.
All scientific imports, preparation, checks and numerical work run on Slurm.
