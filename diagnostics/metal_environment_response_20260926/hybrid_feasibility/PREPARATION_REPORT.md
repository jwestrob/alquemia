# Exact-source mechanics and PQQ cross-LJ preparation now exist

The authorized follow-on produced reusable standard parameters without a new
energy model or molecular calculation. Metal LJ selection and the complete
hybrid interaction ledger remain root-owned and unresolved; this is preparation
progress, not a validated force field or classifier.

## Exact current 1H4I protein

`workspaces/metal_environment_response_20260926/hybrid_preparation_v1/protein_v2/`
contains the full ff19SB `system.xml`, ordered atom charges/sigma/epsilon,
49,903 native nonbonded exceptions, and bonded-term supports with actual indices
and native parameters. The complete XML retains CMAP tables and all values.
Relative to the present 23-atom QM protein subset, crossing terms are 3 bonds,
18 angles, 59 periodic torsions, and zero CMAP terms. These are inventory counts,
not a prescription to add every term to the electronic Hamiltonian.

The exporter calls unchanged `source_protein()` and merely retains the System
it creates. All 9,113 source protein identities match; total charge is −8 to
floating precision. Max coordinate difference from archived physical source:
7.11e−15 Å. Applying the existing exact boundary charge ledger reconstructs all
9,087 current exterior charges with **zero difference**. No source H/terminal
repair or atom addition was performed. Actual export took 12.46 seconds, no
energy/force evaluation.

One initial observer wrapper failed OpenMM's signature-aware argument validation
because it hid the original method signature. `functools.wraps` fixed only that
Python observation behavior; source function and physical arguments stayed
unchanged. The failed `protein/FAILED_ATTEMPT.json` remains beside successful
`protein_v2`.

## PQQ: real graph, standard atom types, LJ only

The primary [RCSB PQQ CCD](https://files.rcsb.org/ligands/download/PQQ.cif)
was downloaded and content-pinned. It describes neutral oxidized PQQ C14H6N2O8.
The current prepared state is already oxidized PQQ3− C14H3N2O8. Its conversion
is explicit: remove the CCD acid protons HOB2/HOB7/HOB9, assign formal −1 to
their O2B/O7B/O9B parents, and leave heavy bond orders/oxidation state unchanged.
All 24 heavy names match actual source PQQ; three existing prepared hydrogens
map uniquely to HN1/H3/H8 and their existing parent atoms. Nothing is optimized.

Installed antechamber ran with GAFF2 and `-j 1` (atom types only), `-nc -3`,
`-m 1`, fixed sequence, and **no charge-calculation option**. Its log shows only
the atomtype executable; no SQM, DFT, AM1-BCC, geometry optimization or molecular
force/energy call ran. All 29 explicit source/CCD-derived bond orders survive.
Nine assigned types (`c ca cc cd ha hn na nb o`) all have direct entries in the
pinned GAFF2 2.2.20 LJ table; no missing parameter was guessed or zero-filled.

The returned scientific artifact is
`.../pqq_types_v1/PQQ_CROSS_LJ.json`: 27 mapped source atoms with exact archived
coordinates, formal state, standard atom types, Rmin/2 and epsilon. Its
`partial_charge_e` fields are **null**, and energy/mixing/metal models are null.
The raw typing-only mol2 has format-default zero charge columns; it is explicitly
not a physical PQQ charge model and must not enter classical Coulomb evaluation.

Antechamber's intermediate AC format rounds coordinates to three decimals, after
MDL's four decimals. An initial extraction check against a single rounding stage
failed. Inspection established both serializations; the final check reproduces
that exact two-stage transformation at 1e−12 Å and verifies atom order and the
complete bond-order graph. The LJ-only artifact then restores exact source
coordinates. No scientific displacement/geometry threshold was loosened and no
energy was evaluated on rounded intermediates.

## Files, checks and next decision

`PREPARATION_RESULT.json` pins the successful outputs and typing receipt.
`EXPORT_PROTEIN.py`, `TYPE_PQQ.py`, and `EXTRACT_PQQ_LJ.py` are the actual bounded
preparation code; installed executables/environment were unchanged. The typing
receipt records command, executable hash, return code and elapsed time.

Executed checks: whole protein source identity/coordinates, exact reconstructed
field charges, force inventory, PQQ heavy/H mapping, explicit proton/redox state,
all graph bonds, typed ordering/serialization and complete LJ lookup. These are
preparation checks, not molecular validation. No scheduler or GPU allocation.

Next: root can select a defensible metal LJ family independently of biological
scores, then define the actual additive boundary terms, exclusions and treatment
of link-dependent electronic contributions. Do not substitute the raw parent
protein System or typing-only PQQ mol2 for that ledger. No full hybrid score,
new default or production promotion has been created.
