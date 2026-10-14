# 1FLG A/B: conditional matched-source preparation

Declared before any new energies. Prepare both deposited crystal chains using
the existing Q9Z4J7 dry oxidized-PQQ3−/singlet, exact canonical 4.3 Å fragment
union. No optimization, scoring, DFT, SCC, folding or calibration is part of this
job. Both chains and all failures remain in the denominator.

The archived 1FLG audit establishes 100% identity of crystal residues 1–582 to
Q9 positions 35–616. Use that complete saved map, not a guessed numbering offset.
The full deposited dimer and all waters remain immutable source evidence. For
each chain, retain its full protein, PQQ and both native Ca atoms, exclude waters
under the existing dry policy, and explicitly map identifiers to Q9's identifiers.
No source coordinate, element, bond or occupancy is edited. Reuse the frozen
crystal alternate-conformer selector (occupancy sum, A then lexical tie-break),
with unresolved PQQ/metal conformers rejected. Keep an atomwise map back to the
actual deposited identities. Record every excluded residue/conformer. The exact
already-audited CCD↔AF3 PQQ atom-name map is reused and its element/bond/donor
equivalence checked; only names change, never coordinates or chemistry.

Declared raw sites: A PQQ701/Ca751/Ca761; B PQQ702/Ca752/Ca762. The first Ca is
catalytic, the second structural. Verify the existing unique strict <4 Å PQQ
association independently; do not assume these site labels establish admission.
Before protonation, inventory source occupancies, B factors, PQQ connections,
waters, nonstandard residues, alternate conformers and selected heavy atoms.
Any covalent PQQ external connection or extra atom is unsupported. Existing
PQQ schema/connectivity and typed donor gates remain unchanged. A deposited PQQ
label is not a measured redox state; no bond is shortened or atom removed to
match an oxidized interpretation.

Primary source record/abstract supports the mature ExaA dimer and structural Ca.
Full article retrieval currently returns HTTP403; its detailed redox discussion
has not been read. Earlier saved analysis shows long C4/O4 and C5/O5 bonds and
nonplanar PQQ. These are not used to infer redox or invent a new admission cutoff.
The allowed experiment is explicitly conditional on the existing oxidized
PQQ3− model; it is not an experimentally verified matched oxidation-state test.

Reuse the unchanged spectator protonator (pH7/RNG20260914/CPU1/no nonstandard H),
`pqq_fast_prepare.prepare_case`, `motion_envelope.prepare_task` and the existing
spectator distance check (>4.3 Å to every retained physical heavy atom). Missing
heavy atoms, unsupported graph/caps, or local state mismatch remain failures.
Require exact Q9 local state signature, paired maps/charges and source-heavy
coordinate preservation. The full protein identity must equal the exact mapped
mature subset of Q9; its different terminal construct is recorded, not hidden by
copying Q9's full-protein hash. Report nearby omitted chain/water atoms using the
existing 4.3 Å radius as a descriptive inventory, not a new score-tuned gate.

Finite execution: one actual Slurm CPU, two sequential source preparations;
standard partition, one node, no memory/time overrides. The usual protonator may
internally place hydrogens, but no new heavy-atom repair or geometry optimization
is authorized. Reuse completed AF3 reference and old1FLG metrics; do not redo them.
