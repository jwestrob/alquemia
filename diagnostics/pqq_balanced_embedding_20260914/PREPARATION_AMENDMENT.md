# Pre-energy preparation amendment: corrupted C5 terminal oxygens

**Freeze status:** this amendment was written before generating point-charge
files and before calculating or viewing any embedded electronic energy. It
does not change the frozen panel, boundary rule, endpoints, or decision
thresholds in `EXPERIMENT.md`.

## Defect found during input audit

The protonated C5 source models contain one geometrically corrupted terminal
oxygen on each donor-bearing chain used by the experiment. OpenMM accepts the
topologies and assigns ff19SB charges, so an uncorrected whole-chain embedding
would silently place a terminal oxygen point charge tens of angstroms from its
carbonyl carbon.

| model and target chain | corrupted atom | C--O distance | sound partner |
|---|---|---:|---|
| C5B120 monomer, chain A | ASN 601 OXT | 35.3438 A | O, 1.2500 A from C |
| C5B120 dimer conformation, chain B | ASN 601 O | 51.4108 A | OXT, 1.2790 A from C |

The defect is present upstream in the predicted structure coordinates; it was
not introduced by this embedding experiment. It did not affect the local QM
carves because ASN 601 is outside their QM regions. The source files remain
immutable. Only isolated embedding-coordinate copies will be repaired and
hashed.

The analogous corrupted OXT on non-target dimer chain A is recorded but is
irrelevant because the frozen environment contains target chain B only.

## Frozen deterministic repair

For each affected terminal ASN, preserve C, CA, and the geometrically sound
terminal oxygen. Let `u` be the unit vector from C to CA and `v` the unit
vector from C to the sound oxygen. Decompose

`v_parallel = dot(v,u) u` and `v_perpendicular = v - v_parallel`.

The repaired oxygen direction is the reflection

`v_repaired = v_parallel - v_perpendicular`,

and its coordinate is

`r_repaired = r_C + |r_sound-r_C| v_repaired`.

Thus the repaired atom has exactly the sound C--O bond length and the
opposite in-plane orientation about the C--CA axis. No other coordinate,
atom identity, protonation state, force-field charge, or QM coordinate is
changed. The repair is identical for the 3.3- and 3.6-A environments and for
La and Ca.

## Mandatory pre-energy checks

- Both terminal C--O distances must be 1.15--1.40 A after repair.
- The repaired atom and complete repaired-chain coordinate files must be
  SHA-256 recorded.
- The ff19SB chain charges must remain -3 e for both C5 chains.
- La and Ca at a fixed model/radius must receive byte-identical point charges.
- The original charge-closure, distance, and boundary-ledger checks in
  `EXPERIMENT.md` remain mandatory.

Any failure stops preparation. No embedded energy may be interpreted from an
input that did not pass these checks.
