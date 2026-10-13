# Two-Ca source preparation: bounded adapter plan

**Read-only design; not implemented or executed.** Root owns the refold scout and
will decide whether to implement after its actual Q9Z4J7 samples are inspected.
Heisenberg owns refold inputs/execution. No structure was prepared in this audit.

## Existing behavior and concrete gaps

1. `accommodation_folds.adapt_ca` already aliases only an explicitly selected
   Ca residue as LA, preserving its native Ca element, atom name, chain/residue
   identity and coordinates. It retains every other atom. Reuse this function;
   LA is a historical residue alias, not a La-conditioned full-source structure.
2. `normalize_af3_cif.normalize` recognizes actual PQQ chemistry and preserves
   other residues. A second AF3 Ca may retain a generic ligand residue name;
   normalize only that verified monatomic spectator residue to CA with an
   explicit old-to-new identity record. Do not delete or move either ion.
3. The pinned `protonate_standard_only.py` hardcodes `{LA:1,PQQ:1}` both in input
   admission and its output checks. The zero-H comparison is separately literal.
   Changing only the accepted input count cannot make this path work correctly.
4. `second_shell_context.discover` filters elements to N/O/S before checking
   nearby unsupported chemistry. Consequently a Ca spectator could be silently
   skipped. Existing carver second-metal rejection covers its buffered core,
   but does not establish safe exclusion from an expanded context.

## Minimal implementation after the scout

Create one research-only `prepare_two_ca_source.py`; preserve production and
all pinned helpers. It accepts explicit source/config/selected-site/spectator/
PQQ/role/union records and a fresh workspace, with a real Slurm guard before
scientific imports. No general protonation or multisite framework is needed.

- Require one structural model, exactly two actual Ca atoms and one schema-valid
  PQQ. Use the already-declared rule: exactly one Ca is within 4.0 Å of PQQ.
  Preserve both complete selectors (chain, residue, insertion code, atom,
  alternate conformer and element); ambiguity or absent association is unavailable.
- Reuse the active-Ca residue alias and PQQ normalizer; add only the spectator CA
  residue alias. Verify source atom identity/element/count and heavy coordinates
  after each step. CIF aliasing changes no coordinates; PDB serialization retains
  the existing 0.001 Å tolerance. Reject missing/added/repaired heavy atoms.
- Add a narrow new standard-only protonation wrapper around the pinned
  `protonate_cif.protonate`, copying the existing small PDBFixer subclass policy
  rather than mutating globals or making the old wrapper pretend to pass.
  Its allowed nonstandard inventory is `{LA:1,CA:1,PQQ:1}`, with both metal atoms
  still actual Ca and zero added H on all three residue types. Retain pH7,
  RNG20260914, OpenMM CPU/one thread, forcefield=None, no nonstandard CCD download,
  and no missing-residue construction. New wrapper/subprotocol IDs and pins make
  the changed occupancy explicit. Verify both ions and the complete heavy-atom map
  again after protonation. PQQ hydrogens remain the pinned carver's responsibility.
- Reuse `pqq_fast_prepare.build_core`, the existing homologous role mapping,
  typed 3.1 Å donor gate (at least six, at most two N), PQQ C14H3N2O8/3− policy,
  dry inventory, and exact original canonical La-source 4.3 Å fragment union.
  Only the selected active-site endpoint is substituted Ca↔La. The spectator
  remains Ca in the full normalized/protonated source and is never used as the
  active-site surrogate. The local endpoint charge contains only modeled atoms.
- Add an explicit spectator admission audit before accepting the expanded
  context. Record distance to active metal, PQQ, each donor anchor, and the nearest
  retained physical heavy atom across the complete core plus fixed union.
  **Spectator within or at 4.3 Å of any retained physical heavy atom is unsupported**
  in this local policy. Keep the existing core rejection too. This is a boundary
  gate, not merely spectator-to-PQQ distance. Verify the spectator is absent from
  source-to-local mappings/caps/selected fragments, present unchanged in full
  source, and never silently omitted by the N/O/S discover filter. Record explicitly
  that no spectator energy/charge is represented in the local Hamiltonian.
- Reuse existing paired-coordinate, electron/charge, source-graph/cap and fixed
  union state checks. Any new occupancy source that fails them stays unavailable.
  Do not widen selection, change water/protons, or alter bands in response.

## Finite qualification, on Slurm only

First use the already-requested real 1FLG chain/site mapping as a preparation
fixture, preserving its chosen chain and both real Ca atoms. No invented molecule
or dummy values. The crystal can check spectator retention/selector/protonation
contracts; it is not a new numerical reference or a classifier test by itself.
Then admit every Q9 scout sample under the same frozen rule, including failures.
Subsequent C5AXV8/Q88JH5 preparation is conditional on root's scout decision.

Required checks: full-source two-Ca identity and heavy-map preservation; no
nonstandard H; unique active selector; fixed core/context pairing and state;
actual spectator boundary distances; explicit rejection of ambiguous/nearby
spectators; unchanged existing single-Ca preparation regression on a pinned real
fixture. Corrupted copies, if needed for negative parser checks, must be marked
as such and never reported as observed structures.

New source IDs remain `conditioned_2Ca_v1`, with the old source ID retained as a
separate provenance field. Two-Ca occupancy is directly motivated for Q9/ExaA;
it remains a tested structural hypothesis for the other two proteins. They do
not replace the original one-Ca benchmark samples or create independent proteins.
A successful preparation establishes neither better classification nor a global
spectator/metal response model. No energy evaluation is part of this adapter plan.
