# One XoxF hydrogen-placement recovery

Prepared one validated La/Ca pair for `PQQSEQ_3cad7f010e2b9d9374c4`, source gene
`PLM2_60_b1_jun17_scaffold_1437_2`, retaining selected AF3 seed101/sample0.
No endpoint was submitted by this preparation agent. The134 ready cases and
running original quantum job are untouched; the two X/UNK targets remain unscored.

Ready manifest:
`workspaces/plm_xoxf_all_fixed_core_20260916/retry_1200785/hydrogen_repair/continuation/prepared_pairs.json`

SHA256: `228ef66075e4a1db0dc74a0093ad49ba8f25e70785ffda9d9af188b8c4633e7e`.
Adjacent `target_outcomes.json` contains exactly one ready target.

The original failure was present ARG347 HD2 lying1.252523Å from CD, outside the
unchanged1.25Å carver attachment cutoff. No atom was missing. Fragment assembly
failed before the batch's post-carve H-recovery hook could execute. All failed
source files, captured System, original PDB and the original outcome remain intact.

The same captured OpenMM System was restarted using only existing H coordinates;
all heavy atoms and atom/protonation identities stayed fixed. The first recovery
used the original internal tolerance1.0 and stopped after169 iterations at
H-only force RMS1.034419kJ/mol/nm. It remains recorded **FAIL**. The parent then
authorized continuing those exact full-precision coordinates with internal
tolerance0.5, retaining the H-only acceptance threshold1.0 and original objective.
Two additional iterations reached **0.714564kJ/mol/nm**; no new H generation or
random draws were performed. The continuation allowed at most831 further steps,
so the two recovery attempts together remained below1,000 iterations.

PDB serialization gives a separately measured H force RMS **3.025895kJ/mol/nm**.
The serialized structure passes geometric checks; it is **not claimed to satisfy
the full-precision force-convergence threshold**. No calibration was rerun or refit.
Numerical preparation changes are explicit; the fixed chemistry, source heavy
geometry, native ORCA method and existing computational bands remain unchanged.

The repaired core contains80 atoms: C27H31N6O15 plus the metal, with La−2/Ca−3
charges. Protein H–parent distances are1.180805–1.194674Å, below the strict1.25Å
attachment cutoff. Minimum H–H distance is1.711714Å. All27 PQQ atoms, five generated
caps,24 protein core heavy atoms and every full-protein heavy coordinate are
unchanged. Full9026 atom identities are preserved;4428 full-protein H coordinates
changed, including23 in the core. Paired La/Ca nuclear coordinates are identical.

The unchanged original carver rebuilt the core using the retained repaired PDB;
PDBFixer was not reexecuted. Its temporary input-replay hook and new protonation
manifest explicitly record that distinction. Recovery and preparation scripts
are pinned independently from the original wrapper/helpers.

A packaging comparison initially compared full-precision generated caps with
their six-decimal metadata and correctly halted. `finalize.py` applies the exact
frozen `atom_record`/XYZ serialization to the comparison; it performs no new
minimization or preparation and changes no cap coordinates. The executed
`prepare_repaired.py` and `continuation.py` remain unchanged and pinned.

`review.py` reproduces the original actual parent-count rejection and confirms
the repaired fragment, whole-protein identities, exact PQQ/caps, strict H gate,
paired coordinates and parent endpoint prerequisites. Receipts under the ready
directory: `recovery_validation.json`, `geometry_validation.json`,
`independent_validation.json`, and `implementation_pins.json`. Both the executed
recovery script and finalizer are pinned. Summary receipt: `validation.json` here.
