# Correct actual stereocenter identity before candidate scoring

The original generic four-neighbor-carbon guard was chemically overrestrictive:
it treated CH2/CH3 atom-name permutations as stereochemistry. This is a preparation/
optimizer implementation defect, not a discovery that the protein cannot respond.
The completed OpenCL `searches_v2` attempts and their force evaluations remain
immutable and will not be scored. Root explicitly approved a uniform correction
and the same 24 starts again, before any candidate classification was inspected.

Concrete evidence from actual final coordinates:

- 4MAE's original-target search stopped at 0.32537 Å heavy movement with the
  tetrahedral volume at **LYS387 CE (CH2)** reduced to 1.000218e−5 Å³.
- Q9Z4J7's original-target search stopped at 0.25718 Å with **PRO440 CG (CH2)**
  reduced to 1.000087e−5 Å³.

Those are not stereocenters. Neither has a reason to retain a signed atom-name
volume. Generic volume blocking must not be interpreted as validated mechanical
stationarity. Source H bond lengths and identities remain retained; this fix
neither normalizes H coordinates nor changes protonation.

## Uniform correction, frozen before new forces

New protocol: `source_ff19SB_collective_matched_donor_target_proposal_v2`.

Protect every actual standard-amino-acid non-Gly CA and Ile/Thr CB. Require the
source graph to establish four substituents, exactly one H, and the correct same-
residue heavy neighbors: CA has N/C/CB; Ile CB has CA/CG1/CG2; Thr CB has CA/OG1/CG2.
Unsupported nonstandard parent residues or malformed centers fail explicitly.
Preserve each actual source signed volume, without imposing ideal L stereochemistry.
The existing 1e−5 Å³ nonplanarity guard now applies only to these actual centers.
All actual peptide cis/trans checks remain unchanged.

Nothing else about the source set, 8 Å shell, 24 starts/targets, serialized ff19SB,
source bond lengths, fixed atoms, cap reconstruction, 0.8 Å heavy limit, overlap
guards, stationarity tests, finite optimizer or admission changes. The previously
documented [OpenCL double recovery](OPENCL_RECOVERY.md) remains in effect. Record
full physical rejection details in addition to each rejection reason; this only
improves diagnostics and changes no numerical decision.

Execute `searches_v3` only after v2 is terminal. This is the root-authorized
technical retry of the intended chirality-preserving model. Count all earlier
force calls/allocations; they are not free or relabeled successful calculations.
No MACE/GFN2 score has been run on v2. A finite feasible v3 endpoint is still only
a proposal; large residual forces or boundary stops remain explicit.
